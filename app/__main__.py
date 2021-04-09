from __future__ import print_function

##########################################################
## OncoMerge:  app.py                                   ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @github: https://github.com/plaisier-lab/mpm_web     ##
## @Author:  Chris Plaisier                             ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

import logging
import json
import csv

from flask import Flask, Response, render_template, request

from constants import HALLMARKS, SELECTABLE_PHENOTYPES, SELECTABLE_PHENOTYPES_BLACKLIST
from database import dbconn
from bicluster import bicluster_page
from causal_analysis import causal_analysis_page
from search import search_page, get_index_locals
from jacks import jacks_page

app = Flask(__name__)
# app.config.from_envvar('MESO_SETTINGS')
app.register_blueprint(bicluster_page)
app.register_blueprint(causal_analysis_page)
app.register_blueprint(search_page)
app.register_blueprint(jacks_page)


@app.errorhandler(Exception)
def unhandled_exception(e):
	app.logger.exception(e)
	return render_template('unknown_error.html')

@app.route('/')
def index():
	hallmarks, selectable_phenotypes = get_index_locals()
	return render_template('index.html', hallmarks=hallmarks, selectable_phenotypes=selectable_phenotypes)

@app.route('/network')
def network():
	return render_template('network.html')

@app.route('/about')
def about():
	return render_template('about.html')

@app.route('/download')
def download():
	return render_template('download.html')

@app.route('/citation')
def citation():
	return render_template('citation.html')

@app.route('/genecompletions')
def genecompletions():
	term = request.args.get('term')
	db = dbconn()
	try:
		c = db.cursor()
		c.execute("""SELECT symbol FROM gene WHERE symbol LIKE %s""", [str(term)+'%'])
		tmpGene = [i[0] for i in c.fetchall()]
		c.execute("""SELECT name FROM mirna WHERE name LIKE %s""", [str(term)+'%'])
		tmpMiRNA = [i[0] for i in c.fetchall()]
		json1 = json.dumps(tmpGene+tmpMiRNA)
	finally:
		db.close()
	return Response(response=json1, status=200, mimetype='application/json')

@app.route('/combinatorial_network')
def combinatorial_network():
	with open(app.config['NODES_FILE'], 'r') as infile:
		csvreader = csv.reader(infile, delimiter=',')
		csvreader.next()
		nodes = {node_id: {'id': node_id, 'tf_ko': tf_ko, 'in_gbm': in_gbm}
				 for node_id, tf_ko, in_gbm in csvreader}

	with open(app.config['EDGES_FILE'], 'r') as infile:
		csvreader = csv.reader(infile, delimiter=',')
		csvreader.next()
		edges = []
		for edge, sig_coocc in csvreader:
			source, edge_type, target = edge.split()
			edges.append({'source': source, 'target': target, 'type': edge_type,
						  'sig_coocc': sig_coocc})

	graph_data = []
	for node_id, node_data in nodes.items():
		classes = []
		if node_id.startswith('hsa-miR'):
			classes.append('mirna')
		else:
			classes.append('gene')

		if node_data['tf_ko'] == 'Yes':
			classes.append('crispr')
		if node_data['in_gbm'] == 'Yes':
			classes.append('in_gbm')
		if 'in_gbm' in classes and 'crispr' in classes:
			classes.append('crispr_gbm')
		graph_data.append({ 'data': { 'id': node_id }, 'classes': ' '.join(classes) })

	for i, edge in enumerate(edges):
		if edge['sig_coocc'] == 'Yes':
			graph_data.append({ 'data': { 'id': 'e%d' % i, 'source': edge['source'], 'target': edge['target'] }, 'classes': 'sigcoocc' })
		else:
			graph_data.append({ 'data': { 'id': 'e%d' % i, 'source': edge['source'], 'target': edge['target'] } })

	return render_template('combinatorial_network.html', **locals())


if __name__ == '__main__':
	handler = logging.StreamHandler()
	handler.setLevel(logging.INFO)
	app.debug = True
	app.secret_key = 'supercalifragilistic'
	app.logger.addHandler(handler)
	app.run(host='0.0.0.0', debug=True)
