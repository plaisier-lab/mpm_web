import json
import numpy as np
from database import dbconn
from flask import Blueprint, render_template, request

jacks_page = Blueprint("jacks_page", __name__, template_folder="templates")

def get_jacks_data(bicluster_id):
	db = dbconn()
	c = db.cursor()

	c.execute(
		"""SELECT g.symbol, score, std, p_value, cl.name
		FROM gene_jacks_result gjr
		JOIN gene g ON gjr.gene_id = g.id
		JOIN cell_line cl ON gjr.cell_line_id = cl.id
		WHERE g.id IN (
			SELECT bg.gene_id FROM bicluster b
			JOIN bic_gene bg ON b.id = bg.bicluster_id
			WHERE b.id = %s
		);""",
		[bicluster_id]
	)
	data = c.fetchall()

	# order by cell line
	output = {}
	for symbol, score, std, p_value, cell_line_name in data:
		if cell_line_name not in output:
			output[cell_line_name] = []
		output[cell_line_name].append((symbol, score, std, p_value))

	db.close()

	return output

def get_jacks_grna_data(gene_name, cell_line_name):
	db = dbconn()
	c = db.cursor()

	c.execute(
		"""SELECT grna.name, gs.sequence, gjr.mean, gjr.p_value
    FROM grna_jacks_result gjr
    JOIN grna ON grna.id = gjr.grna_id
    JOIN gene g ON g.id = grna.gene_id
    JOIN grna_sequence gs ON gs.grna_id = grna.id
    JOIN cell_line cl ON cl.id = gjr.cell_line_id
    WHERE g.symbol = %s AND cl.name = %s;""",
		[gene_name, cell_line_name]
	)
	data = c.fetchall()

	db.close()

	return data

@jacks_page.route('/grna/<gene>/<cell_line>')
def gnra_jacks_data(gene=None, cell_line=None):
	return json.dumps(get_jacks_grna_data(gene, cell_line))

def get_meso1_data(bicluster_id):
	db = dbconn()
	c = db.cursor()

	c.execute(
		"""SELECT g.symbol, gjr.mean, gjr.p_value
    FROM grna_jacks_result gjr
    JOIN grna gr ON gr.id = gjr.grna_id
    JOIN gene g ON g.id = gr.gene_id
    JOIN cell_line cl ON cl.id = gjr.cell_line_id
    WHERE cl.name="MESO1" AND g.id IN (
      SELECT bg.gene_id
      FROM bic_gene bg
      JOIN bicluster b ON bg.bicluster_id = b.id
      WHERE b.id = %s
    )
    AND gjr.p_value = (
      SELECT MIN(gjr2.p_value)
      FROM grna_jacks_result gjr2
      JOIN grna gr2 ON gr2.id = gjr2.grna_id
      JOIN gene g2 ON g2.id = gr2.gene_id
      WHERE g2.id = g.id
      AND gjr2.cell_line_id = cl.id
    );""",
		[bicluster_id]
	)
	genes = c.fetchall()

	db.close()

	return genes

def get_achilles_data(gene):
	db = dbconn()
	c = db.cursor()
	
	c.execute(
		"""SELECT score
    FROM achilles_results ar
    JOIN gene g ON ar.gene_id = g.id
    WHERE g.symbol = %s;""",
		[gene]
	)
	results = [score[0] for score in c.fetchall()]

	db.close()

	return {
		"name": gene,
		"low": np.min(results),
		"q1":  np.percentile(results, q=25.0),
		"median": np.median(results),
		"q3": np.percentile(results, q=75.0),
		"high": np.max(results),
	}

@jacks_page.route('/achilles/<gene>/')
def achilles_data(gene=None):
	return json.dumps(get_achilles_data(gene))

def get_achilles_data_from_bicluster(bicluster):
	db = dbconn()
	c = db.cursor()

	c.execute(
		"""SELECT ar.score, g.symbol
    FROM achilles_results ar
    JOIN gene g ON ar.gene_id = g.id
    WHERE g.id IN (
			SELECT g2.id
			FROM bic_gene bg
			JOIN gene g2 ON bg.gene_id = g2.id
			JOIN bicluster b ON bg.bicluster_id = b.id
			WHERE b.name = %s
		);""",
		[bicluster]
	)
	results = c.fetchall()
	output_dictionary = {}
	for score, gene in results:
		if gene not in output_dictionary:
			output_dictionary[gene] = []
		output_dictionary[gene].append(score)
	
	output = []
	for gene in output_dictionary.keys():
		output.append({
			"name": gene,
			"low": np.min(output_dictionary[gene]),
			"q1":  np.percentile(output_dictionary[gene], q=25.0),
			"median": np.median(output_dictionary[gene]),
			"q3": np.percentile(output_dictionary[gene], q=75.0),
			"high": np.max(output_dictionary[gene]),
		})

	db.close()

	return output

@jacks_page.route('/achilles-bicluster/<bicluster>/')
def achilles_data_bicluster(bicluster=None):
	return json.dumps(get_achilles_data_from_bicluster(bicluster))