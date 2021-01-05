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

#!/usr/bin/env python
import sys
import traceback as tb
import logging
import json
import os
from os import path, walk
import csv
from functools import wraps

# import MySQLdb
from flask import Flask, Response, url_for, redirect, render_template, request, session, flash, jsonify
import flask

import math
import itertools
import gzip
import pandas
import numpy as np
import base64
from sklearn.linear_model import LinearRegression
# import rpy2.robjects as robjects
from scipy import stats
from colour import Color

from constants import NUM_PARTS, HALLMARK_TO_ICON, HALLMARKS, ENRICHMENT_PHENOTYPES, SEX_PHENOTYPES, PHENOTYPES_DICT, GRAPH_COLOR_MAP, GRAPH_COLOR_GRADIENTS, SELECTABLE_PHENOTYPES, SELECTABLE_PHENOTYPES_BLACKLIST, USE_PHENOTYPE_SCATTERPLOT, PHENOTYPE_INDEX_TO_UINAME
from database import dbconn
from bicluster import bicluster_page
from causal_analysis import causal_analysis_page

app = Flask(__name__)
# app.config.from_envvar('MESO_SETTINGS')
app.register_blueprint(bicluster_page)
app.register_blueprint(causal_analysis_page)
# General helpers


# def phyper(q, m, n, k, lower_tail=False):
#    """calls the R function phyper"""
#    r_phyper = robjects.r['phyper']
#    kwargs = {'lower.tail': lower_tail}
#    return float(r_phyper(float(q), float(m), float(n), float(k), **kwargs)[0])

def phyper(x, M, n, bigN, lower_tail=False):
    """uses scipy.stats to compute hypergeometric overlap p-value"""
    return 1 - stats.hypergeom.cdf(float(x), float(M), float(n), float(bigN))

# rewritten function, works to the best of my knowledge
def cluster_data(cursor, cluster_id):
    cursor.execute(
        """SELECT value, p.name FROM gene_expression ge
        JOIN gene g ON ge.gene_id=g.id
        JOIN patient p ON ge.patient_id=p.id
        WHERE g.entrez IN (
            SELECT g.entrez from bic_gene bg
            JOIN gene g ON bg.gene_id=g.id WHERE bicluster_id=%s
        );""",
        [cluster_id]
    )
    results = cursor.fetchall()
    grouped_by_patient = {}
    for value, patient in results:
        if patient not in grouped_by_patient:
            grouped_by_patient[patient] = []    
        grouped_by_patient[patient].append(value)
    
    data = []
    for patient in grouped_by_patient:
        gene_expressions = grouped_by_patient[patient]
        data.append([
            patient,
            np.min(gene_expressions), # minimum
            np.percentile(gene_expressions, q=25.0), # lower quartile
            np.median(gene_expressions), # median
            np.percentile(gene_expressions, q=75.0), # upper quartile
            np.max(gene_expressions), # maximum
        ])
    return sorted(data, key=lambda x: x[3])


def subtype_enrichment(cursor, cluster_id, phenotype_name, phenotype_min_max):
    # now make a phenotype map
    cursor.execute("""SELECT p.name, pd.phenotype_string, pd.phenotype_value FROM pheno_data pd
        JOIN patient p ON pd.patient_id=p.id
        JOIN phenotype pt ON pd.phenotype_id=pt.id
        WHERE pt.name=%s;""",
        [phenotype_name]
    )
    data = cursor.fetchall()

    string_ptmap = {patient: phenotype for patient, phenotype, _ in data} # don't change
    value_ptmap = {patient: value for patient, _, value in data} # don't change
    all_patients = {patient for patient in string_ptmap.keys()} # don't change

    if phenotype_min_max == None:
        phenotypes = {phenotype for phenotype in string_ptmap.values()} # don't change
    else: # build range list if we're dealing with continuous values
        phenotypes = []
        step = (phenotype_min_max[1] - phenotype_min_max[0]) / 5
        for i in range(0, 5):
            phenotypes.append((step * i + phenotype_min_max[0], step * (i + 1) + phenotype_min_max[0]))

    # sorted by median pValue
    sorted_patient_names = [row[0] for row in cluster_data(cursor, cluster_id)]

    # group patients into phenotype groups.
    # NOTE: the ptmap items need to be sorted, otherwise groupby fails to group correctly
    if phenotype_min_max == None:
        pt_patients = itertools.groupby(sorted(string_ptmap.items(), key=lambda pair: pair[1]), key=lambda pair: pair[1])
        pt_patients = {phenotype: set(map(lambda p: p[0], patients)) for phenotype, patients in pt_patients}
    else: # build the phenotype quintiles
        pt_patients = {}
        for min_max in phenotypes:
            pt_patients[min_max] = []
            for patient in value_ptmap.keys():
                value = value_ptmap[patient]
                if value >= min_max[0] and (value < min_max[1] or (value == phenotype_min_max[1] and phenotype_min_max[1] == min_max[1])):
                    pt_patients[min_max].append(patient)

    num_columns = len(sorted_patient_names)
    cols_per_part = int(math.floor(num_columns / NUM_PARTS))
    # print "# columns: %d # cols/part: %d" % (num_columns, cols_per_part)

    pvalues = []
    min_pvalue = 100.0
    max_pvalue = -100.0

    for i in range(NUM_PARTS):
        part_pvalues = {}
        start = cols_per_part * i
        end = (cols_per_part * (i + 1)) - 1

        # adjust end for the last part
        if i == (NUM_PARTS - 1) and end != (num_columns - 1):
            end = num_columns - 1
        cur_patients = [p_i for p_i in sorted_patient_names[start:end + 1]]
        # print "Part %d, %d-%d, # current patients: %d" % (i, start, end, len(sorted_patient_names))
        for phenotype in phenotypes:
            x = len([p for p in cur_patients if p in pt_patients[phenotype]])
            bigN = len(cur_patients)
            n = len(pt_patients[phenotype])
            n = len([p for p in sorted_patient_names if p in pt_patients[phenotype]])
            M = len(sorted_patient_names)
            pvalue = phyper(x, M, n, bigN)
            
            if x != 0 and pvalue != 0.0:
                if pvalue == 0.0:
                    pvalue = 10e-10

                if pvalue <= 0.5:
                    pvalue = -math.log10(2 * pvalue)
                else:
                    pvalue = math.log10(2 * (1.0 - pvalue))

                if math.isinf(pvalue):
                    signum = -1 if pvalue < 0 else 1
                    pvalue = signum * -math.log10(10e-10)

                if pvalue < min_pvalue:
                    min_pvalue = pvalue
                if pvalue > max_pvalue:
                    max_pvalue = pvalue

            part_pvalues[phenotype] = pvalue
        pvalues.append(part_pvalues)

    return pvalues, min_pvalue, max_pvalue

def get_phenotype_min_max(cursor, phenotype_name, bueno_deep_filter=False):
    if bueno_deep_filter == False:
        cursor.execute("""
            SELECT phenotype_value FROM pheno_data pd
                JOIN patient p ON pd.patient_id=p.id
                JOIN phenotype pt ON pd.phenotype_id=pt.id
                WHERE pt.name=%s;
        """, [phenotype_name])
    else:
        cursor.execute("""
            SELECT phenotype_value FROM pheno_data pd
                JOIN patient p ON pd.patient_id=p.id
                JOIN phenotype pt ON pd.phenotype_id=pt.id
                WHERE pt.name=%s
                AND p.id IN (SELECT DISTINCT patient_id FROM bueno_deep_filter);
        """, [phenotype_name])
    data = cursor.fetchall()
    phenotype_data = [item[0] for item in data]

    if phenotype_data[0] == None:
        return None
    else:
        return (np.min(phenotype_data), np.max(phenotype_data))

def get_phenotype_color(phenotype_name, value, min_max):
    min, max = min_max
    color1 = Color(GRAPH_COLOR_GRADIENTS[phenotype_name][0])
    color2 = Color(GRAPH_COLOR_GRADIENTS[phenotype_name][1])
    percent = (value - min) / (max - min)
    hex = Color(rgb=(color1.red * (1 - percent) + color2.red * percent, color1.green * (1 - percent) + color2.green * percent, color1.blue * (1 - percent) + color2.blue * percent)).hex

    if len(hex) == 4:
        return hex + hex[1] + hex[1] + hex[1]
    else:
        return hex

######################################################################
# Available application paths
######################################################################

@app.errorhandler(Exception)
def unhandled_exception(e):
    app.logger.exception(e)
    return render_template('unknown_error.html')

def get_index_locals():
    selectable_phenotypes = [('None', '')]
    for name, value in SELECTABLE_PHENOTYPES:
        if value not in SELECTABLE_PHENOTYPES_BLACKLIST:
            selectable_phenotypes.append((name, value))
    
    hallmarks = [('None', '')]
    for hallmark in HALLMARKS:
        hallmarks.append((hallmark, hallmark))
        
    return [hallmarks, selectable_phenotypes]

@app.route('/')
def index():
    hallmarks, selectable_phenotypes = get_index_locals()
    return render_template('index.html', hallmarks=hallmarks, selectable_phenotypes=selectable_phenotypes)


@app.route('/search')
def search():
    gene = request.args.get('gene')
    phenotype = request.args.get('phenotype')
    hallmark = request.args.get('hallmark')

    db = dbconn()
    c = db.cursor()
    bme_type ='gene'
    if not gene and phenotype == "None" and hallmark == "None":
        hallmarks, selectable_phenotypes = get_index_locals()
        return render_template('index.html', hallmarks = hallmarks, selectable_phenotypes = selectable_phenotypes, invalid_gene = True)

    if gene.find('hsa-')==-1:
        c.execute("""SELECT symbol FROM gene WHERE symbol=%s""", [gene])
        geneData = c.fetchall()
    else:
        c.execute("""SELECT name FROM mirna WHERE name=%s""", [gene])
        geneData = c.fetchall()
        bme_type = 'mirna'
    db.close()

    if len(geneData)==0 and phenotype == "None" and hallmark == "None":
        hallmarks, selectable_phenotypes = get_index_locals()
        return render_template('index.html', hallmarks = hallmarks, selectable_phenotypes = selectable_phenotypes, invalid_gene = True)

    query_statements = []
    if phenotype:
        query_statements.append("phenotype=%s" % phenotype)
    
    if hallmark:
        query_statements.append("hallmark=%s" % hallmark)
    
    query_statement = ""
    if len(query_statements):
        query_statement = "?%s" % "&".join(query_statements)

    symbol = geneData[0][0]
    if bme_type == 'gene':
        return redirect(url_for('gene', symbol=symbol) + query_statement)
    else:
        return redirect(url_for('mirna', symbol=symbol) + query_statement)


def __get_muts(c, gene_pk, symbol, hallmark_filter=None):
    # Get causal flows downstream of mutation in gene
    muts = {}
    muts['name'] = symbol
    muts['flows'] = 0
    muts['regs'] = []
    muts['tfs'] = []
    muts['miRNAs'] = []
    muts['biclusters'] = []
    muts['data'] = []

    c.execute("""SELECT * FROM somatic_mutation WHERE locus_id IS NULL AND ext_id=%s""", [gene_pk]) # discrep. (we never have a mutation type of gene). test: make sure that we only get genes (not loci) from this query
    tmp_muts = c.fetchall()
    if len(tmp_muts)==1:
        c.execute("""SELECT * FROM causal_flow WHERE somatic_mutation_id=%s""", [tmp_muts[0][0]])
        tmp_cf = c.fetchall()
        for cf1 in tmp_cf:
            g1 = ''
            if cf1[3]=='tf':
                # c.execute("""SELECT * FROM tf_regulator WHERE tf_id=%s AND bicluster_id=%s""", [cf1[2], cf1[4]])
                c.execute("""SELECT * FROM tf_regulator WHERE id=%s""", [cf1[2]])
                tf = c.fetchall()
                if len(tf) > 0:
                    c.execute("""SELECT symbol FROM gene WHERE id=%s""", [tf[0][2]])
                    g1 = c.fetchall()[0][0]
                    if not g1 in muts['regs']:
                        muts['regs'].append(g1)
                        muts['tfs'].append(g1)
            else:
                c.execute("""SELECT * FROM mirna_regulator WHERE id=%s""", [cf1[2], cf1[4]])
                mirna = c.fetchall()
                if len(mirna) > 0:
                    c.execute("""SELECT name FROM mirna WHERE id=%s""", [mirna[0][2]])
                    g1 = c.fetchall()[0][0]
                    if not g1 in muts['regs']:
                        muts['regs'].append(g1)
                        muts['miRNAs'].append(g1)
            
            # print the bicluster in the table
            if not g1=='':
                c.execute("""SELECT name, survival, survival_p_value FROM bicluster WHERE id=%s""", [cf1[4]])
                b1 = c.fetchall()[0]
                c.execute(
                    "SELECT hm.name FROM hallmark hm join bic_hal bh on hm.id=bh.hallmark_id  WHERE bh.bicluster_id=%s",
                    [cf1[4]]
                )
                tmp1 = [i[0] for i in c.fetchall()]
                
                has_hallmark = True
                if hallmark_filter != None and len(hallmark_filter) > 0:
                    has_hallmark = False
                    
                    # make sure tmp1 has at least one element in hallmark_filter
                    for name in tmp1:
                        index = HALLMARKS.index(name) + 1
                        if index in hallmark_filter:
                            has_hallmark = True
                            break

                if has_hallmark:
                    h1 = list(set([convert[i] for i in tmp1]))
                    h2 = [(i,convert[i]) for i in tmp1]
                    muts['data'].append([symbol, g1, b1[0], b1[1], b1[2], h2])

                    if not b1 in muts['biclusters']:
                        muts['biclusters'].append(b1[0])
    return muts


def __get_regulators(c, symbol, bme_type, hallmark_filter=None):
    regs = {}
    regs['name'] = symbol
    regs['biclusters'] = 0
    regs['data'] = []
    tmp_regs = c.fetchall()
    if len(tmp_regs) > 0:
        # Collect all biclusters downstream regulated by TF or miRNA
        for reg in tmp_regs:
            action = 'Rep.'
            if bme_type =='gene' and reg[3] == 'activator':
                action = 'Act.'
            c.execute("""SELECT name, survival, survival_p_value FROM bicluster WHERE id=%s""", [reg[1]])
            b1 = c.fetchall()[0]
            c.execute("SELECT hm.name FROM hallmark hm join bic_hal bh on hm.id=bh.hallmark_id WHERE bh.bicluster_id=%s",
                      [reg[1]])
            tmp1 = [i[0] for i in c.fetchall()]

            has_hallmark = True
            if hallmark_filter != None and len(hallmark_filter) > 0:
                has_hallmark = False
                
                # make sure tmp1 has at least one element in hallmark_filter
                for name in tmp1:
                    index = HALLMARKS.index(name) + 1
                    if index in hallmark_filter:
                        has_hallmark = True
                        break

            # filter hallmarks using advanced search parameters
            if has_hallmark:
                h1 = list(set([convert[i] for i in tmp1]))
                h2 = [(i,convert[i]) for i in tmp1]
                regs['data'].append([symbol, action, b1[0], b1[1], b1[2], h2])
                regs['biclusters'] = regs['biclusters'] + 1
    return regs


@app.route('/mirna')
@app.route('/mirna/<symbol>')
def mirna(symbol=None, defaults={'symbol': None}):
    # Get biclusters regulated by mirna
    db = dbconn()
    c = db.cursor()
    c.execute("""SELECT id FROM mirna WHERE name=%s""", [symbol])
    mirna_pk = c.fetchone()[0]
    muts = __get_muts(c, mirna_pk, symbol)
    c.execute("""SELECT * FROM mirna_regulator WHERE mirna_id=%s""", [mirna_pk])
    regs = __get_regulators(c, symbol, "mirna")
    db.close()
    return render_template('search.html', gene=symbol, muts={}, regs=regs, bics={})


@app.route('/gene')
@app.route('/gene/<symbol>')
def gene(symbol=None, defaults={'symbol': None}):
    phenotype = request.args.get('phenotype')
    phenotype = phenotype if phenotype else None
    hallmark = request.args.get('hallmark')
    hallmark = hallmark if hallmark else None

    # decode the hallmark
    hallmarks = []
    if hallmark:
        hallmarks = json.loads(base64.b64decode(hallmark))

        # verify hallmarks
        if type(hallmarks) != list:
            hallmarks = []
        
        verified = True
        for number in hallmarks:
            if type(number) != int:
                verified = False

        if not verified:
            hallmarks = []
    
    db = dbconn()
    c = db.cursor()
    c.execute("""SELECT id FROM gene WHERE symbol=%s""", [symbol])
    gene_pk = c.fetchone()[0]
    muts = __get_muts(c, gene_pk, symbol, hallmark_filter=hallmarks)

    c.execute("""SELECT * FROM tf_regulator WHERE tf_id=%s""", [gene_pk])
    regs = __get_regulators(c, symbol, "gene", hallmark_filter=hallmarks)

    # Get biclusters that gene resides
    bics = {}
    bics['name'] = gene
    bics['biclusters'] = 0
    bics['data'] = []
    c.execute(
        """SELECT * FROM bic_gene bg
        JOIN bicluster b ON bg.bicluster_id=b.id
        WHERE bg.gene_id=%s""",
        [gene_pk]
    )
    tmp_bics = c.fetchall()
    if len(tmp_bics) > 0:
        for bic1 in tmp_bics:
            c.execute("SELECT hm.name FROM bic_hal bh JOIN hallmark hm on bh.hallmark_id=hm.id WHERE bh.bicluster_id=%s",
                      [bic1[3]])
            tmp1 = [i[0] for i in c.fetchall()]

            # filter hallmarks using advanced search parameters
            if hallmark == None or hallmark in tmp1:
                h1 = list(set([convert[i] for i in tmp1]))
                h2 = [(i,convert[i]) for i in tmp1]
                bics['data'].append([bic1[4], bic1[5], bic1[6], bic1[7], bic1[8], h2])
                bics['biclusters'] = bics['biclusters'] + 1
    db.close()
    return render_template('search.html', gene=symbol, muts=muts, regs=regs, bics=bics)


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
