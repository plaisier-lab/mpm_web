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
import csv
from functools import wraps

# import MySQLdb
import mysql.connector
from flask import Flask, Response, url_for, redirect, render_template, request, session, flash, jsonify
import flask

import math
import itertools
import gzip
import pandas
import numpy as np
from sklearn.linear_model import LinearRegression
# import rpy2.robjects as robjects
from scipy import stats
from colour import Color


NUM_PARTS = 5

convert = {'Evading apoptosis': 'cellDeath.gif', 'Evading immune detection': 'avoidImmuneDestruction.gif', 'Genome instability and mutation': 'genomicInstability.gif', 'Insensitivity to antigrowth signals': 'evadeGrowthSuppressors.gif', 'Limitless replicative potential': 'immortality.gif',
    'Reprogramming energy metabolism': 'cellularEnergetics.gif', 'Self sufficiency in growth signals': 'sustainedProliferativeSignalling.gif', 'Sustained angiogenesis': 'angiogenesis.gif', 'Tissue invasion and metastasis': 'invasion.gif', 'Tumor promoting inflammation': 'promotingInflammation.gif'}

print(os.getcwd())

app = Flask(__name__)
# app.config.from_envvar('MESO_SETTINGS')

######################################################################
# General helpers
######################################################################


def dbconn():
    '''
    'user': app.config['USER'],
    'password': app.config['PASS'],
    'host': app.config['HOST'],
    'port': app.config['port'],
    'database': app.config['DB'],
    '''

    '''
    'user': "root",
    'password': "root",
    'host': "db",
    'port': "3306",
    'database': "mpm_web",
    '''

    config = {
        'user': "root",
        'password': "root",
        'host': "db",
        'port': "3306",
        'database': "mpm_web",
    }
    return mysql.connector.connect(**config)


def read_exps():
    with gzip.open(app.config['GENE_EXPR_FILE'], 'rb') as f:
        return pandas.read_csv(f, sep=',', index_col=0, header=0)

######################################################################
# Graph/Visualization functionality
######################################################################


'''GRAPH_COLOR_MAP = {
    'control': '#6abd45',
    'classical': 'black',
    'neural': '#32689b',
    'NA': 'grey',
    'g_cimp': '#8a171a',
    'proneural': '#ed2024',
    'mesenchymal': '#faa41a'
}'''

# order in which the enrichment phenotypes are ordered
'''ENRICHMENT_PHENOTYPES = [
    'g_cimp', 'proneural', 'neural', 'classical', 'mesenchymal', 'control'
]'''

ENRICHMENT_PHENOTYPES = [
    "Sarcomatoid", "Epithelioid", "Desmoplastic", "Biphasic",
]

SEX_PHENOTYPES = ['M', 'F']

PHENOTYPES_DICT = {
    "sex": SEX_PHENOTYPES,
    "histology_WHO": ENRICHMENT_PHENOTYPES,
}

GRAPH_COLOR_MAP = {
    'Sarcomatoid': '#2c6dd4',
    'Epithelioid': '#8a171a',
    'Desmoplastic': '#ed2024',
    'Biphasic': '#faa41a',
    'M': '#6abd45',
    'F': '#2c6dd4',
    'NA': 'grey',
}

GRAPH_COLOR_GRADIENTS = {
    'age_at_surgery': ('#333333', '#2cd46d'),
    'survival': ('#333333', '#2cd46d'),
    'preop_treatment': ('#333333', '#2cd46d'),
    'B.cells.naive': ('#333333', '#2cd46d'),
    'B.cells.memory': ('#333333', '#2cd46d'),
    'Plasma.cells': ('#333333', '#2cd46d'),
    'T.cells.CD8': ('#333333', '#2cd46d'),
    'T.cells.CD4.naive': ('#333333', '#2cd46d'),
    'T.cells.CD4.memory.resting': ('#333333', '#2cd46d'),
    'T.cells.CD4.memory.activated': ('#333333', '#2cd46d'),
    'T.cells.follicular.helper': ('#333333', '#2cd46d'),
    'T.cells.regulatory..Tregs.': ('#333333', '#2cd46d'),
    'T.cells.gamma.delta': ('#333333', '#2cd46d'),
    'NK.cells.resting': ('#333333', '#2cd46d'),
    'NK.cells.activated': ('#333333', '#2cd46d'),
    'Monocytes': ('#333333', '#2cd46d'),
    'Macrophages.M0': ('#333333', '#2cd46d'),
    'Macrophages.M1': ('#333333', '#2cd46d'),
    'Macrophages.M2': ('#333333', '#2cd46d'),
    'Dendritic.cells.resting': ('#333333', '#2cd46d'),
    'Dendritic.cells.activated': ('#333333', '#2cd46d'),
    'Mast.cells.resting': ('#333333', '#2cd46d'),
    'Mast.cells.activated': ('#333333', '#2cd46d'),
    'Eosinophils': ('#333333', '#2cd46d'),
    'Neutrophils': ('#333333', '#2cd46d'),
}

SELECTABLE_PHENOTYPES = [
    ("WHO Histology", "histology_WHO"),
    ("Patient Sex", "sex"),
    ("Age At Surgery", "age_at_surgery"),
    ("Survival", "survival"),
    ("Preop Treatment", "preop_treatment"),
    ("B-Cells Naive", "B.cells.naive"),
    ("B-Cells Memory", "B.cells.memory"),
    ("Plasma Cells", "Plasma.cells"),
    ("T-Cells CD8", "T.cells.CD8"),
    ("T-Cells CD4 Naive", "T.cells.CD4.naive"),
    ("T-Cells CD4 Memory Resting", "T.cells.CD4.memory.resting"),
    ("T-Cells CD4 Memory Activated", "T.cells.CD4.memory.activated"),
    ("T-Cells Follicular Helper", "T.cells.follicular.helper"),
    ("T-Cells Regulatory TRegs", "T.cells.regulatory..Tregs."),
    ("T-Cells Gamma Delta", "T.cells.gamma.delta"),
    ("NK-Cells Resting", "NK.cells.resting"),
    ("NK-Cells Activated", "NK.cells.activated"),
    ("Monocytes", "Monocytes"),
    ("Macrophages M0", "Macrophages.M0"),
    ("Macrophages M1", "Macrophages.M1"),
    ("Macrophages M2", "Macrophages.M2"),
    ("Dendritic Cells Resting", "Dendritic.cells.resting"),
    ("Dendritic Cells Activated", "Dendritic.cells.activated"),
    ("Mast Cells Resting", "Mast.cells.resting"),
    ("Mast Cells Activated", "Mast.cells.activated"),
    ("Eosinophils", "Eosinophils"),
    ("Neutrophils", "Neutrophils"),
]

USE_PHENOTYPE_SCATTERPLOT = True

PHENOTYPE_INDEX_TO_UINAME = {item[1]: item[0] for item in SELECTABLE_PHENOTYPES}

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

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/bicluster-causal-analysis/<mutation>/<regulator>/<bicluster>/<phenotype_name>')
# mutation and regulator are gene symbol names, bicluster is the bicluster name
def bicluster_causal_analysis(mutation=None, regulator=None, bicluster=None, phenotype_name="histology_WHO"):
    db = dbconn()
    c = db.cursor()

    """
    mutation gene expression
    """

    js_mutation_gene_expression_data = []
    js_mutation_gene_expression_colors = ["#AAAAAA", "#ff2b2b"]
    js_mutation_gene_expression_name = regulator + " Expression"

    c.execute("""
        SELECT ge.value, bd.value, pt.name FROM gene_expression ge
            JOIN patient pt on ge.patient_id = pt.id
            JOIN gene g ON ge.gene_id = g.id
            JOIN tf_regulator tf ON tf.tf_id = g.id
            JOIN causal_flow cf ON cf.regulator_id = tf.id AND cf.regulator_type = "tf"
            JOIN bicluster b ON cf.bicluster_id = b.id
            JOIN bueno_deep_filter bd ON bd.somatic_mutation_id = cf.somatic_mutation_id AND bd.patient_id = ge.patient_id
            JOIN somatic_mutation sm ON cf.somatic_mutation_id = sm.id
            LEFT JOIN gene g2 ON sm.ext_id = g2.id
            LEFT JOIN locus l on sm.locus_id = l.id
            WHERE (g2.symbol=%s OR l.locus_name=%s) AND g.symbol=%s AND b.name=%s;
    """, [mutation, mutation, regulator, bicluster])
    gene_expression_data = c.fetchall()
    wt = [item[0] for item in gene_expression_data if item[1] == 1]
    mutated = [item[0] for item in gene_expression_data if item[1] == 0]

    # do statistics on WT
    js_mutation_gene_expression_data.append({
        "name": "WT",
        "low": np.min(wt),
        "q1":  np.percentile(wt, q=25.0),
        "median": np.median(wt),
        "q3": np.percentile(wt, q=75.0),
        "high": np.max(wt),
        "fillColor": js_mutation_gene_expression_colors[0],
    })

    # do statistics on mutated
    js_mutation_gene_expression_data.append({
        "name": "Mutated",
        "low": np.min(mutated),
        "q1":  np.percentile(mutated, q=25.0),
        "median": np.median(mutated),
        "q3": np.percentile(mutated, q=75.0),
        "high": np.max(mutated),
        "fillColor": js_mutation_gene_expression_colors[1],
    })

    # do statistics on the data we have
    js_mutation_gene_expression_stats = stats.ttest_ind(mutated, wt)

    """
    eigengene expression
    """

    js_bicluster_eigengene_expression_data = []
    js_bicluster_eigengene_expression_colors = ["#AAAAAA", "#ff2b2b"]
    js_bicluster_eigengene_expression_name = bicluster + " Expression"

    c.execute("""
        SELECT e.value, bd.value, pt.name  FROM eigengene e
            JOIN patient pt on e.patient_id = pt.id
            JOIN bicluster b ON e.bicluster_id = b.id
            JOIN causal_flow cf ON cf.bicluster_id = b.id
            JOIN bueno_deep_filter bd ON bd.somatic_mutation_id = cf.somatic_mutation_id AND bd.patient_id = e.patient_id
            JOIN tf_regulator tf ON cf.regulator_id = tf.id AND cf.regulator_type = "tf"
            JOIN gene g ON g.id = tf.tf_id
            JOIN somatic_mutation sm ON cf.somatic_mutation_id = sm.id
            LEFT JOIN gene g2 ON sm.ext_id = g2.id
            LEFT JOIN locus l on sm.locus_id = l.id
            WHERE (g2.symbol=%s OR l.locus_name=%s) AND g.symbol=%s AND b.name=%s;
    """, [mutation, mutation, regulator, bicluster])
    eigengene_data = c.fetchall()
    wt = [item[0] for item in eigengene_data if item[1] == 1]
    mutated = [item[0] for item in eigengene_data if item[1] == 0]

    # do statistics on WT
    js_bicluster_eigengene_expression_data.append({
        "name": "WT",
        "low": np.min(wt),
        "q1":  np.percentile(wt, q=25.0),
        "median": np.median(wt),
        "q3": np.percentile(wt, q=75.0),
        "high": np.max(wt),
        "fillColor": js_bicluster_eigengene_expression_colors[0],
    })

    # do statistics on mutated
    js_bicluster_eigengene_expression_data.append({
        "name": "Mutated",
        "low": np.min(mutated),
        "q1":  np.percentile(mutated, q=25.0),
        "median": np.median(mutated),
        "q3": np.percentile(mutated, q=75.0),
        "high": np.max(mutated),
        "fillColor": js_bicluster_eigengene_expression_colors[1],
    })

    # do statistics on the data we have
    js_bicluster_eigengene_expression_stats = stats.ttest_ind(mutated, wt)

    """
    residual aka "bicluster conditioned on mutation expression"
    """

    js_bicluster_residual_data = []
    js_bicluster_residual_colors = ["#AAAAAA", "#ff2b2b"]
    js_bicluster_residual_name = bicluster + " Expression Conditioned on " + regulator + " Expression"

    # we already got the gene expression and eigengene expression data previously, so no sql fetches here
    flattened_gene_expression = [(i[0],) for i in gene_expression_data]
    flattened_eigengene = [(i[0]) for i in eigengene_data]
    regression = LinearRegression().fit(flattened_gene_expression, flattened_eigengene)
    residual = flattened_eigengene - regression.predict(flattened_gene_expression)

    wt = []
    for i in range(0, len(flattened_gene_expression)):
        if eigengene_data[i][1] == 1:
            wt.append(residual[i])
    
    mutated = []
    for i in range(0, len(flattened_gene_expression)):
        if eigengene_data[i][1] == 0:
            mutated.append(residual[i])
    
    # do statistics on WT
    js_bicluster_residual_data.append({
        "name": "WT",
        "low": np.min(wt),
        "q1":  np.percentile(wt, q=25.0),
        "median": np.median(wt),
        "q3": np.percentile(wt, q=75.0),
        "high": np.max(wt),
        "fillColor": js_bicluster_residual_colors[0],
    })

    # do statistics on mutated
    js_bicluster_residual_data.append({
        "name": "Mutated",
        "low": np.min(mutated),
        "q1":  np.percentile(mutated, q=25.0),
        "median": np.median(mutated),
        "q3": np.percentile(mutated, q=75.0),
        "high": np.max(mutated),
        "fillColor": js_bicluster_residual_colors[1],
    })

    # do statistics on the data we have
    js_bicluster_residual_stats = stats.ttest_ind(mutated, wt)

    """
    scatter plot:
    x-axis: TF gene expression
    y-axis: bicluster eigengenes
    gene_expression_data: (float, boolean, Patient)
    eigengene_data: (float, boolean, Patient)
    """

    phenotype_min_max = get_phenotype_min_max(c, phenotype_name, True)

    # intersect gene_expression_data and eigengene_data by patient
    gene_expression_data_by_patient = {}
    for datum in gene_expression_data:
        gene_expression_data_by_patient[datum[2]] = (datum[0], datum[1])
    
    eigengene_data_by_patient = {}
    for datum in eigengene_data:
        eigengene_data_by_patient[datum[2]] = (datum[0], datum[1])
    
    js_scatterplot_data = {}
    if phenotype_min_max == None:
        for phenotype in PHENOTYPES_DICT[phenotype_name]:
            js_scatterplot_data[phenotype] = {}
            
            js_scatterplot_data[phenotype][0] = {
                "name": "{} (Mut)".format(phenotype),
                "data": [],
                "color": GRAPH_COLOR_MAP[phenotype],
                "marker": {
                    "symbol": "circle",
                    "radius": 3,
                },
            }

            js_scatterplot_data[phenotype][1] = {
                "name": "{} (WT)".format(phenotype),
                "data": [],
                "color": GRAPH_COLOR_MAP[phenotype],
                "marker": {
                    "symbol": "cross",
                    "lineColor": None,
                    "lineWidth": 2,
                },
            }
    else:
        js_scatterplot_data[0] = {
            "name": "Mutated",
            "data": [],
            "color": GRAPH_COLOR_GRADIENTS[phenotype_name][1],
            "marker": {
                "symbol": "circle",
                "radius": 3,
            },
        }

        js_scatterplot_data[1] = {
            "name": "WT",
            "data": [],
            "color": GRAPH_COLOR_GRADIENTS[phenotype_name][1],
            "marker": {
                "symbol": "cross",
                "lineColor": None,
                "lineWidth": 2,
            },
        }

    use_table = gene_expression_data_by_patient
    if len(eigengene_data) < len(gene_expression_data):
        use_table = eigengene_data_by_patient
    
    c.execute("""
        SELECT p.name, pd.phenotype_string, pd.phenotype_value FROM pheno_data pd
            JOIN patient p ON pd.patient_id=p.id
            JOIN phenotype pt ON pd.phenotype_id=pt.id
            WHERE pt.name=%s;
    """, [phenotype_name])
    phenotype_data = c.fetchall()

    string_by_patient = {patient: string for patient, string, _ in phenotype_data}
    value_by_patient = {patient: value for patient, _, value in phenotype_data}
    
    all_data = []
    highest_x = -10000
    lowest_x = 10000
    for patient in use_table.keys():
        x = gene_expression_data_by_patient[patient][0]
        y = eigengene_data_by_patient[patient][0]
        boolean = eigengene_data_by_patient[patient][1]

        if phenotype_min_max == None:
            if string_by_patient[patient] not in js_scatterplot_data:
                continue
            
            js_scatterplot_data[string_by_patient[patient]][boolean]["data"].append({
                "x": x,
                "y": y,
                "name": patient,
            })
        else:
            js_scatterplot_data[boolean]["data"].append({
                "x": x,
                "y": y,
                "name": patient + " - " + PHENOTYPE_INDEX_TO_UINAME[phenotype_name] + ": " + "{:.2f}".format(value_by_patient[patient]),
                "color": get_phenotype_color(phenotype_name, value_by_patient[patient], phenotype_min_max),
            })

        all_data.append([x, y])

        if x < lowest_x:
            lowest_x = x
        
        if x > highest_x:
            highest_x = x
    
    js_scatterplot_stats = stats.pearsonr([i[0] for i in all_data], [i[1] for i in all_data])

    js_scatterplot_data_array = []
    if phenotype_min_max == None:
        for phenotype in js_scatterplot_data.keys():
            js_scatterplot_data_array.append(js_scatterplot_data[phenotype][1])
            js_scatterplot_data_array.append(js_scatterplot_data[phenotype][0])
    else:
        js_scatterplot_data_array.append(js_scatterplot_data[0])
        js_scatterplot_data_array.append(js_scatterplot_data[1])

    return json.dumps({
        "mutation_gene_expression": {
            "data": js_mutation_gene_expression_data,
            "stats": js_mutation_gene_expression_stats,
            "colors": js_mutation_gene_expression_colors,
            "name": js_mutation_gene_expression_name,
        },
        "bicluster_eigengene_expression": {
            "data": js_bicluster_eigengene_expression_data,
            "stats": js_bicluster_eigengene_expression_stats,
            "colors": js_bicluster_eigengene_expression_colors,
            "name": js_bicluster_eigengene_expression_name,
        },
        "residual": {
            "data": js_bicluster_residual_data,
            "stats": js_bicluster_residual_stats,
            "colors": js_bicluster_residual_colors,
            "name": js_bicluster_residual_name,
        },
        "scatter": {
            "data": js_scatterplot_data_array,
            "stats": js_scatterplot_stats,
            "regression": [
                [lowest_x, regression.coef_[0] * lowest_x + regression.intercept_],
                [highest_x, regression.coef_[0] * highest_x + regression.intercept_]
            ]
        },
    })

    db.close()

@app.route('/bicluster-expression-graph/<bicluster>/<phenotype_name>')
def bicluster_expression_graph(bicluster=None, phenotype_name="histology_WHO"):
    db = dbconn()
    c = db.cursor()
    c.execute("SELECT id FROM bicluster WHERE name=%s;", [bicluster])
    bc_pk = c.fetchone()[0]

    phenotype_min_max = get_phenotype_min_max(c, phenotype_name)

    # Prepare graph plotting data
    all_boxplot_data = cluster_data(c, bc_pk)

    patients = [item[0] for item in all_boxplot_data]

    c.execute("SELECT p.name, pd.phenotype_string, pd.phenotype_value FROM patient p"
        + " JOIN pheno_data pd ON pd.patient_id=p.id"
        + " JOIN phenotype pt ON pd.phenotype_id=pt.id"
        + " WHERE p.name IN (" + (','.join(map(lambda p: '\'%s\'' % p, patients))) + ") AND pt.name=%s;",
        [phenotype_name]
    )
    values = c.fetchall()
    string_ptmap = {patient: phenotype for patient, phenotype, _ in values}
    value_ptmap = {patient: value for patient, _, value in values}

    js_enrichment_quintiles = None
    js_enrichment_scatter = None
    if USE_PHENOTYPE_SCATTERPLOT == False or phenotype_min_max == None:
        enrichment_pvalues, min_enrichment_pvalue, max_enrichment_pvalue = subtype_enrichment(c, bc_pk, phenotype_name, phenotype_min_max)
        enrichment_data = []
        enrichment_colors = []

        for part in enrichment_pvalues:
            if phenotype_min_max == None:
                for phenotype in PHENOTYPES_DICT[phenotype_name]:
                    enrichment_data.append([phenotype, part[phenotype]])
                    enrichment_colors.append(GRAPH_COLOR_MAP[phenotype])
            else:
                index = 1
                for key in part.keys():
                    key_range = "[" + "{:.2f}".format(key[0]) + "," + "{:.2f}".format(key[1]) + ")"
                    if index == 5:
                        key_range = "[" + "{:.2f}".format(key[0]) + "," + "{:.2f}".format(key[1]) + "]"
                    
                    enrichment_data.append(["Phenotype Quintile " + str(index) + ": " + key_range, part[key]])
                    enrichment_colors.append(get_phenotype_color(phenotype_name, key[0], phenotype_min_max))
                    index = index + 1

        enrichment_upper = -math.log10(0.05/30.0)
        enrichment_lower = math.log10(0.05/30.0)
        enrich_perc20 = len(enrichment_data) / 5
        enrich_quintiles = [enrich_perc20 * i - 0.5 for i in range(1, 6)]

        js_enrichment_quintiles = {
            "quintiles": enrich_quintiles,
            "maxPValue": max_enrichment_pvalue,
            "minPValue": min_enrichment_pvalue,
            "upper": enrichment_upper,
            "lower": enrichment_lower,
            "data": enrichment_data,
            "colors": enrichment_colors,
        }
    else:
        js_scatterplot = {
            "name": "Values",
            "data": [],
            "color": GRAPH_COLOR_GRADIENTS[phenotype_name][1],
            "marker": {
                "symbol": "circle",
                "radius": 3,
            },
        }

        all_data = []
        highest_x = -10000
        lowest_x = 10000
        for i, item in enumerate(all_boxplot_data):
            x = item[3] # get the median
            y = value_ptmap[patients[i]]
            
            js_scatterplot["data"].append({
                "x": x,
                "y": y,
                "name": patients[i],
                "color": get_phenotype_color(phenotype_name, y, phenotype_min_max),
            })

            all_data.append([x, y])

            if x < lowest_x:
                lowest_x = x
            
            if x > highest_x:
                highest_x = x

        flattened_x = [i[0] for i in all_data]
        flattened_y = [i[1] for i in all_data]

        js_scatterplot_stats = stats.pearsonr(flattened_x, flattened_y)
        regression = LinearRegression().fit([(i[0],) for i in all_data], flattened_y)
            
        js_enrichment_scatter = {
            "biclusterName": "{} Median Expression".format(bicluster),
            "phenotypeName": "{} Value".format(PHENOTYPE_INDEX_TO_UINAME[phenotype_name]),
            "scatter": {
                "data": js_scatterplot,
                "stats": js_scatterplot_stats,
                "regression": [
                    [lowest_x, regression.coef_[0] * lowest_x + regression.intercept_],
                    [highest_x, regression.coef_[0] * highest_x + regression.intercept_]
                ]
            },
        }

    ratios_mean = 8.167528228975065 # TODO: have some way to precompute the mean for all the 4.5 million gene_expression cells

    js_boxplot_colors = []
    js_boxplot_legend = []
    if phenotype_min_max == None:
        phenotypes = [string_ptmap[patient] for patient in patients]
        for phenotype in PHENOTYPES_DICT[phenotype_name]:
            js_boxplot_legend.append({
                "name": phenotype,
                "color": GRAPH_COLOR_MAP[phenotype],
            })
    else:
        js_boxplot_legend.append({
            "name": "Lowest Value",
            "color": GRAPH_COLOR_GRADIENTS[phenotype_name][0],
        })

        js_boxplot_legend.append({
            "name": "Highest Value",
            "color": GRAPH_COLOR_GRADIENTS[phenotype_name][1],
        })
    
    js_boxplot_data = []
    for i, item in enumerate(all_boxplot_data):
        color = None
        if phenotype_min_max == None:
            color = GRAPH_COLOR_MAP[string_ptmap[patients[i]]]
            name = patients[i]
        else:
            color = get_phenotype_color(phenotype_name, value_ptmap[patients[i]], phenotype_min_max)
            name = patients[i] + " - " + PHENOTYPE_INDEX_TO_UINAME[phenotype_name] + ": " + "{:.2f}".format(value_ptmap[patients[i]])

        js_boxplot_colors.append(color)
        js_boxplot_data.append({
            "name": name,
            "low": item[1],
            "q1":  item[2],
            "median": item[3],
            "q3": item[4],
            "high": item[5],
            "fillColor": "%sB3" % color, # get the color and append a transparency value
        })
    
    # perc20 = len(in_data) / 5
    perc20 = len(all_boxplot_data) / 5
    quintiles = [perc20 * i for i in range(1, 6)]

    db.close()

    return json.dumps({
        "boxplotData": js_boxplot_data,
        "boxplotLegend": js_boxplot_legend,
        "boxplotColors": js_boxplot_colors,
        "boxplotRatiosMean": ratios_mean,
        "boxplotQuintiles": quintiles,
        "enrichmentQuintiles": js_enrichment_quintiles,
        "enrichmentScatter": js_enrichment_scatter,
    })

@app.route('/bicluster/<bicluster>/info/')
def bicluster_cytoscape(bicluster=None):
    db = dbconn()
    c = db.cursor()
    
    c.execute("SELECT id, name FROM bicluster WHERE name=%s;", [bicluster])
    bc_pk, bc_name = c.fetchone()

    r_cutoff = request.args.get('r_cutoff')
    if r_cutoff != None:
        r_cutoff = float(r_cutoff)
    else:
        r_cutoff = 0.0
    
    # Regulators
    elements = []
    elements.append({'data': { 'id': 'bc%d' % bc_pk, 'name': bc_name}, 'classes': 'bicluster' })

    regulators = []
    regulatorIds = []
    # c.execute("""SELECT g.id, g.symbol, tfr.action FROM tf_regulator tfr join gene g on tfr.gene_id=g.id WHERE tfr.bicluster_id=%s""", [bc_pk])
    c.execute("""SELECT tfr.id, g.symbol, tfr.r_value, tfr.p_value FROM tf_regulator tfr join gene g on tfr.tf_id=g.id WHERE tfr.bicluster_id=%s""", [bc_pk])
    tfs = list(c.fetchall())
    tfList = []
    for tf in tfs:
        if abs(tf[2]) < r_cutoff:
            continue
        
        known = 'No'

        tf_action = "Repressor"
        if tf[2] > 0:
            tf_action = "Activator"

        regulators.append(['TF', tf[0], tf[1], tf_action, known, "{:.5g}".format(tf[3]), "{:.5g}".format(tf[2])])
        regulatorIds.append(tf[0])
        tfList.append(tf[1])
        elements.append({'data': { 'id': 'reg%d' % tf[0], 'name': tf[1] }, 'classes': 'tf' })
        elements.append({'data': { 'id': 'tfbc%d' % tf[0], 'source': 'reg%d' % tf[0], 'target': 'bc%d' % bc_pk }, 'classes': "N/A Action" })

    c.execute("""SELECT mirna.id, mirna.name FROM mirna_regulator mr join mirna on mirna.id=mr.mirna_id WHERE mr.bicluster_id=%s""", [bc_pk])
    mirnas = list(c.fetchall())

    mirnaList = []
    for mirna in mirnas:
        if not mirna[0] in mirnaList:
            # our known is yes if we can find any mirna's given the mirna sql id
            known = 'No'
            c.execute("""SELECT COUNT(*) FROM mirna_prior WHERE mirna_id=%s""", [mirna[0]])
            mirna_count = c.fetchone()[0]
            if mirna_count != 0:
                known = 'Yes'
            
            regulators.append(['miRNA', mirna[0], mirna[1], 'Repressor', known])
            mirnaList.append(mirna[1])
            elements.append({'data': { 'id': 'reg%d' % mirna[0], 'name': mirna[1]}, 'classes': 'mirna' })
            elements.append({'data': { 'id': 'mirnabc%d' % mirna[0], 'source': 'reg%d' % mirna[0], 'target': 'bc%d' % bc_pk }, 'classes': 'repressor' })

    regulators = sorted(regulators, key=lambda name: name[1])

    # Get causal flows with bicluster
    c.execute("""SELECT id,somatic_mutation_id,regulator_id,regulator_type,bicluster_id,leo_nb_atob,mlogp_m_atob
FROM causal_flow WHERE bicluster_id=%s""", [bc_pk])
    tmp_cf = c.fetchall()

    for cf_pk, cf_som_mut_id, cf_reg_id, cf_reg_type, cf_bc_id, cf_leo, cf_mlogp in tmp_cf:
        if cf_reg_id not in regulatorIds:
            continue
        
        if cf_reg_type == 'tf':
            c.execute("""SELECT g.symbol FROM tf_regulator tr JOIN gene g ON tr.tf_id = g.id WHERE tr.id=%s""", [cf_reg_id])
            g1 = c.fetchone()[0]
        else:
            c.execute("""SELECT m.name FROM mirna_regulator mr JOIN mirna m ON mr.mirna_id = m.id WHERE mr.id=%s""", [cf_reg_id])
            g1 = c.fetchone()[0]

        if (cf_reg_type == 'tf' and g1 in tfList) or (cf_reg_type == 'mirna' and g1 in mirnaList):
            c.execute("""SELECT * FROM somatic_mutation WHERE id=%s""", [cf_som_mut_id])
            m1 = c.fetchone()
            if m1[3] == None:
                c.execute("""SELECT symbol FROM gene WHERE id=%s""", [m1[1]])
                mut = c.fetchone()[0]
            else:
                c.execute("""SELECT locus_name FROM locus WHERE id=%s """, [m1[3]]) # test: make sure this gets loci
                mut = c.fetchone()[0]

            elements.append({'data': { 'id': 'mut%d' % cf_som_mut_id, 'name': mut}, 'classes': 'genotype' }) # feature: differnet colors for genes vs loci
            elements.append({'data': { 'id': 'cf%d' % cf_pk, 'source': 'mut%d' % cf_som_mut_id, 'target': 'reg%d' % cf_reg_id } })

    # Hallmarks of Cancer
    c.execute("""SELECT hm.id,hm.name FROM hallmark hm join bic_hal bh on hm.id=bh.hallmark_id
WHERE bh.bicluster_id=%s""", [bc_pk])
    for hm_id, hm_name in c.fetchall():
        elements.append({'data': { 'id': 'hm%d' % hm_id, 'name': hm_name}, 'classes': 'hallmark' })
        elements.append({'data': { 'id': 'bchm%d' % hm_id, 'source': 'bc%d' % bc_pk, 'target': 'hm%d' % hm_id } })
    
    return json.dumps({
        "cytoscape": elements,
        "regulators": regulators,
    })

@app.route('/bicluster/<bicluster>')
def bicluster(bicluster=None):
    selectable_phenotypes = SELECTABLE_PHENOTYPES
    
    db = dbconn()
    c = db.cursor()
    c.execute("""SELECT id,name,var_exp_fpc,var_exp_fpc_p_value,survival,survival_p_value
FROM bicluster WHERE name=%s""", [bicluster])
    bc_pk, bc_name, bc_varexp_fpc, bc_varexp_fpc_pval, bc_survival, bc_survival_pval = c.fetchone()
    bic_info = {
        'pk': bc_pk,
        'name': bc_name,
        'varexp_fpc': bc_varexp_fpc,
        'varexp_fpc_pval': bc_varexp_fpc_pval,
        'survival': bc_survival,
        'survival_pval': bc_survival_pval,
        'varexp_flag': bc_varexp_fpc_pval <= 0.05,
        'survival_flag': bc_survival_pval <= 0.05
        }

    c.execute("""SELECT g.id, g.symbol, g.entrez FROM bic_gene bg join gene g on bg.gene_id=g.id where bg.bicluster_id=%s order by g.symbol""", [bc_pk])
    genes = list(c.fetchall())
    c.execute("""SELECT p.id, p.name FROM bic_pat bp join patient p on p.id=bp.patient_id where bp.bicluster_id=%s order by p.name""", [bc_pk])
    tumors = list(c.fetchall())
    # Replication
    c.execute("""SELECT * FROM replication WHERE bicluster_id=%s""", [bc_pk])
    tmp = list(c.fetchall())

    # replication datasets (add MESO TCGA in correctly)
    repConvert = {"MESO TCGA": "Hmeljak, et al. 2018"}
    repPubmed = {"MESO TCGA": "30322867"}
    replication = []

    replicated = [0, 0]
    bic_info['repl_coexp'] = False
    bic_info['repl_survival'] = False

    for i in tmp:
        tmp1 = [0,0]
        if bic_info['varexp_fpc_pval'] <= 0.05 and float(i[4])<=0.05:
            tmp1[0] = 1
            replicated[0] = 1
            bic_info['repl_coexp'] = True

        if (( bic_info['survival'] > 0 and (float(i[5]) > 0 and float(i[6])<=0.05)) or
            ( bic_info['survival'] < 0 and (float(i[5]) < 0 and float(i[6])<=0.05))):
            tmp1[1] = 1
            replicated[1] = 1
            bic_info['repl_survival'] = True

        replication.append(list(i)+[repConvert[i[2]], repPubmed[i[2]]]+tmp1)

    regulators = []
    # c.execute("""SELECT g.id, g.symbol, tfr.action FROM tf_regulator tfr join gene g on tfr.gene_id=g.id WHERE tfr.bicluster_id=%s""", [bc_pk])
    c.execute("""SELECT tfr.id, g.symbol, tfr.r_value, tfr.p_value FROM tf_regulator tfr join gene g on tfr.tf_id=g.id WHERE tfr.bicluster_id=%s""", [bc_pk])
    tfs = list(c.fetchall())
    tfList = []
    for tf in tfs:
        tfList.append(tf[1])

    c.execute("""SELECT mirna.id, mirna.name FROM mirna_regulator mr join mirna on mirna.id=mr.mirna_id WHERE mr.bicluster_id=%s""", [bc_pk])
    mirnas = list(c.fetchall())

    mirnaList = []
    for mirna in mirnas:
        if not mirna[0] in mirnaList:
            mirnaList.append(mirna[1])

    regulators = sorted(regulators, key=lambda name: name[1])

    # Get causal flows with bicluster
    c.execute("""SELECT id,somatic_mutation_id,regulator_id,regulator_type,bicluster_id,leo_nb_atob,mlogp_m_atob
FROM causal_flow WHERE bicluster_id=%s""", [bc_pk])
    tmp_cf = c.fetchall()
    causalFlows = []

    for cf_pk, cf_som_mut_id, cf_reg_id, cf_reg_type, cf_bc_id, cf_leo, cf_mlogp in tmp_cf:
        if cf_reg_type == 'tf':
            c.execute("""SELECT g.symbol FROM tf_regulator tr JOIN gene g ON tr.tf_id = g.id WHERE tr.id=%s""", [cf_reg_id])
            g1 = c.fetchone()[0]
        else:
            c.execute("""SELECT m.name FROM mirna_regulator mr JOIN mirna m ON mr.mirna_id = m.id WHERE mr.id=%s""", [cf_reg_id])
            g1 = c.fetchone()[0]

        if (cf_reg_type == 'tf' and g1 in tfList) or (cf_reg_type == 'mirna' and g1 in mirnaList):
            c.execute("""SELECT * FROM somatic_mutation WHERE id=%s""", [cf_som_mut_id])
            m1 = c.fetchone()
            if m1[3] == None:
                c.execute("""SELECT symbol FROM gene WHERE id=%s""", [m1[1]])
                mut = c.fetchone()[0]
            else:
                c.execute("""SELECT locus_name FROM locus WHERE id=%s """, [m1[3]]) # test: make sure this gets loci
                mut = c.fetchone()[0]
            
            c.execute("""
                SELECT COUNT(*) > 0 FROM bueno_deep_filter bd
                    JOIN somatic_mutation sm ON bd.somatic_mutation_id = sm.id
                    LEFT JOIN gene g2 ON sm.ext_id = g2.id
                    LEFT JOIN locus l on sm.locus_id = l.id
                    WHERE g2.symbol=%s OR l.locus_name=%s;
            """, [mut, mut])
            has_graphs = c.fetchone()[0]
            if has_graphs > 0:
                graph_button_style = "inline"
            else:
                graph_button_style = "none"

            causalFlows.append([mut, g1, graph_button_style])

    causalFlows = sorted(causalFlows, key=lambda mutation: mutation[0])

    # Hallmarks of Cancer
    c.execute("""SELECT hm.id,hm.name FROM hallmark hm join bic_hal bh on hm.id=bh.hallmark_id
WHERE bh.bicluster_id=%s""", [bc_pk])
    hallmarks = []
    hallmark_to_image = {
        "Sustained angiogenesis": 5,
        "Insensitivity to antigrowth signals": 10,
        "Evading apoptosis": 3,
        "Evading immune detection": 9,
        "Tissue invasion and metastasis": 6,
        "Self sufficiency in growth signals": 1,
        "Tumor promoting inflammation": 7,
        "Reprogramming energy metabolism": 2,
        "Genome instability and mutation": 4,
        "Limitless replicative potential": 8,
    }
    hallmark_image_class = {
        1: "disabled",
        2: "disabled",
        3: "disabled",
        4: "disabled",
        5: "disabled",
        6: "disabled",
        7: "disabled",
        8: "disabled",
        9: "disabled",
        10: "disabled",
    }
    for hm_id, hm_name in c.fetchall():
        hallmarks.append([hm_name, convert[hm_name] ])
        hallmark_image_class[hallmark_to_image[hm_name]] = ""

    # GO
    c.execute(
        """SELECT go_bp.id, go_bp.go_id, go_bp.name FROM bic_go
        JOIN go_bp ON go_bp.id = bic_go.go_bp_id
        WHERE bic_go.bicluster_id = %s""",
        [bc_pk]
    )
    tmps = list(c.fetchall())
    gobps = []
    for gobp in tmps:
        c.execute(
            """SELECT DISTINCT gene.symbol FROM go_gene, gene, bic_gene
            WHERE go_gene.go_bp_id = %s AND bic_gene.bicluster_id = %s
            AND go_gene.gene_id = gene.id
            AND gene.id = bic_gene.gene_id
            ORDER BY gene.symbol""",
            [gobp[0], bc_pk]
        )
        gobps.append(list(gobp) + [[row[0] for row in c.fetchall()]])

    db.close()

    return render_template('bicluster.html', **locals())


@app.route('/search')
def search():
    gene = request.args.get('gene')
    db = dbconn()
    c = db.cursor()
    bme_type ='gene'
    if not gene:
        return render_template('index.html')
    if gene.find('hsa-')==-1:
        c.execute("""SELECT symbol FROM gene WHERE symbol=%s""", [gene])
        geneData = c.fetchall()
    else:
        c.execute("""SELECT name FROM mirna WHERE name=%s""", [gene])
        geneData = c.fetchall()
        bme_type = 'mirna'
    db.close()

    if len(geneData)==0:
        return render_template('index.html')

    symbol = geneData[0][0]
    if bme_type == 'gene':
        return redirect(url_for('gene', symbol=symbol))
    else:
        return redirect(url_for('mirna', symbol=symbol))


def __get_muts(c, gene_pk, symbol):
    # Get causal flows downstream of mutation in gene
    muts = {}
    c.execute("""SELECT * FROM somatic_mutation WHERE locus_id IS NULL AND ext_id=%s""", [gene_pk]) # discrep. (we never have a mutation type of gene). test: make sure that we only get genes (not loci) from this query
    tmp_muts = c.fetchall()
    if len(tmp_muts)==1:
        muts['name'] = symbol
        c.execute("""SELECT * FROM causal_flow WHERE somatic_mutation_id=%s""", [tmp_muts[0][0]])
        tmp_cf = c.fetchall()
        muts['flows'] = 0
        muts['regs'] = []
        muts['tfs'] = []
        muts['miRNAs'] = []
        muts['biclusters'] = []
        muts['data'] = []
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
            if not g1=='':
                c.execute("""SELECT name, survival, survival_p_value FROM bicluster WHERE id=%s""", [cf1[4]])
                b1 = c.fetchall()[0]
                if not b1 in muts['biclusters']:
                    muts['biclusters'].append(b1[0])
                c.execute("SELECT hm.name FROM hallmark hm join bic_hal bh on hm.id=bh.hallmark_id  WHERE bh.bicluster_id=%s",
                          [cf1[4]])

                tmp1 = c.fetchall()
                h1 = list(set([convert[i[0]] for i in tmp1]))
                h2 = [[i[0],convert[i[0]]] for i in tmp1]
                muts['data'].append([symbol, g1, b1[0], b1[1], b1[2], h2])
    return muts


def __get_regulators(c, symbol, bme_type):
    regs = {}
    tmp_regs = c.fetchall()
    if len(tmp_regs)>0:
        regs['name'] = symbol
        regs['biclusters'] = len(set([i[1] for i in tmp_regs]))
        regs['data'] = []
        # Collect all biclusters downstream regulated by TF or miRNA
        for reg in tmp_regs:
            action = 'Rep.'
            if bme_type =='gene' and reg[3] == 'activator':
                action = 'Act.'
            c.execute("""SELECT name, survival, survival_p_value FROM bicluster WHERE id=%s""", [reg[1]])
            b1 = c.fetchall()[0]
            c.execute("SELECT hm.name FROM hallmark hm join bic_hal bh on hm.id=bh.hallmark_id WHERE bh.bicluster_id=%s",
                      [reg[1]])
            tmp1 = c.fetchall()
            h1 = list(set([convert[i[0]] for i in tmp1]))
            h2 = [[i[0],convert[i[0]]] for i in tmp1]
            regs['data'].append([symbol, action, b1[0], b1[1], b1[2], h2])
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
    db = dbconn()
    c = db.cursor()
    c.execute("""SELECT id FROM gene WHERE symbol=%s""", [symbol])
    gene_pk = c.fetchone()[0]
    muts = __get_muts(c, gene_pk, symbol)
    c.execute("""SELECT * FROM tf_regulator WHERE tf_id=%s""", [gene_pk])
    regs = __get_regulators(c, symbol, "gene")

    # Get biclusters that gene resides
    bics = {}
    c.execute("SELECT * FROM bic_gene bg join bicluster b on bg.bicluster_id=b.id where bg.gene_id=%s", [gene_pk])
    tmp_bics = c.fetchall()
    if len(tmp_bics) > 0:
        bics['name'] = gene
        bics['biclusters'] = len(tmp_bics)
        bics['data'] = []
        for bic1 in tmp_bics:
            c.execute("SELECT hm.name FROM bic_hal bh JOIN hallmark hm on bh.hallmark_id=hm.id WHERE bh.bicluster_id=%s",
                      [bic1[3]])
            tmp1 = c.fetchall()
            h1 = list(set([convert[i[0]] for i in tmp1]))
            h2 = [[i[0],convert[i[0]]] for i in tmp1]
            bics['data'].append([bic1[4], bic1[5], bic1[6], bic1[7], bic1[8], h2])
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
