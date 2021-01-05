"""
handles the endpoints:

/bicluster-causal-analysis/<mutation>/<regulator>/<bicluster>/<phenotype_name>
"""

import json
import numpy as np
import math
from scipy import stats
from sklearn.linear_model import LinearRegression
from phenotype import get_phenotype_min_max, get_phenotype_color

from database import dbconn
from flask import Blueprint, request

from constants import GRAPH_COLOR_MAP, GRAPH_COLOR_GRADIENTS, PHENOTYPE_INDEX_TO_UINAME, PHENOTYPES_DICT

causal_analysis_page = Blueprint("causal_analysis_page", __name__, template_folder="templates")

@causal_analysis_page.route('/bicluster-causal-analysis/<mutation>/<regulator>/<bicluster>/<phenotype_name>')
# mutation and regulator are gene symbol names, bicluster is the bicluster name
def bicluster_causal_analysis_page(mutation=None, regulator=None, bicluster=None, phenotype_name="histology_WHO"):
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
	js_scatterplot_data[0] = {
		"name": "Mutated",
		"data": [],
	}

	js_scatterplot_data[1] = {
		"name": "WT",
		"data": [],
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
			js_scatterplot_data[boolean]["data"].append({
				"x": x,
				"y": y,
				"phenotype": string_by_patient[patient],
				"name": patient,
				"color": GRAPH_COLOR_MAP[string_by_patient[patient]],
			})
		else:
			js_scatterplot_data[boolean]["data"].append({
				"x": x,
				"y": y,
				"z": value_by_patient[patient],
				"phenotype": PHENOTYPE_INDEX_TO_UINAME[phenotype_name],
				"name": patient,
				"color": get_phenotype_color(phenotype_name, value_by_patient[patient], phenotype_min_max),
			})

		all_data.append([x, y])

		if x < lowest_x:
			lowest_x = x
		
		if x > highest_x:
			highest_x = x
	
	js_scatterplot_stats = stats.pearsonr([i[0] for i in all_data], [i[1] for i in all_data])

	js_scatterplot_data_array = []
	js_scatterplot_data_array.append(js_scatterplot_data[0])
	js_scatterplot_data_array.append(js_scatterplot_data[1])

	js_scatterplot_legend = []
	if phenotype_min_max == None:
		for phenotype in PHENOTYPES_DICT[phenotype_name]:
			js_scatterplot_legend.append({
				"name": phenotype,
				"color": GRAPH_COLOR_MAP[phenotype],
			})
	else:
		js_scatterplot_legend.append({
			"name": "Lowest Phenotype",
			"color": GRAPH_COLOR_GRADIENTS[phenotype_name][0],
		})

		js_scatterplot_legend.append({
			"name": "Highest Phenotype",
			"color": GRAPH_COLOR_GRADIENTS[phenotype_name][1],
		})

	db.close()

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
			"legend": js_scatterplot_legend,
			"stats": js_scatterplot_stats,
			"regression": [
				[lowest_x, regression.coef_[0] * lowest_x + regression.intercept_],
				[highest_x, regression.coef_[0] * highest_x + regression.intercept_]
			]
		},
	})