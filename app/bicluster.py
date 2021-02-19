"""
handles the endpoints:

/bicluster/<bicluster>
/bicluster/<bicluster>/info
/bicluster-expression-graph/<bicluster>/<phenotype_name>
"""

import json
import math

from database import dbconn
from scipy import stats
from sklearn.linear_model import LinearRegression
from jacks import get_jacks_data, get_meso1_data

from graph_common import cluster_data, subtype_enrichment
from phenotype import get_phenotype_min_max, get_phenotype_color

from constants import HALLMARK_TO_ICON, SELECTABLE_PHENOTYPES, SELECTABLE_PHENOTYPES_BLACKLIST, USE_PHENOTYPE_SCATTERPLOT, PHENOTYPE_INDEX_TO_UINAME, GRAPH_COLOR_MAP, GRAPH_COLOR_GRADIENTS, PHENOTYPES_DICT
from database import dbconn
from flask import Blueprint, render_template, request

bicluster_page = Blueprint("bicluster_page", __name__, template_folder="templates")

@bicluster_page.route('/bicluster/<bicluster>')
def bicluster(bicluster=None):
	selectable_phenotypes = SELECTABLE_PHENOTYPES
	
	db = dbconn()
	c = db.cursor()
	c.execute(
		"""SELECT id, name, var_exp_fpc, var_exp_fpc_p_value, survival, survival_p_value
		FROM bicluster
		WHERE name=%s""",
		[bicluster]
	)
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
	
	# select phenotype/bicluster relationships where the p-value is at 95% significance
	c.execute(
		"""SELECT bps.p_value, p.name FROM bicluster_phenotype_significance bps
		JOIN bicluster b ON b.id = bps.bicluster_id
		JOIN phenotype p ON p.id = bps.phenotype_id
		WHERE b.id = %s AND p_value < 0.05;""",
		[bc_pk]
	)
	data = c.fetchall()
	phenotype_to_p_value = {datum[1]: datum[0] for datum in data}

	expression_selectable_phenotypes = []
	for uiName, dbName in SELECTABLE_PHENOTYPES:
		if (
			dbName in SELECTABLE_PHENOTYPES_BLACKLIST
			or dbName in phenotype_to_p_value
		):
			expression_selectable_phenotypes.append((uiName, dbName))

	c.execute(
		"""SELECT g.id, g.symbol, g.entrez
		FROM bic_gene bg
		JOIN gene g ON bg.gene_id=g.id
		WHERE bg.bicluster_id=%s
		ORDER BY g.symbol""",
		[bc_pk]
	)
	genes = list(c.fetchall())
	c.execute(
		"""SELECT p.id, p.name
		FROM bic_pat bp
		JOIN patient p ON p.id=bp.patient_id
		WHERE bp.bicluster_id=%s
		ORDER BY p.name""",
		[bc_pk]
	)
	tumors = list(c.fetchall())
	# Replication
	c.execute(
		"""SELECT *
		FROM replication
		WHERE bicluster_id=%s""",
		[bc_pk]
	)
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

	info = get_bicluster_info(bicluster, 0)
	# set variables for template to use
	causal_flows = info["causal_flows"]
	regulators = info["regulators"]
	hallmark_image_class = info["hallmark_image_class"]
	hallmarks = info["hallmarks"]

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
	
	# fetch jacks data
	jacks_data = get_jacks_data(bc_pk)

	# fetch MESO1 data
	meso1_data = get_meso1_data(bc_pk)

	db.close()

	return render_template('bicluster.html', **locals())

@bicluster_page.route('/bicluster/<bicluster>/info/')
def bicluster_cytoscape(bicluster=None):
	db = dbconn()
	c = db.cursor()

	r_cutoff = request.args.get('r_cutoff')
	if r_cutoff != None:
		r_cutoff = float(r_cutoff)
	else:
		r_cutoff = 0.0
	
	info = get_bicluster_info(bicluster, r_cutoff)
	
	return json.dumps({
		"cytoscape": info["cytoscape"],
		"regulators": info["regulators"],
	})

def get_bicluster_info(bicluster, r_cutoff):
	db = dbconn()
	c = db.cursor()
	
	c.execute("SELECT id, name FROM bicluster WHERE name=%s;", [bicluster])
	bc_pk, bc_name = c.fetchone()
	
	# Regulators
	elements = []
	elements.append({'data': { 'id': 'bc%d' % bc_pk, 'name': bc_name}, 'classes': 'bicluster' })

	regulators = []
	c.execute(
		"""SELECT tfr.id, g.symbol, tfr.r_value, tfr.p_value
		FROM tf_regulator tfr
		JOIN gene g ON tfr.tf_id = g.id
		WHERE tfr.bicluster_id=%s""",
		[bc_pk]
	)
	tfs = list(c.fetchall())
	tfList = []
	for tf in tfs:
		if abs(tf[2]) < r_cutoff:
			continue
		
		known = 'No'
		'''
		c.execute("""SELECT * FROM tf_crispr WHERE gene_id=%s""", [tf[0]]) # we will have this table soon
		for crispr in c.fetchall():
			if float(crispr[4])<=0.05:
				known = 'Yes'
		'''

		tf_action = "Repressor"
		if tf[2] > 0:
			tf_action = "Activator"

		regulators.append(['TF', tf[0], tf[1], tf_action, known, "{:.5g}".format(tf[3]), "{:.5g}".format(tf[2])])
		tfList.append(tf[1])
		elements.append({'data': { 'id': 'reg%d' % tf[0], 'name': tf[1] }, 'classes': 'tf' })
		elements.append({'data': { 'id': 'tfbc%d' % tf[0], 'source': 'reg%d' % tf[0], 'target': 'bc%d' % bc_pk }, 'classes': "N/A Action" })

	c.execute(
		"""SELECT mirna.id, mirna.name, mr.id, mr.p_value, mr.r_value
    FROM mirna_regulator mr
    JOIN mirna ON mirna.id = mr.mirna_id
    WHERE mr.bicluster_id = %s;""",
		[bc_pk]
	)
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
			
			regulators.append(['miRNA', mirna[0], mirna[1], 'Repressor', known, mirna[3], mirna[4]])
			mirnaList.append(mirna[1])
			elements.append({'data': { 'id': 'reg%d' % mirna[2], 'name': mirna[1]}, 'classes': 'mirna' })
			elements.append({'data': { 'id': 'mirnabc%d' % mirna[2], 'source': 'reg%d' % mirna[2], 'target': 'bc%d' % bc_pk }, 'classes': 'repressor' })

	regulators = sorted(regulators, key=lambda name: name[1])

	# Get causal flows with bicluster
	c.execute(
		"""SELECT id, somatic_mutation_id, regulator_id, regulator_type, bicluster_id, leo_nb_atob, mlogp_m_atob
    FROM causal_flow WHERE bicluster_id = %s;""",
		[bc_pk]
	)
	tmp_cf = c.fetchall()
	causal_flows = []
	for cf_pk, cf_som_mut_id, cf_reg_id, cf_reg_type, cf_bc_id, cf_leo, cf_mlogp in tmp_cf:
		if cf_reg_type == 'tf':
			c.execute(
				"""SELECT g.symbol
				FROM tf_regulator tr
				JOIN gene g ON tr.tf_id = g.id
				WHERE tr.id=%s""",
				[cf_reg_id]
			)
			g1 = c.fetchone()[0]
		else:
			c.execute(
				"""SELECT m.name
				FROM mirna_regulator mr
				JOIN mirna m ON mr.mirna_id = m.id
				WHERE mr.id = %s""",
				[cf_reg_id]
			)
			g1 = c.fetchone()[0]

		if (cf_reg_type == 'tf' and g1 in tfList) or (cf_reg_type == 'mirna' and g1 in mirnaList):
			c.execute(
				"""SELECT *
				FROM somatic_mutation
				WHERE id = %s""",
				[cf_som_mut_id]
			)
			m1 = c.fetchone()
			if m1[3] == None:
				c.execute(
					"""SELECT symbol
					FROM gene
					WHERE id = %s""",
					[m1[1]]
				)
				mut = c.fetchone()[0]
			else:
				c.execute(
					"""SELECT locus_name
					FROM locus WHERE id = %s """,
					[m1[3]]
				) # test: make sure this gets loci
				mut = c.fetchone()[0]
			
			c.execute(
				"""SELECT COUNT(*) > 0
				FROM bueno_deep_filter bd
				JOIN somatic_mutation sm ON bd.somatic_mutation_id = sm.id
				LEFT JOIN gene g2 ON sm.ext_id = g2.id
				LEFT JOIN locus l on sm.locus_id = l.id
				WHERE g2.symbol = %s OR l.locus_name = %s;""",
				[mut, mut]
			)
			has_graphs = c.fetchone()[0]
			if has_graphs > 0 and cf_reg_type == "tf":
				graph_button_style = "inline"
			else:
				graph_button_style = "none"

			causal_flows.append([mut, g1, graph_button_style])

			elements.append({'data': { 'id': 'mut%d' % cf_som_mut_id, 'name': mut}, 'classes': 'genotype' }) # feature: differnet colors for genes vs loci
			elements.append({'data': { 'id': 'cf%d' % cf_pk, 'source': 'mut%d' % cf_som_mut_id, 'target': 'reg%d' % cf_reg_id } })

	causal_flows = sorted(causal_flows, key=lambda mutation: mutation[0])

	# Hallmarks of Cancer
	c.execute(
		"""SELECT hm.id, hm.name
		FROM hallmark hm
		JOIN bic_hal bh ON hm.id = bh.hallmark_id
		WHERE bh.bicluster_id = %s""",
		[bc_pk]
	)
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
		hallmarks.append([hm_name, HALLMARK_TO_ICON[hm_name] ])
		elements.append({'data': { 'id': 'hm%d' % hm_id, 'name': hm_name}, 'classes': 'hallmark' })
		elements.append({'data': { 'id': 'bchm%d' % hm_id, 'source': 'bc%d' % bc_pk, 'target': 'hm%d' % hm_id } })
		hallmark_image_class[hallmark_to_image[hm_name]] = ""
	
	return {
		"cytoscape": elements,
		"regulators": regulators,
		"hallmarks": hallmarks,
		"hallmark_image_class": hallmark_image_class,
		"causal_flows": causal_flows,
	}

@bicluster_page.route('/bicluster-expression-graph/<bicluster>/<phenotype_name>')
def bicluster_expression_graph(bicluster=None, phenotype_name="histology_WHO"):
	db = dbconn()
	c = db.cursor()
	c.execute("SELECT id FROM bicluster WHERE name=%s;", [bicluster])
	bc_pk = c.fetchone()[0]

	phenotype_min_max = get_phenotype_min_max(c, phenotype_name)

	# Prepare graph plotting data
	all_boxplot_data = cluster_data(c, bc_pk)

	patients = [item[0] for item in all_boxplot_data]

	c.execute("""
		SELECT p.name, pd.phenotype_string, pd.phenotype_value FROM patient p
		JOIN pheno_data pd ON pd.patient_id=p.id
		JOIN phenotype pt ON pd.phenotype_id=pt.id
		WHERE pt.name=%s;""",
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
		c.execute("""
			SELECT e.value, pt.name FROM eigengene e
			JOIN patient pt on e.patient_id = pt.id
			JOIN bicluster b ON e.bicluster_id = b.id
			WHERE b.name = %s;""",
			[bicluster]
		)
		eigengene_values = c.fetchall()
		
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
		for eigengene_value, patient in eigengene_values:
			if patient not in value_ptmap:
				continue

			x = eigengene_value
			y = value_ptmap[patient]
			
			js_scatterplot["data"].append({
				"x": x,
				"y": y,
				"name": patient,
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
			"biclusterName": "{} Eigengene Expression".format(bicluster),
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