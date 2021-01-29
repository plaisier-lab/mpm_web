import base64
import json

from constants import HALLMARKS, HALLMARK_TO_ICON
from database import dbconn
from flask import Blueprint, render_template, request, redirect, url_for

search_page = Blueprint("search_page", __name__, template_folder="templates")

@search_page.route('/search')
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

	if gene.find('hsa-') == -1:
		c.execute("""SELECT symbol FROM gene WHERE symbol=%s""", [gene])
		geneData = c.fetchall()
	else:
		c.execute("""SELECT name FROM mirna WHERE name=%s""", [gene])
		geneData = c.fetchall()
		bme_type = 'mirna'
	db.close()

	if len(geneData) == 0 and phenotype == "None" and hallmark == "None":
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

	if len(geneData) != 0:
		symbol = geneData[0][0]
		if bme_type == 'gene':
			return redirect(url_for('.gene', symbol=symbol) + query_statement)
		else:
			return redirect(url_for('.mirna', symbol=symbol) + query_statement)
	else:
		return redirect(url_for('.advanced_search') + query_statement)


def __get_muts(c, gene_pk, gene, hallmark_search=None, phenotype_search=None):
	# Get causal flows downstream of mutation in gene
	mutations = {
		"name": gene,
		"flows": 0,
		"regulators": [],
		"tfs": [],
		"miRNAs": [],
		"biclusters": [],
		"causal_flows": [],
	}

	c.execute(
		"""SELECT *
		FROM somatic_mutation
		WHERE locus_id IS NULL AND ext_id=%s""",
		[gene_pk]
	) # discrep. (we never have a mutation type of gene). test: make sure that we only get genes (not loci) from this query
	tmp_muts = c.fetchall()
	if len(tmp_muts)==1:
		c.execute(
			"""SELECT *
			FROM causal_flow
			WHERE somatic_mutation_id=%s""",
			[tmp_muts[0][0]]
		)
		causal_flow_results = c.fetchall()
		for causal_flow in causal_flow_results:
			regulator = ''
			if causal_flow[3] == "tf": # try finding tf regulator
				c.execute(
					"""SELECT *
					FROM tf_regulator
					WHERE id=%s""",
					[causal_flow[2]]
				)
				tf = c.fetchall()
				if len(tf) > 0:
					c.execute(
						"""SELECT symbol
						FROM gene
						WHERE id=%s""",
						[tf[0][2]]
					)
					regulator = c.fetchall()[0][0]
			else: # try finding mirna regulator
				c.execute(
					"""SELECT *
					FROM mirna_regulator
					WHERE id=%s""",
					[causal_flow[2]]
				)
				mirna = c.fetchall()
				if len(mirna) > 0:
					c.execute(
						"""SELECT name
						FROM mirna
						WHERE id=%s""",
						[mirna[0][2]]
					)
					regulator = c.fetchall()[0][0]
			
			# print the bicluster in the table
			if not regulator == '':
				c.execute(
					"""SELECT id, name, survival, survival_p_value
					FROM bicluster
					WHERE id=%s""",
					[causal_flow[4]]
				)
				bicluster = c.fetchall()[0]

				# hallmark handling
				c.execute(
					"""SELECT hm.name
					FROM hallmark hm
					JOIN bic_hal bh ON hm.id=bh.hallmark_id
					WHERE bh.bicluster_id=%s""",
					[causal_flow[4]]
				)
				ordered_hallmarks = [i[0] for i in c.fetchall()]
				
				# check if the searched phenotype is significant
				phenotype_result = []
				if phenotype_search != None:
					c.execute(
						"""SELECT bps.r_value, bps.p_value
						FROM bicluster b
						JOIN bicluster_phenotype_significance bps ON bps.bicluster_id = b.id
						JOIN phenotype p ON bps.phenotype_id = p.id
						WHERE b.id = %s AND p.long_name = %s AND bps.p_value < 0.05;""",
						[bicluster[0], phenotype_search]
					)
					phenotype_result = c.fetchone()
					phenotype_result = [] if phenotype_result == None else phenotype_result

				# filter hallmarks using advanced search parameters
				if (
					has_hallmark_search(ordered_hallmarks, hallmark_search)
					and (
						(
							phenotype_search != None
							and len(phenotype_result) > 0
						)
						or phenotype_search == None
					)
				):
					if causal_flow[3] == "tf":
						if not regulator in mutations['regulators']:
							mutations['regulators'].append(regulator)
							mutations['tfs'].append(regulator)
					else:	
						if not regulator in mutations['regulators']:
							mutations['regulators'].append(regulator)
							mutations['miRNAs'].append(regulator)
					
					mutations['causal_flows'].append({
						"hallmarks": [(i, HALLMARK_TO_ICON[i]) for i in ordered_hallmarks],
						"data": [gene, regulator, bicluster[1], bicluster[2], bicluster[3]],
						"phenotype": phenotype_result,
					})

					if bicluster[1] not in mutations['biclusters']:
						mutations['biclusters'].append(bicluster[1])
	return mutations

def __get_regulators(c, symbol, bme_type, hallmark_search=None, phenotype_search=None):
	regulators = []
	result = c.fetchall()
	if len(result) > 0:
		# Collect all biclusters downstream regulated by TF or miRNA
		for regulator in result:
			action = 'Rep.'
			if bme_type =='gene' and regulator[3] == 'activator':
				action = 'Act.'
			c.execute(
				"""SELECT id, name, survival, survival_p_value
				FROM bicluster
				WHERE id=%s""",
				[regulator[1]]
			)
			bicluster = c.fetchall()[0]

			c.execute(
				"""SELECT hm.name
				FROM hallmark hm
				JOIN bic_hal bh ON hm.id=bh.hallmark_id
				WHERE bh.bicluster_id=%s""",
				[regulator[1]]
			)
			ordered_hallmarks = [i[0] for i in c.fetchall()]

			# check if the searched phenotype is significant
			phenotype_result = []
			if phenotype_search != None:
				c.execute(
					"""SELECT bps.r_value, bps.p_value
					FROM bicluster b
					JOIN bicluster_phenotype_significance bps ON bps.bicluster_id = b.id
					JOIN phenotype p ON bps.phenotype_id = p.id
					WHERE b.id = %s AND p.long_name = %s AND bps.p_value < 0.05;""",
					[bicluster[0], phenotype_search]
				)
				phenotype_result = c.fetchone()
				phenotype_result = [] if phenotype_result == None else phenotype_result

			# filter hallmarks using advanced search parameters
			if (
				has_hallmark_search(ordered_hallmarks, hallmark_search)
				and (
					(
						phenotype_search != None
						and len(phenotype_result) > 0
					)
					or phenotype_search == None
				)
			):
				regulators.append({
					"hallmarks": [(i, HALLMARK_TO_ICON[i]) for i in ordered_hallmarks],
					"bicluster": [symbol, action, bicluster[1], bicluster[2], bicluster[3]],
					"phenotype": phenotype_result,
				})
	return regulators

def has_hallmark_search(ordered_hallmarks, hallmark_search):
	has_hallmark = True # default to true if we have no hallmark search
	if hallmark_search != None and len(hallmark_search) > 0:
		has_hallmark = False

		# make sure ordered_hallmarks has at least one element in the hallmark_search
		for name in ordered_hallmarks:
			index = HALLMARKS.index(name) + 1
			if index in hallmark_search:
				has_hallmark = True
				break
	
	return has_hallmark

def decode_hallmark_query(query):
	# decode the hallmark
	hallmarks = []
	if query:
		hallmarks = json.loads(base64.b64decode(query))

		# verify hallmarks
		if type(hallmarks) != list:
			hallmarks = []
		
		verified = True
		for number in hallmarks:
			if type(number) != int or (
				number <= 0
				or number >= 11
			):
				verified = False

		if not verified:
			hallmarks = []
	
	return hallmarks

@search_page.route('/mirna')
@search_page.route('/mirna/<symbol>')
def mirna(symbol=None, defaults={'symbol': None}):
	phenotype_query = request.args.get('phenotype')
	phenotype_query = phenotype_query if phenotype_query else None
	hallmarks_query = request.args.get('hallmark')
	hallmarks_query = hallmarks_query if hallmarks_query else None

	hallmarks = decode_hallmark_query(hallmarks_query)
	
	# Get biclusters regulated by mirna
	db = dbconn()
	c = db.cursor()
	c.execute("""SELECT id FROM mirna WHERE name=%s""", [symbol])
	mirna_pk = c.fetchone()[0]
	mutations = __get_muts(c, mirna_pk, symbol)
	c.execute("""SELECT * FROM mirna_regulator WHERE mirna_id=%s""", [mirna_pk])
	regulators = __get_regulators(c, symbol, "mirna", hallmark_search=hallmarks, phenotype_search=phenotype_query)
	db.close()
	return render_template('search.html', gene=symbol, mutations={}, regulators=regulators, biclusters={})


@search_page.route('/gene')
@search_page.route('/gene/<symbol>')
def gene(symbol=None, defaults={'symbol': None}):
	phenotype_query = request.args.get('phenotype')
	phenotype_query = phenotype_query if phenotype_query else None
	hallmarks_query = request.args.get('hallmark')
	hallmarks_query = hallmarks_query if hallmarks_query else None

	hallmarks = decode_hallmark_query(hallmarks_query)
	
	db = dbconn()
	c = db.cursor()
	c.execute("""SELECT id FROM gene WHERE symbol=%s""", [symbol])
	gene_pk = c.fetchone()[0]
	mutations = __get_muts(c, gene_pk, symbol, hallmark_search=hallmarks, phenotype_search=phenotype_query)

	c.execute("""SELECT * FROM tf_regulator WHERE tf_id=%s""", [gene_pk])
	regulators = __get_regulators(c, symbol, "gene", hallmark_search=hallmarks, phenotype_search=phenotype_query)

	# Get biclusters that gene resides
	biclusters = []
	c.execute(
		"""SELECT * FROM bic_gene bg
		JOIN bicluster b ON bg.bicluster_id=b.id
		WHERE bg.gene_id=%s""",
		[gene_pk]
	)
	tmp_bics = c.fetchall()
	if len(tmp_bics) > 0:
		for bicluster in tmp_bics:
			c.execute(
				"""SELECT hm.name
				FROM bic_hal bh
				JOIN hallmark hm ON bh.hallmark_id=hm.id
				WHERE bh.bicluster_id=%s""",
				[bicluster[3]]
			)
			ordered_hallmarks = [i[0] for i in c.fetchall()]

			# check if the searched phenotype is significant
			phenotype_result = []
			if phenotype_query != None:
				c.execute(
					"""SELECT bps.r_value, bps.p_value
					FROM bicluster b
					JOIN bicluster_phenotype_significance bps ON bps.bicluster_id = b.id
					JOIN phenotype p ON bps.phenotype_id = p.id
					WHERE b.id = %s AND p.long_name = %s AND bps.p_value < 0.05;""",
					[bicluster[3], phenotype_query]
				)
				phenotype_result = c.fetchone()
				phenotype_result = [] if phenotype_result == None else phenotype_result

			# filter hallmarks using advanced search parameters
			if (
				has_hallmark_search(ordered_hallmarks, hallmarks)
				and (
					(
						phenotype_query != None
						and len(phenotype_result) > 0
					)
					or phenotype_query == None
				)
			):
				biclusters.append({
					"hallmarks": [(i, HALLMARK_TO_ICON[i]) for i in ordered_hallmarks],
					"data": bicluster,
					"phenotype": phenotype_result,
				})
	db.close()
	return render_template('search.html', gene=symbol, mutations=mutations, regulators=regulators, biclusters=biclusters)

@search_page.route('/advanced-search')
def advanced_search():
	phenotype = request.args.get('phenotype')
	phenotype = phenotype if phenotype else None
	hallmarks_query = request.args.get('hallmark')
	hallmarks_query = hallmarks_query if hallmarks_query else None

	# return to search index if bad query
	if not phenotype and not hallmarks_query:
		return redirect(url_for(".search"))

	hallmarks = decode_hallmark_query(hallmarks_query)

	db = dbconn()
	c = db.cursor()

	results = []
	# translate the hallmark UI ids to SQL ids
	hallmark_results = None
	searched_hallmarks = []
	if len(hallmarks) != 0:
		c.execute("SELECT name, id FROM hallmark;")
		results = {hallmark[0]: hallmark[1] for hallmark in c.fetchall()}
		sql_hallmarks = []
		for hallmark_id in hallmarks:
			searched_hallmarks.append((HALLMARKS[hallmark_id - 1], HALLMARK_TO_ICON[HALLMARKS[hallmark_id - 1]]))
			sql_hallmarks.append(results[HALLMARKS[hallmark_id - 1]])

		placeholders = ", ".join(["%s"] * len(sql_hallmarks))
		c.execute(
			"""SELECT b.id, b.name, b.var_exp_fpc, b.var_exp_fpc_p_value, b.survival, b.survival_p_value
			FROM bicluster b
			JOIN bic_hal bh ON bh.bicluster_id = b.id
			JOIN hallmark h ON bh.hallmark_id = h.id
			WHERE h.id IN (%s);""" % placeholders,
			sql_hallmarks
		)
		hallmark_results = c.fetchall()
	
	phenotype_results = None
	searched_phenotype = phenotype
	if phenotype:
		c.execute(
			"""SELECT b.id, b.name, b.var_exp_fpc, b.var_exp_fpc_p_value, b.survival, b.survival_p_value, bps.r_value, bps.p_value
			FROM bicluster b
			JOIN bicluster_phenotype_significance bps ON bps.bicluster_id = b.id
			JOIN phenotype p ON bps.phenotype_id = p.id
			WHERE p.long_name = %s AND bps.p_value < 0.05;""",
			[phenotype]
		)
		phenotype_results = c.fetchall()

	# if we have results for both, then intersect the two based off of first key, taking the datum that has the most information available
	results = []
	if phenotype_results and hallmark_results:
		hallmark_results_by_id = {i[0]: i for i in hallmark_results}
		
		for result in phenotype_results:
			if result[0] in hallmark_results_by_id:
				results.append(result)
	elif phenotype_results:
		results = phenotype_results
	elif hallmark_results:
		results = hallmark_results

	# construct the results
	biclusters = []
	for bicluster in results:
		c.execute(
			"""SELECT h.name
			FROM bic_hal bh
			JOIN hallmark h ON bh.hallmark_id = h.id
			WHERE bh.bicluster_id=%s;""",
			[bicluster[0]]
		)
		bicluster_hallmarks = [i[0] for i in c.fetchall()]
		bicluster_hallmarks = [(i, HALLMARK_TO_ICON[i]) for i in bicluster_hallmarks]
		biclusters.append({
			"hallmarks": bicluster_hallmarks,
			"data": bicluster,
		})

	db.close()
	
	return render_template(
		"advanced-search.html",
		biclusters=biclusters,
		searched_hallmarks=searched_hallmarks, searched_phenotype=searched_phenotype
	)