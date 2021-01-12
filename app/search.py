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
		return redirect(url_for('.gene', symbol=symbol) + query_statement)
	else:
		return redirect(url_for('.mirna', symbol=symbol) + query_statement)


def __get_muts(c, gene_pk, symbol, hallmark_search=None):
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
				biclusters = c.fetchall()[0]
				c.execute(
					"SELECT hm.name FROM hallmark hm join bic_hal bh on hm.id=bh.hallmark_id  WHERE bh.bicluster_id=%s",
					[cf1[4]]
				)
				ordered_hallmarks = [i[0] for i in c.fetchall()]
				
				if has_hallmark_search(ordered_hallmarks, hallmark_search):
					hallmark = [(i, HALLMARK_TO_ICON[i]) for i in ordered_hallmarks]
					muts['data'].append([symbol, g1, biclusters[0], biclusters[1], biclusters[2], hallmark])

					if biclusters[0] not in muts['biclusters']:
						muts['biclusters'].append(biclusters[0])
	return muts


def __get_regulators(c, symbol, bme_type, hallmark_search=None):
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
			biclusters = c.fetchall()[0]
			c.execute("SELECT hm.name FROM hallmark hm join bic_hal bh on hm.id=bh.hallmark_id WHERE bh.bicluster_id=%s", [reg[1]])
			ordered_hallmarks = [i[0] for i in c.fetchall()]

			# filter hallmarks using advanced search parameters
			if has_hallmark_search(ordered_hallmarks, hallmark_search):
				hallmark = [(i, HALLMARK_TO_ICON[i]) for i in ordered_hallmarks]
				regs['data'].append([symbol, action, biclusters[0], biclusters[1], biclusters[2], hallmark])
				regs['biclusters'] = regs['biclusters'] + 1
	return regs

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

@search_page.route('/mirna')
@search_page.route('/mirna/<symbol>')
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


@search_page.route('/gene')
@search_page.route('/gene/<symbol>')
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
	muts = __get_muts(c, gene_pk, symbol, hallmark_search=hallmarks)

	c.execute("""SELECT * FROM tf_regulator WHERE tf_id=%s""", [gene_pk])
	regs = __get_regulators(c, symbol, "gene", hallmark_search=hallmarks)

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
				h1 = [(i, HALLMARK_TO_ICON[i]) for i in tmp1]
				bics['data'].append([bic1[4], bic1[5], bic1[6], bic1[7], bic1[8], h1])
				bics['biclusters'] = bics['biclusters'] + 1
	db.close()
	return render_template('search.html', gene=symbol, muts=muts, regs=regs, bics=bics)