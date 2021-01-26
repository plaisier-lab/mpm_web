import json
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
		"""SELECT g.symbol
		FROM grna_jacks_result gjr
		JOIN grna ON grna.id = gjr.grna_id
		JOIN gene g ON g.id = grna.gene_id
		JOIN cell_line cl ON cl.id = gjr.cell_line_id
		WHERE cl.name="MESO1" AND g.id IN (
			SELECT bg.gene_id
			FROM bic_gene bg
			JOIN bicluster b ON bg.bicluster_id = b.id
			WHERE b.id = %s
		)
		GROUP BY g.id;""",
		[bicluster_id]
	)
	genes = [gene[0] for gene in c.fetchall()]

	db.close()

	return genes