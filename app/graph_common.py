"""
contains some functions that are used in a few different analyses
"""

import numpy as np
import itertools
import math
from scipy import stats

from database import dbconn
from constants import NUM_PARTS

def phyper(x, M, n, bigN, lower_tail=False):
	"""uses scipy.stats to compute hypergeometric overlap p-value"""
	return 1 - stats.hypergeom.cdf(float(x), float(M), float(n), float(bigN))

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