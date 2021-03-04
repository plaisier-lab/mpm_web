##########################################################
## OncoMerge:  populateDatabase.py                      ##
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


'''
SELECT ge.value, cf.somatic_mutation_id, tf.id, b.id FROM gene_expression ge
	JOIN gene g ON ge.gene_id = g.id
	JOIN tf_regulator tf ON tf.tf_id = g.id
	JOIN causal_flow cf ON cf.regulator_id = tf.id AND cf.regulator_type = "tf"
	JOIN bicluster b ON cf.bicluster_id = b.id
	JOIN bueno_deep_filter bd ON bd.somatic_mutation_id = cf.somatic_mutation_id AND bd.patient_id = ge.patient_id
	WHERE cf.somatic_mutation_id = 8 AND tf.id = 3127 AND b.id = 355 AND bd.value = 1;

SELECT e.value FROM eigengene e
	JOIN bicluster b ON e.bicluster_id = b.id
	JOIN causal_flow cf ON cf.bicluster_id = b.id
	JOIN bueno_deep_filter bd ON bd.somatic_mutation_id = cf.somatic_mutation_id AND bd.patient_id = e.patient_id
	WHERE cf.somatic_mutation_id = 8 AND cf.regulator_id = 3127 AND cf.regulator_type = "tf" AND b.id = 355 AND bd.value = 1;
'''

import re
import cPickle
import os
from copy import deepcopy
from decimal import Decimal
import pandas as pd
import json
import time
import math

from flask import Flask
from flask_sqlalchemy import SQLAlchemy

from scipy import stats

from helpers import miRNAInDict, compareMiRNANames, hallmark_name_to_csv_key, parse_mutation_name, parse_locus_name, get_bicluster_id, parse_tf_matrix, parse_miRNA_array, isfloat
import logger

from models import *

#####################################
## Helpers                         ##
#####################################

def miRNAInDict(miRNA, dict1):
	retMe = []
	for i in dict1.keys():
		if compareMiRNANames(miRNA, i):
			retMe.append(miRNAIDs[i])
	if len(retMe) == 0:
		miRNA_shorter = '-'.join(miRNA.split('-')[0:3])
		for i in dict1.keys():
			if compareMiRNANames(miRNA_shorter, i):
				retMe.append(miRNAIDs[i])
	return retMe

def compareMiRNANames(a, b):
	if a == b:
		return 1
	if len(a) < len(b):
		if a[-3:] == '-3p':
			re1 = re.compile(a+'[a-oq-z]?(-\d)?-3p$')
		else:
			re1 = re.compile(a+'[a-oq-z]?(-\d)?-5p$')
		if re1.match(b):
			return 1
	else:
		if b[-3:] == '-3p':
			re1 = re.compile(b+'[a-oq-z]?(-\d)?-3p$')
		else:
			re1 = re.compile(b+'[a-oq-z]?(-\d)?-5p$')
		if re1.match(a):
			return 1
	return 0

def hallmark_name_to_csv_key(hallmark):
	split = hallmark.split(" ")
	output = ""
	for element in split:
		output += "{}{}".format(element[0].upper(), element[1:])
	return output

def parse_mutation_name(input):
	regex = re.compile(r"(X|x|)([0-9]+)_(.+)")
	match = regex.match(input)

	if match:
		return (int(match.group(2)), match.group(3))
	else:
		return None

def parse_locus_name(input):
	regex = re.compile(r"(X|x|)([^\s]+)_([A-Za-z]+)")
	match = regex.match(input)

	if match:
		return (match.group(2), match.group(3))
	else:
		return None

def get_bicluster_id(input):
	regex = re.compile(r"([a-z_]+)_([0-9]+)")
	match = regex.match(input)

	if match:
		return int(match.group(2))
	else:
		return None


def parse_tf_matrix(bicluster, input):
	if input == "NA":
		return None
	
	output = [bicluster]
	rows = input.split(" ")
	for col in rows:
		cols = col.split(":")
		entrez_id = int(cols[0])
		r_value = Decimal(cols[1])
		p_value = Decimal(cols[2])

		output.append((entrez_id, r_value, p_value))
	return output

def parse_miRNA_array(bicluster, input):
	if input == "NA":
		return None
	
	output = [bicluster]
	output = output + input.split(" ")
	return output

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

#####################################
## Load up conversion dictionaries ##
#####################################

start_time = time.time()

# Genes
symbol2entrez = {}
entrez2symbol = {}
print("gene2entrezId.csv")
with open('data/gene2entrezId.csv', 'r') as inFile:
	while 1:
		inLine = inFile.readline()
		if not inLine:
			break
		splitUp = inLine.strip().split(',')
		symbol2entrez[splitUp[0]] = splitUp[1]
		entrez2symbol[splitUp[1]] = splitUp[0]

# miRNAs
miRNAIDs = {}
miRNAIDs_rev = {}
mature_sequence_ids = []
print("hsa.mature.fa")
with open('data/hsa.mature.fa', 'r') as inFile:
	while 1:
		inLine = inFile.readline()
		if not inLine:
			break
		splitUp = inLine.split(' ')
		mature_sequence_ids.append(splitUp[1])
		if not splitUp[1] in miRNAIDs_rev:
			miRNAIDs_rev[splitUp[1]] = splitUp[0].lower()
		if not splitUp[0].lower() in miRNAIDs:
			miRNAIDs[splitUp[0].lower()] = splitUp[1]
		else:
			print('Uh oh!', splitUp)

####################################################
## Python script to populate the mpm_web database ##
####################################################

# Load up Genentech Gene Expression dataset

print("mesothelioma_norm.csv...")
# Read in gene expression matrix (Genentech)
gexp1 = pd.read_csv('data/mesothelioma_norm.csv', header=0, index_col=0, sep=',')

# Load up patient information (Genentech)
patients = {}
patient_check = re.compile('M?[0-9]+PT')
for patient1 in gexp1.columns:
	if patient_check.match(patient1):
		p1 = Patient(name=patient1)
		db.session.add(p1)
		patients[patient1] = p1

# Load up gene information (Genentech)
genes = {}
for gene1 in gexp1.index:
	if str(gene1) in entrez2symbol:
		g1 = Gene(symbol=entrez2symbol[str(gene1)], entrez=int(gene1))
		db.session.add(g1)
		genes[int(gene1)] = g1

db.session.commit()

# Add new expresssion dataset
ed1 = ExpDataset(name='Genentech_MPM', type='gexp')
db.session.add(ed1)
db.session.commit()

# Load up gene expression data
gene_expression_count = 0
for patient1 in gexp1.columns:
	for gene1 in gexp1.index:
		if gene1 in genes and patient_check.match(patient1):
			ge1 = GeneExpression(
				exp_dataset_id=ed1.id, patient_id=patients[patient1].id, gene_id=genes[int(gene1)].id, value=float(gexp1[patient1][gene1]))
			db.session.add(ge1)
			gene_expression_count = gene_expression_count + 1
		
		if gene_expression_count > 200000:
			print("committing {} gene expressions...".format(gene_expression_count))
			db.session.commit()
			gene_expression_count = 0
db.session.commit()

# Add hallmarks of Cancer
hallmark_names = ['Sustained angiogenesis', 'Insensitivity to antigrowth signals', 'Evading apoptosis', 'Limitless replicative potential', 'Evading immune detection',
	'Tissue invasion and metastasis', 'Self sufficiency in growth signals', 'Tumor promoting inflammation', 'Reprogramming energy metabolism', 'Genome instability and mutation']
hallmarks = {}
for hn1 in hallmark_names:
	hm1 = Hallmark(name=hn1)
	db.session.add(hm1)
	hallmarks[hn1] = hm1

# Add gene ontology and GO to gene mappings
# GO_BP
geneOntology = {}

print("go.obo...")
with open('data/go.obo', 'r') as inFile:
	while 1:
		line = inFile.readline()
		if not line:
			break

		if line.strip() == '[Term]':
			id1 = inFile.readline().strip().split('id: ')[1]
			name = inFile.readline().strip().split('name: ')[1]
			category = inFile.readline().strip().split('namespace: ')[1]
			if category == 'biological_process':
				go1 = GoBp(go_id=id1, name=name)
				geneOntology[id1] = go1
				db.session.add(go1)

		if line.strip()[0:6] == 'alt_id' and category == 'biological_process':
			id1 = line.strip().split('alt_id: ')[1]
			go1 = GoBp(go_id=id1, name=name)
			geneOntology[id1] = go1
			db.session.add(go1)

db.session.commit()

# Gene -> GO_BP
print("gene2go.hsa...")
with open('data/gene2go.hsa', 'r') as inFile:
	while 1:
		line = inFile.readline()
		if not line:
			break
		splitUp = line.strip().split('\t')
		if int(splitUp[1]) in genes and splitUp[2] in geneOntology:
			gg1 = GoGene(
				go_bp_id=geneOntology[splitUp[2]].id, gene_id=genes[int(splitUp[1])].id)
			db.session.add(gg1)

db.session.commit()

# miRNAs
print("miRNAs...")
miRNAs = {}
for miR1 in miRNAIDs:
	m1 = Mirna(mature_seq_id=miRNAIDs[miR1], name=miR1)
	db.session.add(m1)
	miRNAs[miRNAIDs[miR1]] = m1

db.session.commit()

# Prior information
# Load DisGeNET
print("all_gene_disease_pmid_associations.csv...")
mesos = ['Mesothelioma','Malignant mesothelioma','Malignant Pleural Mesothelioma',
	'Mesothelioma (malignant, clinical disorder) (disorder)','Malignant Mesothelioma of Peritoneum','Pleural Mesothelioma','Mesothelioma, biphasic, malignant','Sarcomatoid Mesothelioma']
with open('data/all_gene_disease_pmid_associations.tsv','r') as inFile:
	inFile.readline() # get rid of header
	while 1:
		inLine = inFile.readline()
		if not inLine:
			break
		splitUp = inLine.strip().split('\t')
		if splitUp[5] in mesos and int(splitUp[0]) in genes:
			pmid = splitUp[13]
			if not pmid=='':
				pmid = int(pmid)
			else:
				pmid = 0

			gp1 = GenePrior(gene_id=genes[int(
				splitUp[0])].id, source='DisGeNET', description=splitUp[5], pmid=pmid)
			db.session.add(gp1)

# Load HMDD (miRNA_Name,PMID,Description)
print("hmdd_mpm.csv...")
with open('data/hmdd_mpm.csv','r') as inFile:
   inFile.readline() # Get rid of header
   while 1:
	   inLine = inFile.readline()
	   if not inLine:
		   break
	   splitUp = inLine.strip().split(',')
	   miR1 = miRNAInDict(splitUp[0],miRNAIDs)
	   # print(splitUp[0],miR1)
	   for m1 in miR1:
		   mp1 = MirnaPrior(
			   mirna_id=miRNAs[m1].id, source='hmdd', description=splitUp[2], pmid=int(splitUp[1]))
		   db.session.add(mp1)

# Load dbDEMC
print("dbDEMC.csv...")
with open('data/dbDEMC.csv','r') as inFile:
   inFile.readline() # Get rid of header
   while 1:
	   inLine = inFile.readline()
	   if not inLine:
		   break
	   splitUp = inLine.strip().split(',')
	   miR1 = miRNAInDict(splitUp[0].lower(),miRNAIDs)
	   # print(splitUp[0],miR1)
	   for m1 in miR1:
		   mp1 = MirnaPrior(
			   mirna_id=miRNAs[m1].id, source='dbDEMC', description=splitUp[1], pmid=int(splitUp[2]))
		   db.session.add(mp1)

db.session.commit()

# load up TfTargets
print("humanTFs_All.csv...")
with open("data/humanTFs_All.csv") as in_file:
	in_file.readline() # consume headers

	motif_to_entrez = {}

	line = in_file.readline()
	while line:
		split = line.split(",")
		motif_to_entrez[split[0].lower()] = int(split[2])
		line = in_file.readline()

	in_file.close()

	# open up json file and start figuring out various TFTargets
	json_file = open("data/tfbsDb_plus_and_minus_5000_entrez.json")
	json_data = json.loads(json_file.readline())
	json_file.close()

	entries = {}

	not_found = []
	entries_made = 0
	entries_max = 200000
	for motif in json_data:
		if motif.lower() in motif_to_entrez:
			# entrez id associated with the motif name
			tf_entrez = motif_to_entrez[motif.lower()]
			# array of all the gene entrez ids associated with the motif
			gene_entrez = json_data[motif]

			for entrez in gene_entrez:
				entrez = int(entrez)

				if tf_entrez not in genes and tf_entrez not in not_found:
					print("tftarget warning: tf_entrez {} not found".format(tf_entrez))
					not_found.append(tf_entrez)

				if entrez not in genes and entrez not in not_found:
					print("tftarget warning: gene_entrez {} not found".format(entrez))
					not_found.append(entrez)

				if tf_entrez in genes and entrez in genes:
					# add an element to our dictionary
					if tf_entrez not in entries:
						entries[tf_entrez] = []

					# make sure we don't add any duplicates
					if entrez not in entries[tf_entrez]:
						db.session.add(TfTarget( \
							tf_id = genes[tf_entrez].id, \
							gene_id = genes[entrez].id \
						))
						
						entries[tf_entrez].append(entrez)
						entries_made = entries_made + 1
			
			# commit entries in 200k batches, so we don't eat up all memory
			if entries_made > entries_max:
				print("commiting {} tftarget entries...".format(entries_made))
				db.session.commit()
				entries_made = 0
	
	# final commit
	print("committing {} tftarget entries, final commit...".format(entries_made))
	db.session.commit()

# open up json file and start figuring out various miRNATargets
def handle_miRNA_json(file_name, source_name):
	json_file = open("data/{}".format(file_name))
	json_data = json.loads(json_file.readline())
	json_file.close()

	entries = {}

	not_found = []
	entries_made = 0
	entries_max = 200000
	for mimat in json_data:
		if mimat in miRNAs:
			# array of all the gene entrez ids associated with the mimat
			gene_entrez = json_data[mimat]

			for entrez in gene_entrez:
				entrez = int(entrez)

				if entrez not in genes and entrez not in not_found:
					print("miRNA target warning ({}): gene_entrez {} not found".format(source_name, entrez))
					not_found.append(entrez)

				if entrez in genes:
					# add an element to our dictionary
					if mimat not in entries:
						entries[mimat] = []

					# make sure we don't add any duplicates
					if entrez not in entries[mimat]:
						db.session.add(MirnaTarget( \
							mirna_id = miRNAs[mimat].id, \
							gene_id = genes[entrez].id, \
							source = source_name, \
						))
						
						entries[mimat].append(entrez)
						entries_made = entries_made + 1
			
			# commit entries in 200k batches, so we don't eat up all memory
			if entries_made > entries_max:
				print("commiting {} miRNA target entries ({})...".format(entries_made, source_name))
				db.session.commit()
				entries_made = 0
	
	# final commit
	print("committing {} miRNA target entries ({}), final commit...".format(entries_made, source_name))
	db.session.commit()

handle_miRNA_json("targetScan_miRNA_sets_entrez_hsa.json", "Target Scan")
handle_miRNA_json("pita_miRNA_sets_entrez_hsa.json", "Pita")

sif_mutation_to_regulators = {} # maps a mutation to an array of regualtors
sif_regulator_to_mutations = {} # maps a regulator to an array of mutations
sif_regulator_to_biclusters = {} # maps a regulator to an array of biclusters
sif_bicluster_to_regulators = {} # maps a bicluster to an array of regulators
def interpret_sif(filename):
	file = open(filename)
	line = file.readline()
	while line:
		split = line.split(" ")
		first = split[0].strip()
		command = split[1].strip()
		last = split[2].strip()

		line = file.readline()

		if command == "g2r": # mutation to regulator
			# first is a mutation
			# last is a regulator
			
			if last not in symbol2entrez and last.lower() not in miRNAIDs:
				continue

			if last.lower() in miRNAIDs:
				last = miRNAIDs[last.lower()] # translate readable name to MIMAT
			else:
				last = int(symbol2entrez[last])
			
			if first not in sif_mutation_to_regulators:
				sif_mutation_to_regulators[first] = []
			
			if last not in sif_regulator_to_mutations:
				sif_regulator_to_mutations[last] = []
			
			sif_mutation_to_regulators[first].append(last)
			sif_regulator_to_mutations[last].append(first)
		elif command == "r2b": # regulator to bicluster
			# first is a regulator
			# last is a bicluster
			
			if first not in symbol2entrez and first.lower() not in miRNAIDs:
				continue

			if "mir" in first.lower() and first.lower() not in miRNAIDs:
				print("failed to find mirna " + first + " in the dictionary")
				continue
			elif "mir" in first.lower():
				print("found " + first.lower() + " as " + miRNAIDs[first.lower()])
			
			if first.lower() in miRNAIDs:
				first = miRNAIDs[first.lower()] # translate readable name to MIMAT
			else:
				first = int(symbol2entrez[first])
			
			if first not in sif_regulator_to_biclusters:
				sif_regulator_to_biclusters[first] = []
			
			if last not in sif_bicluster_to_regulators:
				sif_bicluster_to_regulators[last] = []
			
			sif_regulator_to_biclusters[first].append(last)
			sif_bicluster_to_regulators[last].append(first)
	
	file.close()

somatic_mutations = {}
locus_map = {}
locus_array = []
'''
this might get the causal_flow structure

SELECT CONCAT_WS(
		'',
		CONCAT(l.locus_name, '_', l.mutation_type),
		CONCAT(g2.entrez, '_', sm.mutation_type)
	) AS mutation,
	g.symbol AS regulator,
	b.name AS bicluster
FROM causal_flow cf
	JOIN bicluster b ON cf.bicluster_id=b.id
	JOIN tf_regulator tf ON cf.regulator_id=tf.id
	JOIN gene g ON tf.tf_id=g.id 
	JOIN somatic_mutation sm ON cf.somatic_mutation_id=sm.id
	LEFT JOIN gene g2 ON sm.ext_id=g2.id
	LEFT JOIN locus l ON sm.locus_id=l.id;
'''
def interpret_causality_summary(filename, bicluster_prefix):
	file = open(filename)
	header = file.readline() # absorb the header

	mimat_regex = re.compile(r'^MIMAT')
	bicluster_number_regex = re.compile(r'[0-9]+')
	bicluster_prefix_regex = re.compile(r'^[A-Za-z_]+(?=_)')

	for line in file:
		split = [line.strip() for line in line.split(",")]

		# mutations
		mutation = parse_mutation_name(split[0]) # parse TF mutation
		is_locus = False
		if mutation == None: # parse locus
			mutation = parse_locus_name(split[0])
			sif_mutation_name = mutation
			is_locus = True
		else:
			sif_mutation_name = "{}_{}".format(entrez2symbol[str(mutation[0])], mutation[1])

		mutation_name = "{}_{}".format(mutation[0], mutation[1])

		# regulators
		entrez_or_mimat = split[1]
		regulator = ""
		if mimat_regex.match(entrez_or_mimat):
			regulator = entrez_or_mimat # MIMAT#######
		else:
			regulator = int(split[1]) # entrez id

		# biclusters		
		bicluster = int(bicluster_number_regex.findall(split[2])[0]) # digits at end of bicluster name
		# some files have bicluster prefixes, so we need to account for that
		if bicluster_prefix_regex.match(split[2]):
			bicluster_prefix = bicluster_prefix_regex.match(split[2]).group(0)
		bicluster_name = "{}_{}".format(bicluster_prefix, bicluster)
		
		leo_nb_atob = float(split[3])
		mlogp_m_atob = float(split[5])

		if (
			( # check mutation -> regulator edge
				sif_mutation_name in sif_mutation_to_regulators
				and regulator in sif_mutation_to_regulators[sif_mutation_name]
			)
			and ( # check regulator -> bicluster
				regulator in sif_regulator_to_biclusters
				and bicluster_name in sif_regulator_to_biclusters[regulator]
			)
		):
			bicluster_id = biclusters[bicluster_name].id
			regulator_id = ""
			regulator_type = ""

			# determine regulator type
			if regulator in genes:
				regulator_key = (bicluster_id, genes[regulator].id)

				if regulator_key not in tf_regulators_dict:
					continue

				regulator_id = tf_regulators_dict[regulator_key].id
				regulator_type = "tf"
			elif regulator in miRNAs:
				regulator_key = (bicluster_id, miRNAs[regulator].id)

				if regulator_key not in mirna_regulators_dict:
					continue

				regulator_id = mirna_regulators_dict[regulator_key].id
				regulator_type = "mirna"
			
			# create somatic mutation if we don't have one
			if mutation[0] in genes and is_locus == False and mutation_name not in somatic_mutations:
				somatic_mutations[mutation_name] = SomaticMutation(
					ext_id = genes[mutation[0]].id,
					mutation_type = mutation[1], # PAM, LoF, etc
					mutation_name = mutation_name,
				)
				db.session.add(somatic_mutations[mutation_name])
			# create locus and somatic mmutation if we don't have one
			elif is_locus == True and mutation_name not in locus_map and mutation_name not in somatic_mutations:
				locus_map[mutation_name] = Locus(
					locus_name = mutation[0],
					mutation_type = mutation[1],
				)
				db.session.add(locus_map[mutation_name])

				somatic_mutations[mutation_name] = SomaticMutation(
					locus = locus_map[mutation_name],
					mutation_name = mutation_name,
				)
				db.session.add(somatic_mutations[mutation_name])
			elif mutation_name not in somatic_mutations: # quit if we can't find a somatic mutation
				continue
			
			# create causal flow
			causal_flow = CausalFlow(
				somatic_mutation = somatic_mutations[mutation_name],
				regulator_id = regulator_id,
				regulator_type = regulator_type,
				bicluster_id = bicluster_id,
				leo_nb_atob = leo_nb_atob,
				mlogp_m_atob = mlogp_m_atob,
			)
			db.session.add(causal_flow)
	db.session.commit()

# parse phenotypes
phenotypes = {}
with open('data/phenotypes_meso_noFilter.csv') as in_file:
	# read header
	header = in_file.readline().split(",")
	for i in range(1, len(header)):
		name = header[i].strip()
		phenotypes[name] = Phenotype( \
			name = name,
			long_name = name,
		)
		db.session.add(phenotypes[name])
	
	db.session.commit()

	# for every other line, read information out and assoiate it with a header name
	line = in_file.readline()
	while line:
		parsed = line.split(",")
		current_patient = None
		for i in range(0, len(parsed)):
			if i != 0 and current_patient != None:
				phenotype_name = header[i].strip()
				phenotype_value = parsed[i].strip()
				if isfloat(phenotype_value): # handle floats
					db.session.add(PhenoDatum( \
						phenotype_id = phenotypes[phenotype_name].id,
						patient_id = current_patient.id,
						phenotype_value = float(phenotype_value),
					))
				else: # handle strings
					db.session.add(PhenoDatum( \
						phenotype_id = phenotypes[phenotype_name].id,
						patient_id = current_patient.id,
						phenotype_string = phenotype_value,
					))
			else:
				patient_name = parsed[i].strip()
				if patient_name in patients:
					current_patient = patients[parsed[i]]
				else:
					current_patient = None

		line = in_file.readline()
	
	db.session.commit()

def read_eigengenes(filename, prefix):
	eigengenes = pd.read_csv(filename, header=0, index_col=0)
	eigengene_count = 0
	for patient in eigengenes.columns:
		old_patient = patient
		if patient[0].lower() == "x":
			patient = patient[1:]
		
		if patient_check.match(patient) == None:
			continue

		for bicluster_number in eigengenes.index:
			bicluster_name = "{}_{}".format(prefix, bicluster_number)
			if bicluster_name in biclusters and patient in patients:
				eigengene = Eigengene(
					patient_id = patients[patient].id,
					bicluster_id = biclusters[bicluster_name].id,
					value = float(eigengenes[old_patient][bicluster_number]))
				db.session.add(eigengene)
				eigengene_count = eigengene_count + 1
			
			if eigengene_count > 200000:
				print("committing {} eigengenes...".format(eigengene_count))
				db.session.commit()
				eigengene_count = 0
	db.session.commit()

mimat_r_values = {}
mimat_p_values = {}

def read_miRNA_p_r_values(r_values_file, p_values_file):
	r_values_file = open(r_values_file, "r")
	p_values_file = open(p_values_file, "r")

	# figure out what mimats we're dealing with
	r_values_mimat = [mimat.strip()[1:-1] for mimat in r_values_file.readline().split(",")[1:]]
	p_values_mimat = [mimat.strip()[1:-1] for mimat in p_values_file.readline().split(",")[1:]]

	# read in r-values
	for line in r_values_file:
		split = line.split(",")
		bicluster = split[0][1:-1]

		index = 0
		for value in split[1:]:
			mimat = r_values_mimat[index]
			value = float(value)
			mimat_r_values[(bicluster, mimat)] = value
			index = index + 1
	
	# read in p-values
	for line in p_values_file:
		split = line.split(",")
		bicluster = split[0][1:-1]

		index = 0
		for value in split[1:]:
			mimat = p_values_mimat[index]
			value = float(value)
			mimat_p_values[(bicluster, mimat)] = value
			index = index + 1

	r_values_file.close()
	p_values_file.close()

read_miRNA_p_r_values(
	"./data/miRNA_cor_mesoTCGA/eig_miR_cor_12112020.csv",
	"./data/miRNA_cor_mesoTCGA/eig_miR_pv_12112020.csv",
)

biclusters = {}
tf_regulators_dict = {}
mirna_regulators_dict = {}
with open('data/postProcessed_clustersOfBiclusters_CNA_CNVkit_01212021.csv') as in_file:
	# we need to cache the values so we can create the classes in order. we need to do this so the database can automatically complete the relationships defined in models.py
	bic_gene_values = []
	bic_go_values = []
	bic_pat_values = []
	bic_halmark_values = []
	replication_values = []
	
	header = in_file.readline()
	header_dict = {}
	# make a dictionary of header -> column id names
	header_split = header.strip().split(',')
	index = 0
	for header_name in header_split:
		if header_name not in header_dict:
			header_dict[header_name] = index
		else: # handle duplicate names
			suffix = 2
			while "{}{}".format(header_name, suffix) in header_dict:
				suffix = suffix + 1
			
			# for duplicate names, keep adding 1 to a suffix as long as the name is taken. once we find an untaken suffix, we will stop.
			# for instance, a duplicate colmun Genes will become Genes2, and another duplicate of the same name will become Genes3
			header_dict["{}{}".format(header_name, suffix)] = index
		
		index = index + 1

	tf_regulators = []
	miRNA_regulators = []

	# now read the file in full
	line = in_file.readline()
	while line:
		split = line.strip().split(',')

		bicluster_name = split[header_dict["bicluster"]]

		# create a bicluster
		if bicluster_name not in biclusters:
			bicluster = Bicluster(
				name = bicluster_name,
				var_exp_fpc = split[header_dict["Var. Exp. First PC"]],
				var_exp_fpc_p_value = split[header_dict["Var. Exp. First PC Perm. P-Value"]],
				survival = split[header_dict["OS.covAgeSex"]],
				survival_p_value = split[header_dict["OS.covAgeSex.p"]]
			)
			db.session.add(bicluster)

			biclusters[bicluster_name] = bicluster
		
		# create a BicGene
		bic_gene_values.append((bicluster_name, split[header_dict["Genes2"]])) # fix, should be split by spaces

		# create a BicGo
		bic_go_values.append((bicluster_name, split[header_dict["GO_Term_BP"]]))

		# create a BicPat
		bic_pat_values.append((bicluster_name, split[header_dict["Conditions"]]))

		# create a Replication
		replication_values.append((bicluster_name, split[header_dict["mesoTCGA_var.exp"]], split[header_dict["mesoTCGA_pc1.perm.p"]], split[header_dict["mesoTCGA_OS.age.sex"]], split[header_dict["mesoTCGA_OS.age.sex.p"]]))

		# handle TF regulators
		meme_motif1_matches = parse_tf_matrix(bicluster_name, split[header_dict["Up.MEME Motif1 Correlated Matches"]])
		if meme_motif1_matches != None:
			tf_regulators.append(meme_motif1_matches)

		meme_motif2_matches = parse_tf_matrix(bicluster_name, split[header_dict["Up.MEME Motif2 Correlated Matches"]])
		if meme_motif2_matches != None:
			tf_regulators.append(meme_motif2_matches)

		weeder_motif1_matches = parse_tf_matrix(bicluster_name, split[header_dict["Up.WEEDER Motif1 Correlated Matches"]])
		if weeder_motif1_matches != None:
			tf_regulators.append(weeder_motif1_matches)

		weeder_motif2_matches = parse_tf_matrix(bicluster_name, split[header_dict["Up.WEEDER Motif2 Correlated Matches"]])
		if weeder_motif2_matches != None:
			tf_regulators.append(weeder_motif2_matches)

		tfbsdb_matches = parse_tf_matrix(bicluster_name, split[header_dict["TFBS_DB.Correlated Matches"]])
		if tfbsdb_matches != None:
			tf_regulators.append(tfbsdb_matches)

		# handle miRNA regulators
		weeder_motif1_miRNA = parse_miRNA_array(bicluster_name, split[header_dict["3pUTR.WEEDER Motif1 Matches"]])
		if weeder_motif1_miRNA != None:
			miRNA_regulators.append(weeder_motif1_miRNA)
		
		weeder_motif2_miRNA = parse_miRNA_array(bicluster_name, split[header_dict["3pUTR.WEEDER Motif2 Matches"]])
		if weeder_motif2_miRNA != None:
			miRNA_regulators.append(weeder_motif2_miRNA)
		
		pita_miRNA = parse_miRNA_array(bicluster_name, split[header_dict["3pUTR_pita.miRNAs"]])
		if pita_miRNA != None:
			miRNA_regulators.append(pita_miRNA)
		
		target_scan_miRNA = parse_miRNA_array(bicluster_name, split[header_dict["3pUTR_targetScan.miRNAs"]])
		if target_scan_miRNA != None:
			miRNA_regulators.append(target_scan_miRNA)
		
		# figure out all the different hallmarks (there's a ton of them)
		for hallmark in hallmark_names:
			csv_key = hallmark_name_to_csv_key(hallmark)
			key = header_dict['"{}"'.format(csv_key)] # for some reason i have to surround this in quote literals

			if split[key] != "NA":
				bic_halmark_values.append((bicluster_name, hallmark, float(split[key])))

		line = in_file.readline()
	
	db.session.commit()

	# create all the BicGenes now that we have loaded bicluster information
	for tuple_val in bic_gene_values:
		if tuple_val[0] not in biclusters:
			print("bic_gene warning: could not find bicluster name {}".format(tuple_val[0]))
		else:
			split = tuple_val[1].split(" ")
			for entrez in split:
				entrez = int(entrez)
				if entrez not in genes:
					print("bic_gene warning: could not find gene id {}".format(entrez))
				else:
					db.session.add(BicGene(
						bicluster_id = biclusters[tuple_val[0]].id,
						gene_id = genes[entrez].id
					))

		'''
		elif tuple_val[1] not in genes:
			print("bic_gene warning: could not find gene id {}".format(tuple_val[1])) # report: sometimes a gene just cannot be associated with a bicluster, which screws with results in the search. what to do?
		else:
			db.session.add(BicGene(
				bicluster_id = biclusters[tuple_val[0]].id,
				gene_id = genes[tuple_val[1]].id
			))
		'''
	
	# create all the BicGos
	for tuple_val in bic_go_values:
		if tuple_val[1] != "NA":
			split = tuple_val[1].split(";")
			for go_val in split:
				db.session.add(BicGo(
					bicluster_id = biclusters[tuple_val[0]].id,
					go_bp_id = geneOntology[go_val].id
				))
	
	# create all the BicPat's
	for tuple_val in bic_pat_values:
		split = tuple_val[1].split(" ")
		for patient in split:
			if patient_check.match(patient):
				db.session.add(BicPat(
					bicluster_id = biclusters[tuple_val[0]].id,
					patient_id = patients[patient].id
				))

	# create all the BicHals
	for tuple_val in bic_halmark_values:
		bicluster_name = tuple_val[0]
		hallmark_name = tuple_val[1]
		hallmark_value = tuple_val[2]

		if hallmark_value >= 0.8:
			db.session.add(BicHal(
				bicluster_id = biclusters[bicluster_name].id,
				hallmark_id = hallmarks[hallmark_name].id
			))
	
	# create all the Replications
	for tuple_val in replication_values:
		bicluster_name = tuple_val[0]
		var_exp_fpc = tuple_val[1]
		var_exp_fpc_p_value = tuple_val[2]
		survival = tuple_val[3]
		survival_p_value = tuple_val[4]

		db.session.add(Replication(
			bicluster_id = biclusters[bicluster_name].id,
			study = "MESO TCGA",
			var_exp_fpc = var_exp_fpc,
			var_exp_fpc_p_value = var_exp_fpc_p_value,
			survival = survival,
			survival_p_value = survival_p_value
		))
	
	# create all TF regulators
	for tf_regulator_array in tf_regulators:
		bicluster_name = tf_regulator_array[0]
		for index in range(1, len(tf_regulator_array)):
			tf_regulator = tf_regulator_array[index]
			key = (biclusters[bicluster_name].id, genes[tf_regulator[0]].id)
			if key not in tf_regulators_dict: # have to add this so we do not add duplicate tf regulator entries
				regulator_object = TfRegulator(
					bicluster_id = biclusters[bicluster_name].id,
					tf_id = genes[tf_regulator[0]].id,
					r_value = tf_regulator[1],
					p_value = tf_regulator[2],
				)
				key = (biclusters[bicluster_name].id, genes[tf_regulator[0]].id)
				tf_regulators_dict[key] = regulator_object
				db.session.add(regulator_object)
	
	# create all miRNA regulators
	for miRNA_regulator_array in miRNA_regulators:
		bicluster_name = miRNA_regulator_array[0]
		for index in range(1, len(miRNA_regulator_array)):
			mimat = miRNA_regulator_array[index]
			if mimat in miRNAs:
				miRNA_entry = miRNAs[mimat]
				key = (biclusters[bicluster_name].id, miRNA_entry.id)
				if key not in mirna_regulators_dict: # have to make sure we don't add duplicate mirna regulator entries
					p_value = None
					r_value = None
					if (bicluster_name, mimat) in mimat_p_values:
						p_value = mimat_p_values[(bicluster_name, mimat)]
						r_value = mimat_r_values[(bicluster_name, mimat)]
					
					regulator_object = MirnaRegulator(
						bicluster_id = biclusters[bicluster_name].id,
						mirna_id = miRNA_entry.id,
						p_value = p_value,
						r_value = r_value,
					)
					mirna_regulators_dict[key] = regulator_object
					db.session.add(regulator_object)
			else:
				print("missing MIMAT {}, can't link to bicluster table".format(mimat))

	print("committing post processed biclusters...")
	db.session.commit()
	print("committed")


# handle somatic mutations and causal flows
interpret_sif("./data/sifs/causalAndMechanistic_network_CNA_CNVkit_01_21_2021.sif")
interpret_causality_summary("./data/causality_CNA_final_8_13_2019/causalitySummary_pita_12112020.csv", "pita")
interpret_causality_summary("./data/causality_CNA_final_8_13_2019/causalitySummary_targetscan_12112020.csv", "targetscan")
interpret_causality_summary("./data/causality_CNA_final_8_13_2019/causalitySummary_tfbs_db_12112020.csv", "tfbs_db")
interpret_causality_summary("./data/causal_v9/summaryCausality_CNV_11_1_2020_0.3_0.05_cleanedUp_12112020.csv", None)

with open("./data/oncoMerged_MESO/oncoMerge_mergedMuts_12112020.csv") as file:
	header = file.readline().split(',')[1:]
	for line in file:
		somatic_mutation = line.split(',')[0]
		booleans = line.split(',')[1:]
		for index in range(0, len(booleans)):
			patient = header[index].strip()
			if somatic_mutation in somatic_mutations and patient in patients:
				db.session.add(BuenoDeepFilter(
					somatic_mutation_id = somatic_mutations[somatic_mutation].id,
					patient_id = patients[patient].id,
					value = booleans[index].strip() == '1',
				))
	db.session.commit()

# handle eigengenes
read_eigengenes("./data/eigengenes/biclusterEigengenes_pita.csv", "pita")
read_eigengenes("./data/eigengenes/biclusterEigengenes_targetscan.csv", "targetscan")
read_eigengenes("./data/eigengenes/biclusterEigengenes_tfbs_db.csv", "tfbs_db")

# find bicluster phenotype significances
for bicluster in biclusters.values():
	bicluster_name = bicluster.name
	for phenotype in phenotypes.values():
		values = db.engine.execute(
			"SELECT p.name, pd.phenotype_value FROM patient p "
			+ "JOIN pheno_data pd ON pd.patient_id=p.id "
			+ "JOIN phenotype pt ON pd.phenotype_id=pt.id "
			+ "WHERE pt.name = '%s';" % phenotype.name
		)
		value_ptmap = {patient: value for patient, value in values}

		result = db.engine.execute(
			"SELECT e.value, pt.name FROM eigengene e "
			+ "JOIN patient pt on e.patient_id = pt.id "
			+ "JOIN bicluster b ON e.bicluster_id = b.id "
			+ "WHERE b.name = '%s';" % bicluster_name
		)

		eigengene_values = []
		phenotype_values = []
		for eigengene_value, patient in result:
			if patient not in value_ptmap or value_ptmap[patient] == None:
				continue

			eigengene_values.append(eigengene_value)
			phenotype_values.append(value_ptmap[patient])

		if len(eigengene_values) > 0:
			found_stats = stats.pearsonr(eigengene_values, phenotype_values)

			db.session.add(
				BiclusterPhenotypeSignificance(
					bicluster_id = bicluster.id,
					phenotype_id = phenotype.id,
					r_value = found_stats[0],
					p_value = found_stats[1],
				)
			)

db.session.commit()

def interpret_jacks_results(subtypes_file, results_file, std_file, p_value_file):
	subtypes_file = open(subtypes_file, "r")
	results_file = open(results_file, "r")
	std_file = open(std_file, "r")
	p_value_file = open(p_value_file, "r")
	
	# read subtypes
	subtypes = {}
	for line in subtypes_file:
		line = [i.strip() for i in line.split("\t")]
		subtypes[line[0]] = line[1]

	first_line = results_file.readline().split("\t")
	cell_lines = []
	# create cell lines
	for index in range(1, len(first_line)):
		element = first_line[index].strip()
		cell_line = CellLine(
			name=element,
			subtype=subtypes[element],
		)
		cell_lines.append(cell_line)
		db.session.add(cell_line)
	
	db.session.commit() # commit cell lines

	results = {}

	# handle results in results file
	for line in results_file:
		line = [i.strip() for i in line.split("\t")]
		gene_name = line[0]

		if gene_name not in symbol2entrez or int(symbol2entrez[gene_name]) not in genes:
			print("could not find %s" % gene_name)
			continue

		# loop through the cell lines we found in the file
		for index in range(1, len(line)):
			value = float(line[index])
			cell_line = cell_lines[index - 1]

			if gene_name not in results:
				results[gene_name] = {}

			results[gene_name][cell_line] = GeneJACKSResult(
				gene_id=genes[int(symbol2entrez[gene_name])].id,
				cell_line_id=cell_line.id,
				score=value,
			)
	
	# handle the std file
	std_file.readline() # absorb header
	for line in std_file:
		line = [i.strip() for i in line.split("\t")]
		gene_name = line[0]

		if gene_name not in symbol2entrez or int(symbol2entrez[gene_name]) not in genes:
			continue

		# loop through the cell lines we found in the file
		for index in range(1, len(line)):
			value = float(line[index])
			cell_line = cell_lines[index - 1]
			results[gene_name][cell_line].std = value
	
	# handle the p-value file
	p_value_file.readline() # absorb header
	for line in p_value_file:
		line = [i.strip() for i in line.split("\t")]
		gene_name = line[0]

		if gene_name not in symbol2entrez or int(symbol2entrez[gene_name]) not in genes:
			continue

		# loop through the cell lines we found in the file
		for index in range(1, len(line)):
			value = float(line[index])
			cell_line = cell_lines[index - 1]
			results[gene_name][cell_line].p_value = value
	
	# commit the results
	for cell_line_dict in results.values():
		for data in cell_line_dict.values():
			db.session.add(data)

	db.session.commit()

	subtypes_file.close()
	results_file.close()
	std_file.close()
	p_value_file.close()

	return cell_lines

def interpret_jacks_grna_results(means_file, std_file, sequences_file, cell_lines=[]):
	means_file = open(means_file, "r")
	std_file = open(std_file, "r")
	sequences_file = open(sequences_file, "r")

	results = {}
	grna = {}

	# read through the means file
	means_file.readline() # absorb header
	for line in means_file:
		line = [i.strip() for i in line.split("\t")]
		grna_name = line[0]
		gene_name = line[1]

		if gene_name not in symbol2entrez or int(symbol2entrez[gene_name]) not in genes:
			continue

		# loop through the cell lines we found in the file
		for index in range(2, len(line)):
			value = float(line[index])
			cell_line = cell_lines[index - 2]

			if grna_name not in results:
				results[grna_name] = {}
			
			if grna_name not in grna:
				grna[grna_name] = GRNA(
					name=grna_name,
					gene_id=genes[int(symbol2entrez[gene_name])].id,
				)
			
			results[grna_name][cell_line] = GRNAJACKSResult(
				cell_line_id=cell_line.id,
				grna=grna[grna_name],
				mean=value,
			)
	
	# read through the std file
	std_file.readline() # absorb header
	for line in std_file:
		line = [i.strip() for i in line.split("\t")]
		grna_name = line[0]
		gene_name = line[1]

		if gene_name not in symbol2entrez or int(symbol2entrez[gene_name]) not in genes:
			continue
	
		# loop through the cell lines we found in the file
		for index in range(2, len(line)):
			value = float(line[index])
			cell_line = cell_lines[index - 2]

			if not math.isnan(value):
				results[grna_name][cell_line].std = value
	
	# read through the sequences file
	sequences_file.readline() # absorb header
	for line in sequences_file:
		split = [i.strip() for i in line.split("\t")]
		sequence = split[0]
		grna_name = split[1]

		if grna_name in grna:
			db.session.add(GRNASequence(
				sequence=sequence,
				grna=grna[grna_name],
			))
	
	# commit the results
	for cell_line_dict in results.values():
		for data in cell_line_dict.values():
			db.session.add(data)

	db.session.commit()

	means_file.close()
	std_file.close()
	sequences_file.close()

# we need something specific for this cell line
def interpret_meso1_results(results_file, sequences_file):
	results_file = open(results_file, "r")
	sequences_file = open(sequences_file, "r")
	
	# create the MESO1 cell line
	cell_line = CellLine(
		name="MESO1"
	)

	grnas = {}

	# read the results file
	results_file.readline() # absorb header
	for line in results_file:
		split = [i.strip() for i in line.split(",")]
		gene_name, grna_id = split[0].split("_")
		logfoldchange = float(split[7])
		p_value = float(split[9])

		if gene_name not in symbol2entrez or int(symbol2entrez[gene_name]) not in genes:
			print("could not find %s for MESO1" % gene_name)
			continue

		grnas[grna_id] = GRNA(
			name=grna_id,
			gene_id=genes[int(symbol2entrez[gene_name])].id
		)

		db.session.add(GRNAJACKSResult(
			grna=grnas[grna_id],
			mean=logfoldchange,
			p_value=p_value,
			cell_line=cell_line
		))
	
	# read the sequences file
	sequences_file.readline() # absorb header
	for line in sequences_file:
		split = [i.strip() for i in line.split(",")]
		grna_id = split[2]
		sequence = split[4]

		if grna_id not in grnas:
			continue

		db.session.add(GRNASequence(
			grna=grnas[grna_id],
			sequence=sequence
		))
	
	print("committing MESO1...")
	db.session.commit()

	results_file.close()
	sequences_file.close()

def interpret_achilles_common_essential(file):
	file = open(file, "r")

	file.readline() # skip header
	for line in file:
		gene_name, entrez = line.strip().split(" ")

		if gene_name not in symbol2entrez or int(symbol2entrez[gene_name]) not in genes:
			continue
		
		db.session.add(AchillesCommonEssential(
			gene_id=genes[int(symbol2entrez[gene_name])].id,
		))
	
	db.session.commit()

	file.close()

def interpret_achilles(cell_line_file, subtype_file, effect_file, dependency_file):
	cell_line_file = open(cell_line_file, "r")
	subtype_file = open(subtype_file, "r")
	effect_file = open(effect_file, "r")
	dependency_file = open(dependency_file, "r")

	cell_lines = {}
	gene_results = {}

	# read in cell lines
	cell_line_file.readline() # absorb header
	for line in cell_line_file:
		split = [i.strip() for i in line.split(",")]
		achilles_id = split[0] # achilles id is ACH-######
		name = split[1] # actual name
		cancer_type = split[23]

		if cancer_type == "mesothelioma":
			cell_lines[achilles_id] = CellLine(
				name=name,
			)
		
	# read in subtypes
	for line in subtype_file:
		split = [i.strip() for i in line.split("\t")]
		achilles_id = split[0]
		subtype = split[1]

		cell_lines[achilles_id].subtype = subtype

	# read in effect (essentiality score)
	print("reading achilles essentiallity scores...")
	header = effect_file.readline().split(",")
	found_genes = [int(gene.strip().split(" ")[1][1:-1]) for gene in header[1:]]
	for line in effect_file:
		# trick: to avoid premature splitting (which will have high cpu/memory usage), i know the
		# cell line IDs are fixed length. just gonna grab them with some array manipulation
		cell_line = line[0:10]

		# discard line if we couldn't find the cell line in the lookup we made earlier from the summary
		if cell_line not in cell_lines:
			continue
		
		split = line.split(",")
		for index in range(1, len(split)):
			entrez_id = found_genes[index - 1]
			if entrez_id not in genes: # skip genes that we don't know about
				continue

			data = float(split[index])

			if entrez_id not in gene_results:
				gene_results[entrez_id] = {}

			# add result based on data we found
			gene_results[entrez_id][cell_line] = GeneAchillesResult(
				gene_id=genes[entrez_id].id,
				score=data,
				cell_line=cell_lines[cell_line],
			)
		
	# read in dependency (p-value)
	print("reading achilles p-values...")
	header = dependency_file.readline().split(",")
	found_genes = [int(gene.strip().split(" ")[1][1:-1]) for gene in header[1:]]
	for line in dependency_file:
		# trick: to avoid premature splitting (which will have high cpu/memory usage), i know the
		# cell line IDs are fixed length. just gonna grab them with some array manipulation
		cell_line = line[0:10]

		# discard line if we couldn't find the cell line in the lookup we made earlier from the summary
		if cell_line not in cell_lines:
			continue
		
		split = line.split(",")
		for index in range(1, len(split)):
			entrez_id = found_genes[index - 1]
			if entrez_id not in genes: # skip genes that we don't know about
				continue

			data = float(split[index])

			# update p-value
			if entrez_id in gene_results:
				gene_results[entrez_id][cell_line].p_value = data
	
	# add the cell lines
	for cell_line in cell_lines.values():
		db.session.add(cell_line)

	# add the results
	for cell_line_dict in gene_results.values():
		for data in cell_line_dict.values():
			db.session.add(data)
	
	print("committing achilles results...")
	db.session.commit()

	cell_line_file.close()
	subtype_file.close()
	effect_file.close()
	dependency_file.close()

cell_lines = interpret_jacks_results(
	"./data/jacks/subtypes.txt",
	"./data/jacks/MPM_6_1_20_gene_JACKS_results.txt",
	"./data/jacks/MPM_6_1_20_gene_std_JACKS_results.txt",
	"./data/jacks/MPM_6_1_20_gene_pval_JACKS_results.txt"
)

interpret_jacks_grna_results(
	"./data/jacks/MPM_6_1_20_logfoldchange_means.txt",
	"./data/jacks/MPM_6_1_20_logfoldchange_std.txt",
	"./data/jacks/barcode-counts.txt",
	cell_lines=cell_lines
)

interpret_meso1_results(
	"./data/MESO1/MESO1_CRISPR_Results.csv",
	"./data/MESO1/MESO_sgRNA_seqs.csv"
)

interpret_achilles_common_essential(
	"./data/Achilles_common_essentials.csv"
)

interpret_achilles(
	"./data/achilles/sample_info.csv",
	"./data/achilles/subtype.txt",
	"./data/achilles/Achilles_gene_effect.csv",
	"./data/achilles/Achilles_gene_dependency.csv"
)

print("program took {} seconds to complete".format(time.time() - start_time))