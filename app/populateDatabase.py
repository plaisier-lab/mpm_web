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

import re
import cPickle
import os
from copy import deepcopy
from decimal import Decimal
# from cMonkeyWrapper import cMonkeyWrapper
# from bicluster import bicluster
# from pssm import pssm
import pandas as pd
import json

from flask import Flask
from flask_sqlalchemy import SQLAlchemy

"""
###########################
## Open database session ##
###########################
app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+mysqlconnector://root:root@mpm_web_db_1:3306/mpm_web'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = 'FALSE'

db = SQLAlchemy(app)
"""

from models import *

########################################
### Functions to convert miRNA names ###
########################################


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


#####################################
## Load up conversion dictionaries ##
#####################################

# Genes
symbol2entrez = {}
entrez2symbol = {}
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

print("mesothelioma_norm.txt...")
# Read in gene expression matrix (Genentech)
gexp1 = pd.read_csv('data/mesothelioma_norm.txt',
                    header=0, index_col=0, sep=' ')

# Load up patient information (Genentech)
patients = {}
for patient1 in gexp1.columns:
    p1 = Patient(name=patient1)
    db.session.add(p1)
    patients[patient1] = p1

# Load up gene information (Genentech)
genes = {}
for gene1 in gexp1.index:
    if str(gene1) in entrez2symbol:
        g1 = Gene(symbol=entrez2symbol[str(gene1)], entrez=int(gene1))
        db.session.add(g1)
        genes[gene1] = g1

db.session.commit()

# Add new expresssion dataset
ed1 = ExpDataset(name='Genentech_MPM', type='gexp')
db.session.add(ed1)
gt_gexp = deepcopy(ed1)

# Load up gene expression data
'''
gene_expression_count = 0
for patient1 in gexp1.columns:
    for gene1 in gexp1.index:
        if gene1 in genes:
            ge1 = GeneExpression(
                exp_dataset_id=gt_gexp.id, patient_id=patients[patient1].id, gene_id=genes[gene1].id, value=float(gexp1[patient1][gene1]))
            db.session.add(ge1)
            gene_expression_count = gene_expression_count  +1
        
        if gene_expression_count > 200000:
            print("committing {} gene expressions...".format(gene_expression_count))
            db.session.commit()
            gene_expression_count = 0

db.session.commit()
'''

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

print("gene2go.hsa...")

# Gene -> GO_BP
with open('data/gene2go.hsa', 'r') as inFile:
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split('\t')
        if splitUp[1] in genes and splitUp[2] in geneOntology:
            gg1 = GoGene(
                go_bp_id=geneOntology[splitUp[2]].id, gene_id=genes[splitUp[1]].id)
            db.session.add(gg1)

db.session.commit()

# miRNAs
miRNAs = {}
for miR1 in miRNAIDs:
    m1 = Mirna(mature_seq_id=miRNAIDs[miR1], name=miR1)
    db.session.add(m1)
    miRNAs[miRNAIDs[miR1]] = m1

# Prior information

# Load DisGeNET
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
with open('data/hmdd_mpm.csv','r') as inFile:
   inFile.readline() # Get rid of header
   while 1:
       inLine = inFile.readline()
       if not inLine:
           break
       splitUp = inLine.strip().split(',')
       miR1 = miRNAInDict(splitUp[0],miRNAIDs)
       print(splitUp[0],miR1)
       for m1 in miR1:
           mp1 = MirnaPrior(
               mirna_id=miRNAs[m1].id, source='hmdd', description=splitUp[2], pmid=int(splitUp[1]))
           db.session.add(mp1)


# Load dbDEMC
with open('data/dbDEMC.csv','r') as inFile:
   inFile.readline() # Get rid of header
   while 1:
       inLine = inFile.readline()
       if not inLine:
           break
       splitUp = inLine.strip().split(',')
       miR1 = miRNAInDict(splitUp[0].lower(),miRNAIDs)
       print(splitUp[0],miR1)
       for m1 in miR1:
           mp1 = MirnaPrior(
               mirna_id=miRNAs[m1].id, source='dbDEMC', description=splitUp[1], pmid=int(splitUp[2]))
           db.session.add(mp1)

db.session.commit()


def hallmark_name_to_csv_key(hallmark):
    split = hallmark.split(" ")
    output = ""
    for element in split:
        output += "{}{}".format(element[0].upper(), element[1:])
    return output

# load up TfTargets
'''
with open("data/humanTFs_All.csv") as in_file:
    in_file.readline()  # consume headers

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
'''

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

# handle_miRNA_json("targetScan_miRNA_sets_entrez_hsa.json", "Target Scan")
# handle_miRNA_json("pita_miRNA_sets_entrez_hsa.json", "Pita")

def parse_mutation_name(input):
    regex = re.compile(r"(X|x|)([0-9]+)_(.+)")
    match = regex.match(input)

    if match:
        return (int(match.group(2)), match.group(3))
    else:
        return None

def parse_locus_name(input):
    regex = re.compile(r"(X|x|)(.*)_(.*)")
    match = regex.match(input)

    if match:
        return (match.group(2), match.group(3))
    else:
        return None

somatic_mutations = {}
locus_map = {}
locus_array = []
def interpret_causality_summary(filename, somatic_mutations):
    file = open(filename)
    file.readline() # consume header

    line = file.readline()
    while line:
        split = line.split(",")
        
        mutation = split[0]
        regulator = split[1]
        bicluster = split[2]

        parsed = parse_mutation_name(mutation)

        # create a somatic mutation based on the parsed name
        if parsed and parsed[0] in genes and parsed[0] not in somatic_mutations:
            somatic_mutations[parsed[0]] = SomaticMutation( \
                ext_id = genes[parsed[0]].id, \
                mutation_type = parsed[1], \
            )
            db.session.add(somatic_mutations[parsed[0]])
        elif parsed:
            print("unable to find entrez id for somatic mutation {} (tried searching for id '{}')".format(mutation, parsed[0]))
        else:
            parsed_locus = parse_locus_name(mutation)
            real_name = "{}_{}".format(parsed_locus[0], parsed_locus[1])

            if real_name not in locus_map:
                locus = Locus(\
                    locus_name = parsed_locus[0], \
                    mutation_type = parsed_locus[1], \
                )
                locus_map[real_name] = len(locus_array)
                locus_array.append(locus)
                db.session.add(locus)

                print("{} {}".format(parsed_locus[0], parsed_locus[1]))
                print("found locus {}, will not have an entrez id".format(mutation))
        
        line = file.readline()

    db.session.commit()

    # handle somatic mutations that have a locus
    for locus in locus_array:
        name = "{}_{}".format(locus.locus_name, locus.mutation_type)
        if name not in somatic_mutations:
            somatic_mutations[name] = SomaticMutation( \
                locus_id = locus.id, \
            )
            db.session.add(somatic_mutations[name])
    
    db.session.commit() # we have to commit in order to populate the somatic_mutations.id column


    # re-read the file to actually put in the causal flows
    file = open(filename)
    file.readline()

    line = file.readline()
    while line:
        split = line.split(",")

        mutation = split[0][1:] # string name
        regulator_entrez = int(split[1]) # entrez id
        bicluster_entrez = int(split[2]) # entrez id

        print(bicluster_entrez)
        print(mutation)

        leo_nb_atob = float(split[3])
        mlogp_m_atob = float(split[4])

        if bicluster_entrez in biclusters and mutation in somatic_mutations:
            causal_flow = CausalFlow(\
                somatic_mutation_id = somatic_mutations[mutation].id, \
                regulator_id = regulator_entrez, \
                regulator_type = "", \
                bicluster_id = biclusters[bicluster_entrez], \
                leo_nb_atob = leo_nb_atob, \
                mlogp_m_atob = mlogp_m_atob, \
            )
            db.session.add(causal_flow)
        
        line = file.readline()

    db.session.commit()

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

# parse phenotypes
with open('data/phenotypes_meso_noFilter.csv') as in_file:
    # read header
    header = in_file.readline().split(",")
    header_to_phenotype = {}
    for i in range(1, len(header)):
        name = header[i].strip()
        header_to_phenotype[name] = Phenotype( \
            name = name,
            long_name = name,
        )
        db.session.add(header_to_phenotype[name])
    
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
                        phenotype_id = header_to_phenotype[phenotype_name].id,
                        patient_id = current_patient.id,
                        phenotype_value = float(phenotype_value),
                    ))
                else: # handle strings
                    db.session.add(PhenoDatum( \
                        phenotype_id = header_to_phenotype[phenotype_name].id,
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

biclusters = {}
with open('data/postProcessed_clustersOfBiclusters_CNA_CNVkit.csv') as in_file:
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
        bicluster = Bicluster( \
            name = bicluster_name, \
            var_exp_fpc = split[header_dict["Var. Exp. First PC"]], \
            var_exp_fpc_p_value = split[header_dict["Var. Exp. First PC Perm. P-Value"]], \
            survival = split[header_dict["OS.covAgeSex"]], \
            survival_p_value = split[header_dict["OS.covAgeSex.p"]] \
        )
        db.session.add(bicluster)

        biclusters[bicluster_name] = bicluster
        
        # create a BicGene
        bic_gene_values.append((bicluster_name, int(split[header_dict["Genes"]])))

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
        elif tuple_val[1] not in genes:
            print("bic_gene warning: could not find gene id {}".format(tuple_val[1]))
        else:
            db.session.add(BicGene( \
                bicluster_id = biclusters[tuple_val[0]].id, \
                gene_id = genes[tuple_val[1]].id \
            ))
    
    # create all the BicGo's
    for tuple_val in bic_go_values:
        if tuple_val[1] != "NA":
            split = tuple_val[1].split(";")
            for go_val in split:
                db.session.add(BicGo( \
                    bicluster_id = biclusters[tuple_val[0]].id, \
                    go_bp_id = geneOntology[go_val].id \
                ))
    
    # create all the BicPat's
    for tuple_val in bic_pat_values:
        split = tuple_val[1].split(" ")
        for patient in split:
            db.session.add(BicPat( \
                bicluster_id = biclusters[tuple_val[0]].id, \
                patient_id = patients[patient].id \
            ))

    # create all the BicHal's
    for tuple_val in bic_halmark_values:
        bicluster_name = tuple_val[0]
        hallmark_name = tuple_val[1]
        hallmark_value = tuple_val[2]

        if hallmark_value >= 0.8:
            db.session.add(BicHal( \
                bicluster_id = biclusters[bicluster_name].id, \
                hallmark_id = hallmarks[hallmark_name].id \
            ))
    
    # create all the Replication's
    for tuple_val in replication_values:
        bicluster_name = tuple_val[0]
        var_exp_fpc = tuple_val[1]
        var_exp_fpc_p_value = tuple_val[2]
        survival = tuple_val[3]
        survival_p_value = tuple_val[4]

        db.session.add(Replication( \
            bicluster_id = biclusters[bicluster_name].id, \
            study = "MESO TCGA", \
            var_exp_fpc = var_exp_fpc, \
            var_exp_fpc_p_value = var_exp_fpc_p_value, \
            survival = survival, \
            survival_p_value = survival_p_value \
        ))
    
    # create all TF regulators
    for tf_regulator_array in tf_regulators:
        bicluster_name = tf_regulator_array[0]
        for index in range(1, len(tf_regulator_array)):
            tf_regulator = tf_regulator_array[index]
            db.session.add(TfRegulator( \
                bicluster_id = biclusters[bicluster_name].id, \
                tf_id = genes[tf_regulator[0]].id, \
                r_value = tf_regulator[1], \
                p_value = tf_regulator[2], \
            ))
    
    # create all miRNA regulators
    for miRNA_regulator_array in miRNA_regulators:
        bicluster_name = miRNA_regulator_array[0]
        for index in range(1, len(miRNA_regulator_array)):
            mimat = miRNA_regulator_array[index]
            if mimat in miRNAs:
                miRNA_entry = miRNAs[mimat]
                db.session.add(MirnaRegulator(
                    bicluster_id = biclusters[bicluster_name].id, \
                    mirna_id = miRNA_entry.id, \
                ))
            else:
                print("missing MIMAT {}, can't link to bicluster table".format(mimat))
    
    db.session.commit()


# handle somatic mutations and causal flows
# interpret_causality_summary("./data/causality_CNA_final_8_13_2019/causalitySummary_pita.csv", somatic_mutations)
# interpret_causality_summary("./data/causality_CNA_final_8_13_2019/causalitySummary_targetscan.csv", somatic_mutations)
# interpret_causality_summary("./data/causality_CNA_final_8_13_2019/causalitySummary_tfbs_db.csv", somatic_mutations)