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

import re, cPickle, os
from copy import deepcopy
#from cMonkeyWrapper import cMonkeyWrapper
#from bicluster import bicluster
#from pssm import pssm
import pandas as pd

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
    if len(retMe)==0:
        miRNA_shorter = '-'.join(miRNA.split('-')[0:3])
        for i in dict1.keys():
            if compareMiRNANames(miRNA_shorter, i):
                retMe.append(miRNAIDs[i])
    return retMe

def compareMiRNANames(a, b):
    if a==b:
        return 1
    if len(a)<len(b):
        if a[-3:]=='-3p':
            re1 = re.compile(a+'[a-oq-z]?(-\d)?-3p$')
        else:
            re1 = re.compile(a+'[a-oq-z]?(-\d)?-5p$')
        if re1.match(b):
            return 1
    else:
        if b[-3:]=='-3p':
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
with open('data/gene2entrezId.csv','r') as inFile:
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
with open('data/hsa.mature.fa','r') as inFile:
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
            print('Uh oh!',splitUp)

####################################################
## Python script to populate the mpm_web database ##
####################################################

### Load up Genentech Gene Expression dataset

## Read in gene expression matrix (Genentech)
gexp1 = pd.read_csv('data/mesothelioma_norm.txt', header=0, index_col=0, sep=' ')

## Load up patient information (Genentech)
patients = {}
for patient1 in gexp1.columns:
    p1 = Patient(name=patient1)
    db.session.add(p1)
    patients[patient1] = deepcopy(p1)

## Load up gene information (Genentech)
genes = {}
for gene1 in gexp1.index:
    if str(gene1) in entrez2symbol:
        g1 = Gene(symbol=entrez2symbol[str(gene1)],entrez=int(gene1))
        db.session.add(g1)
        genes[gene1] = deepcopy(g1)

## Add new expresssion dataset
ed1 = ExpDataset(name='Genentech_MPM', type='gexp')
db.session.add(ed1)
gt_gexp = deepcopy(ed1)

## Load up gene expression data
for patient1 in gexp1.columns:
    for gene1 in gexp1.index:
        if str(gene1) in genes:
            ge1 = GeneExpression(exp_dataset_id=gt_gexp.id, patient_id=patients[patient1].id, gene_id=genes[gene1].id, value=gexp1[patient1][gene1])
            db.session.add(ge1)

db.session.commit()

### Add hallmarks of Cancer
hallmark_names = ['Sustained angiogenesis', 'Insensitivity to antigrowth signals', 'Evading apoptosis', 'Limitless replicative potential', 'Evading immune detection', 'Tissue invasion and metastasis', 'Self sufficiency in growth signals', 'Tumor promoting inflammation', 'Reprogramming energy metabolism', 'Genome instability and mutation']
hallmarks = {}
for hn1 in hallmark_names:
    hm1 = Hallmark(name=hn1)
    db.session.add(hm1)
    hallmarks[hn1] = deepcopy(hm1)

### Add gene ontology and GO to gene mappings
## GO_BP
geneOntology = {}
with open('data/go.obo','r') as inFile:
    while 1:
        line = inFile.readline()
        if not line:
            break
        
        if line.strip()=='[Term]':
            id1 = inFile.readline().strip().split('id: ')[1]
            name = inFile.readline().strip().split('name: ')[1]
            category = inFile.readline().strip().split('namespace: ')[1]
            if category=='biological_process':
                go1 = GoBp(go_id=id1, name=name)
                geneOntology[id1] = deepcopy(go1)
        
        if line.strip()[0:6]=='alt_id' and category=='biological_process':
            id1 = line.strip().split('alt_id: ')[1]
            go1 = GoBp(go_id=id1, name=name)
            geneOntology[id1] = deepcopy(go1)

## Gene -> GO_BP
with open('data/gene2go.hsa','r') as inFile:
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split('\t')
        if splitUp[1] in genes and splitUp[2] in geneOntology:
            gg1 = GoGene(go_bp_id=geneOntology[splitUp[2]].id, gene_id=genes[splitUp[1]].id)
            db.session.add(gg1)

db.session.commit()

### miRNAs
miRNAs = {}
for miR1 in miRNAIDs:
    m1 = Mirna(mature_seq_id=miRNAIDs[miR1], name=miR1)
    db.session.add(m1)
    miRNAs[miRNAIDs[miR1]] = deepcopy(m1)

### Prior information

## Load DisGeNET
mesos = ['Mesothelioma','Malignant mesothelioma','Malignant Pleural Mesothelioma','Mesothelioma (malignant, clinical disorder) (disorder)','Malignant Mesothelioma of Peritoneum','Pleural Mesothelioma','Mesothelioma, biphasic, malignant','Sarcomatoid Mesothelioma']
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
            gp1 = GenePrior(gene_id=genes[int(splitUp[0])].id, source='DisGeNET', description=splitUp[5], pmid=pmid)
            db.session.add(gp1)

## Load HMDD (miRNA_Name,PMID,Description)
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
           mp1 = MirnaPrior(mirna_id=miRNAs[m1].id, source='hmdd', description=splitUp[2], pmid=int(splitUp[1]))
           db.session.add(mp1)


## Load dbDEMC
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
           mp1 = MirnaPrior(mirna_id=miRNAs[m1].id, source='dbDEMC', description=splitUp[1], pmid=int(splitUp[2]))
           db.session.add(mp1)

db.session.commit()

