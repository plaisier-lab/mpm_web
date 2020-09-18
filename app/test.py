import gzip
import pandas
import numpy as np
import mysql.connector
from flask import Flask, Response, url_for, redirect, render_template, request, session, flash, jsonify
import flask
import itertools
import math
from scipy.stats import hypergeom

NUM_PARTS = 5

convert = {'Evading apoptosis': 'cellDeath.gif', 'Evading immune detection': 'avoidImmuneDestruction.gif', 'Genome instability and mutation': 'genomicInstability.gif', 'Insensitivity to antigrowth signals': 'evadeGrowthSuppressors.gif', 'Limitless replicative potential': 'immortality.gif',
    'Reprogramming energy metabolism': 'cellularEnergetics.gif', 'Self sufficiency in growth signals': 'sustainedProliferativeSignalling.gif', 'Sustained angiogenesis': 'angiogenesis.gif', 'Tissue invasion and metastasis': 'invasion.gif', 'Tumor promoting inflammation': 'promotingInflammation.gif'}

def phyper(x, M, n, bigN, lower_tail=False):
    """uses scipy.stats to compute hypergeometric overlap p-value"""
    return 1 - hypergeom.cdf(float(x), float(M), float(n), float(bigN))

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
        'host': "localhost",
        'port': "32000",
        'database': "mpm_web",
    }
    return mysql.connector.connect(**config)


db = dbconn()
c = db.cursor()

def submat_data(submat, col_indexes):
    """given a sub matrix and a list of column indexes
    that specify the columns, of the matrix, return a list
    of (col_idx, median, min, max, lower_quartile, upper_quartile)
    tuples
    """
    col_medians = np.median(submat, axis=0)
    col_mins = np.min(submat, axis=0)
    col_maxs = np.max(submat, axis=0)
    col_upper_quarts = np.percentile(submat, q=75.0, axis=0)
    col_lower_quarts = np.percentile(submat, q=25.0, axis=0)

    print(col_medians)

    data = [[idx,
             col_mins[i],
             col_lower_quarts[i],
             col_medians[i],
             col_upper_quarts[i],
             col_maxs[i]
             ]
            for i, idx in enumerate(col_indexes)]
    return sorted(data, key=lambda x: x[3])

def subtype_enrichment(cursor, cluster_id):
    cursor.execute("""SELECT g.entrez from bic_gene bg
        JOIN gene g ON bg.gene_id=g.id WHERE bicluster_id=%s""", [cluster_id])
    genes = [row[0] for row in cursor.fetchall()]
    gene_count = len(genes)
    genes_string = "(%s)" % ','.join(map(lambda p: '%s' % p, genes))

    # get patients based on the given bicluster id
    cursor.execute("""SELECT name FROM bic_pat bp
        JOIN patient p ON bp.patient_id=p.id WHERE bicluster_id=%s""",
        [cluster_id])
    included_patients = [row[0] for row in cursor.fetchall()]
    included_patients_string = "(%s)" % ','.join(map(lambda p: '\'%s\'' % p, included_patients))

    # get excluded patients based on the given bicluster id
    cursor.execute("""SELECT name FROM patient WHERE id NOT IN
        (SELECT patient_id FROM bic_pat WHERE bicluster_id=%s)""",
        [cluster_id])
    excluded_patients = [row[0] for row in cursor.fetchall()]
    excluded_patients_string = "(%s)" % ','.join(map(lambda p: '\'%s\'' % p, excluded_patients))

    # above is exactly like boxplot
    # now make a phenotype map
    # cursor.execute("""SELECT p.name, pt.name FROM patient p JOIN phenotypes pt ON p.phenotype_id=pt.id
    #    WHERE pt.name <> 'NA'""")
    cursor.execute("""SELECT p.name, pd.phenotype_string FROM pheno_data pd
        JOIN patient p ON pd.patient_id=p.id
        JOIN phenotype pt ON pd.phenotype_id=pt.id
        WHERE pt.name="histology_WHO" 
            AND pd.phenotype_string <> 'NA'""")
    # histology_WHO phenotype lookup ^^^

    ptmap = {patient: phenotype for patient, phenotype in cursor.fetchall()} # don't change
    all_patients = {patient for patient in ptmap.keys()} # don't change
    phenotypes = {phenotype for phenotype in ptmap.values()} # don't change

    # we use the submat_data function to sort our patients
    # lol part 1. also, i hate this sql statement. despite the obvious injection risk, this is fine. genes_string and included_patients_string is generated from an unexploitable source
    cursor.execute("SELECT value FROM gene_expression ge " # at the end of this, we're getting values from the gene_expression table
        + "JOIN gene g ON ge.gene_id=g.id " # we want to search based on gene entrez, so we JOIN the gene table with the gene expression table
        + "JOIN patient p ON ge.patient_id=p.id " # we want to also search based on patient name, so we JOIN the patient table with the gene expression table
        + "WHERE g.entrez IN " + genes_string + " " # we're looking for gene entrez ids that we found in the bic_gene table based on the input bicluster id
        + "AND p.name IN " + included_patients_string + ";") # we're looking for patient names that we found in the patient table based on the given bicluster idnames that we found in the patient table based on the given bicluster id
    included_patients_result = [row[0] for row in cursor.fetchall()]
    in_submat = np.ndarray((gene_count, len(included_patients)), dtype=float, buffer=np.array(included_patients_result), offset=0, strides=None, order=None) # we want a submatrix (whatever that is) of all the different gene expression values we found
    in_data = submat_data(in_submat, included_patients) # throw the submatrix at this function, does statistics on the data
    sorted_in_names = [row[0] for row in in_data] # NOTE: due to my changes to account for the gene_expression table, row[0] is now a PATIENT NAME, not a PATIENT INDEX.

    # lol part 2. also, i hate this sql statement. despite the obvious injection risk, this is fine. genes_string and excluded_patients_string is generated from an unexploitable source
    cursor.execute("SELECT value FROM gene_expression ge " # at the end of this, we're getting values from the gene_expression table
        + "JOIN gene g ON ge.gene_id=g.id " # we want to search based on gene entrez, so we JOIN the gene table with the gene expression table
        + "JOIN patient p ON ge.patient_id=p.id " # we want to also search based on patient name, so we JOIN the patient table with the gene expression table
        + "WHERE g.entrez IN " + genes_string + " " # we're looking for gene entrez ids that we found in the bic_gene table based on the input bicluster id
        + "AND p.name IN " + excluded_patients_string + ";") # we're looking for patient names that we DIDN'T find in the patient table based on the given bicluster idfor patient names that we DIDN'T find in the patient table based on the given bicluster id
    excluded_patients_result = [row[0] for row in cursor.fetchall()]
    ex_submat = np.ndarray((gene_count, len(excluded_patients)), dtype=float, buffer=np.array(excluded_patients_result), offset=0, strides=None, order=None) # we want a submatrix (whatever that is) of all the different gene expression values we found 
    ex_data = submat_data(ex_submat, excluded_patients)    
    sorted_ex_names = [row[0] for row in ex_data] # NOTE: due to my changes to account for the gene_expression table, row[0] is now a PATIENT NAME, not a PATIENT INDEX.

    # sorted by median pValue
    sorted_patient_names = sorted_in_names + sorted_ex_names

    # group patients into phenotype groups.
    # NOTE: the ptmap items need to be sorted, otherwise groupby fails to group correctly
    pt_patients = itertools.groupby(sorted(ptmap.items(), key=lambda pair: pair[1]), key=lambda pair: pair[1])
    pt_patients = {phenotype: set(map(lambda p: p[0], patients)) for phenotype, patients in pt_patients}

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

            print start, end
            print "part %d, phenotype: %s, x=%d, M=%d, n=%d, N=%d, pvalue=%f, phyper=%f" % (i, phenotype, x, M, n, bigN, pvalue, phyper(x, M, n, bigN))

            part_pvalues[phenotype] = pvalue
        pvalues.append(part_pvalues)

    return pvalues, min_pvalue, max_pvalue

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
                if tf > 0:
                    print(tf)
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

pvalues, min_pvalue, max_pvalue = subtype_enrichment(c, 21)
# print(pvalues)

# print(__get_muts(c, 12442, "BAP1"))



'''SELECT value FROM gene_expression ge 
    JOIN gene g ON ge.gene_id=g.id
    JOIN patient p ON ge.patient_id=p.id 
    WHERE g.entrez IN (3569,2865,407013,4327,3491,7852,9308,7832,8013,3164,9592,6004,1831,9021,4929,1839,1960,387763,51655,7538,4616,5008,124599,467,23764,9507,64651,6446,5997,5996,6615,3725,9314,1959,1958,1880,407018,8870,345274,5187,25976,407010,4739,28984,84807,1844,2354,153020,2353,728) AND p.name IN ("602PT","603PT","607PT","610PT","612PT","613PT","617PT","624PT","628PT","634PT","636PT","651PT","657PT","658PT","660PT","661PT","667PT","M101PT","M12PT","M13PT","M15PT","M16PT","M18PT","M19PT","M23PT","M24PT","M25PT","M27PT","M29PT","M30PT","M31PT","M33PT","M36PT","M37PT","M38PT","M39PT","M41PT","M42PT","M43PT","M44PT","M47PT","M4PT","M50PT","M51PT","M53PT","M54PT","M57PT","M58PT","M59PT","M604PT","M605PT","M608PT","M616PT","M619PT","M621PT","M622PT","M627PT","M62PT","M630PT","M631PT","M632PT","M637PT","M63PT","M641PT","M646PT","M64PT","M653PT","M656PT","M659PT","M65PT","M663PT","M665PT","M66PT","M670PT","M671PT","M674PT","M675PT","M678PT","M67PT","M681PT","M682PT","M683PT","M687PT","M690PT","M693PT","M695PT","M696PT","M701PT","M702PT","M703PT","M706PT","M70PT","M710PT","M711PT","M717PT","M718PT","M71PT","M76PT","M7PT","M83PT","M8PT","M92PT","M94PT","M95PT","M97PT","M99PT","M9PT")


(3569,2865,407013,4327,3491,7852,9308,7832,8013,3164,9592,6004,1831,9021,4929,1839,1960,387763,51655,7538,4616,5008,124599,467,23764,9507,64651,6446,5997,5996,6615,3725,9314,1959,1958,1880,407018,8870,345274,5187,25976,407010,4739,28984,84807,1844,2354,153020,2353,728)'''