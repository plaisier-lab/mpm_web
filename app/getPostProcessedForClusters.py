#import statsmodels.formula.api as sm
#import numpy as np
from sklearn.linear_model import LinearRegression
from scipy import stats
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns
import json

# Read in clusters of biclusters
"""
clusters = {}
inFile = open('clustersOfBiclusters_VALL.csv','r')
inFile.readline() # Get rid of header
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    if not splitUp[1] in clusters:
        clusters[splitUp[1]] = []
    clusters[splitUp[1]].append(splitUp[0].strip('"').lower())
inFile.close()
"""
def getSignificantHallmarks(b1):
    # SustainedAngiogenesis,InsensitivityToAntigrowthSignals,EvadingApoptosis,LimitlessReplicativePotential,EvadingImmuneDetection,TissueInvasionAndMetastasis,SelfSufficiencyInGrowthSignals,TumorPromotingInflammation,ReprogrammingEnergyMetabolism,GenomeInstabilityAndMutation
    hallmarks = []
    for hm1 in ['SustainedAngiogenesis','InsensitivityToAntigrowthSignals','EvadingApoptosis','LimitlessReplicativePotential','EvadingImmuneDetection','TissueInvasionAndMetastasis','SelfSufficiencyInGrowthSignals','TumorPromotingInflammation','ReprogrammingEnergyMetabolism','GenomeInstabilityAndMutation']:
        if not b1[hm1]=='NA' and float(b1[hm1])>=0.8:
            hallmarks.append(hm1)
    return hallmarks

def getMutName(name1):
    tmp = name1.split('_')
    #if name1.find('p')>0 or name1.find('q')>0:
    #    #tmp = [tmp[0]+'_'+tmp[1], tmp[2]]
    #    tmp = [tmp[0], tmp[2]]
    return tmp

def mapMut2Loci(mut1, loci):
    if not (mut1.find('p')>0 or mut1.find('q')>0 or mut1.find('PAM')>0):
        tmp = mut1.split('_')
        tmp2 = [i for i in loci if tmp[0] in loci[i]]
        if len(tmp2)>0:
            return tmp2[0].split('_')[0]
        else:
            return False
    else:
        return False

# Calculate percentage of target genes by using the maximum for the TF families expansion
def percentTargets(tf, bicGenes, tfbsDb, id2family, family2id, id2motif):
    motifs = []
    if tf in id2motif:
        motifs += id2motif[tf]
    if tf in id2family:
        for tf2 in family2id[id2family[tf]]:
            if tf2 in id2motif:
                motifs += id2motif[tf2]
    percTarget = 0
    for motif in motifs:
        if motif in tfbsDb:
            tmp = float(len(list(set(bicGenes).intersection(tfbsDb[motif]))))/float(len(bicGenes))
            if tmp > percTarget:
                percTarget = tmp
    return percTarget

# Reading TF regulator to TF target gene relationships into Python
# The json library we import takes care of most of the work
with open('tfbsDb_plus_and_minus_5000_entrez.json', 'r') as f:
    tfbsDb = json.load(f)

# Load up id to motif thesaurus
id2motif = {}
with open('humanTFs_All.CSV','r') as inFile:
    header = inFile.readline().strip().split(',')
    while 1:
        inLine = inFile.readline()
        if not inLine:
            break
        split = inLine.strip().split(',')
        if not split[2] in id2motif:
            id2motif[split[2]] = []
        id2motif[split[2]].append(split[0])

# TF family expansion
family2id = {}
id2family = {}
with open('tfFamilies.CSV','r') as inFile:
    header = inFile.readline()
    while 1:
        inLine = inFile.readline()
        if not inLine:
            break
        split = inLine.split(',')
        split[2] = split[2].replace(' ',',').strip().split(',')
        family2id[split[0]] = split[2]
        for splitId in split[2]:
            id2family[splitId] = split[0]

# Read in entrez ID to gene symbol translator
gene2entrezId = {}
entrezId2gene = {}
inFile = open('gene2entrezId.csv','r')
while 1:
    line = inFile.readline()
    if not line:
        break
    splitUp = line.strip().split(',')
    gene2entrezId[splitUp[0]] = splitUp[1]
    entrezId2gene[splitUp[1]] = splitUp[0]
inFile.close()

# Determine which mutations are solid AMPs, DELs, Acts, and LoFs
loci_Bueno = {}
with open('../oncoMerge_final/CNA_loci_0.05.csv','r') as inFile:
    inFile.readline() # get rid of header
    while 1:
        inLine = inFile.readline()
        if not inLine:
            break
        tmp = inLine.strip().split(',')
        loci_Bueno[tmp[0]] = tmp[1].split(' ')

loci_TCGA = {}
with open('../../mesothelioma_TCGA/oncoMerged_MESO/CNA_loci_0.05.csv','r') as inFile:
    inFile.readline() # get rid of header
    while 1:
        inLine = inFile.readline()
        if not inLine:
            break
        tmp = inLine.strip().split(',')
        loci_TCGA[tmp[0]] = tmp[1].split(' ')

# Mapping loci from TCGA to Bueno
loci_mapping = {}
loci_mapping_rev = {}
for loci1 in loci_Bueno:
    for loci2 in loci_TCGA:
        if len(set(loci_Bueno[loci1]).intersection(loci_TCGA[loci2]))>0:
            tmp1 = loci1 #.split('_')[0]
            tmp2 = loci2 #.split('_')[0]
            if not tmp1 in loci_mapping:
                loci_mapping[tmp1] = []
            if not tmp2 in loci_mapping[tmp1]:
                loci_mapping[tmp1].append(tmp2)
            if not tmp2 in loci_mapping_rev:
                loci_mapping_rev[tmp2] = []
            if not tmp1 in loci_mapping_rev[tmp2]:
                loci_mapping_rev[tmp2].append(tmp1)

genMuts = []
genesInMuts = []
with open('../oncoMerge_final/oncoMerged_0.05.csv','r') as inFile:
    inFile.readline() # get rid of header
    while 1:
        inLine = inFile.readline()
        if not inLine:
            break
        tmp = inLine.strip().split(',')[0]
        genMuts.append(tmp)
        tmp1 = tmp.split('_')
        tmp2 = tmp1[0]+'_'+tmp1[1]
        if tmp.find('p')>0 or tmp.find('q')>0 and not tmp2 in entrezId2gene:
            #print tmp,tmp2
            #entrezId2gene[tmp2] = tmp2
            tmp3 = tmp2.split('_')[0]
            entrezId2gene[tmp3] = tmp3
        else:
            genesInMuts.append(tmp1[0])

genMuts_TCGA = []
genesInMuts_TCGA = []
with open('../../mesothelioma_TCGA/oncoMerged_MESO_8_13_2019/MESO_oncoMerged_0.05.csv','r') as inFile:
    inFile.readline() # get rid of header
    while 1:
        inLine = inFile.readline()
        if not inLine:
            break
        tmp = inLine.strip().split(',')[0]
        genMuts_TCGA.append(tmp)
        if tmp in loci_mapping_rev and not loci_mapping_rev[tmp]:
            genMuts.append(loci_mapping_rev[tmp])
        tmp1 = tmp.split('_')
        tmp2 = tmp1[0]+'_'+tmp1[1]
        if tmp.find('p')>0 or tmp.find('q')>0 and not tmp2 in entrezId2gene:
            #entrezId2gene[tmp2] = tmp2
            tmp3 = tmp2.split('_')[0]
            entrezId2gene[tmp3] = tmp3
        else:
            genesInMuts_TCGA.append(tmp1[0])

# Choose which to keep: either a PAM or also in tcgaMuts
keepMut = []
#for mut1 in genMuts:
#    if (mut1.find('p')>0 or mut1.find('q')>0):
#        #print(mut1, len([i for i in loci_Bueno[mut1] if (i+'_Act' in genMuts or i+'_LoF' in genMuts)]),[i for i in loci_Bueno[mut1] if (i+'_Act' in genMuts or i+'_LoF' in genMuts)])
#        if mut1 in loci_mapping and loci_mapping[mut1][0] in genMuts_TCGA:
#            #print(loci_mapping[mut1][0])
#            #print(mut1,len([i for i in loci_Bueno[mut1] if (i+'_Act' in genMuts or i+'_LoF' in genMuts)]),'in')
#            if len([i for i in loci_Bueno[mut1] if (i+'_Act' in genMuts or i+'_LoF' in genMuts)])==0:
#                keepMut.append(mut1)
#    # Include PAM mutations
#    elif not (mut1.find('LoF')>0 or mut1.find('Act')>0):
#        keepMut.append(mut1)

keepMut = sorted(list(set(genMuts+genMuts_TCGA)))

# Read in miRBase to miRNA name
id2miR = {}
with open('hsa.mature.fa','r') as inFile:
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(' ')
        id2miR[splitUp[1]] = splitUp[0]

# Read in post processed files
print 'Load up cMonkey SYGNAL output...'
biclusters = {}
biclusterHeader = ''
for i in ['pita','targetscan','tfbs_db']:
    inFile = open('postprocessed_CNA/postProcessed_meso_'+i+'.csv','r')
    biclusterHeader = inFile.readline() # Get rid of header
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.split(',')
        biclusters[i+'_'+splitUp[0]] = line.strip()
    inFile.close()

    inFile = open('genes_conds/cluster.members.genes_'+i+'.txt','r')
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(' ')
        j = splitUp.pop(0)
        #biclusters[i+'_'+j] = biclusters[i+'_'+j] + ',' + ';'.join(splitUp)
    inFile.close()

    inFile = open('genes_conds/cluster.members.conditions_'+i+'.txt','r')
    while 1:
        line = inFile.readline()
        if not line:
            break
        splitUp = line.strip().split(' ')
        j = splitUp.pop(0)
        biclusters[i+'_'+j] = biclusters[i+'_'+j] + ',' + ' '.join(splitUp)
    inFile.close()
print 'Done.'

# Load up Genentech causality TFs
leo_nb_AtoB_cutoff = 0.3
print 'Load up Genentech causal relationships...'
causalTfs = {}
#causalMiRNAs = {}
biMutAndTf = {}
for run1 in ['pita','targetscan','tfbs_db']:
    with open('causality_CNA_final_8_13_2019/causalitySummary_'+run1+'.csv','r') as inFile:
        header = inFile.readline().strip().split(',')
        biMutAndMiRNA = {}
        while 1:
            line = inFile.readline().strip().split(',')
            if not line or len(line)==1:
                break
            line = dict(zip(header,line))
            if line['Mutation'].lstrip('X') in keepMut and float(line['leo.nb.AtoB'])>=leo_nb_AtoB_cutoff:
                line['Bicluster'] = run1+'_'+line['Bicluster']
                if not line['Bicluster'] in causalTfs:
                    causalTfs[line['Bicluster']] = []
                #if not line['Bicluster'] in causalMiRNAs:
                #    causalMiRNAs[line['Bicluster']] = []
                causalTfs[line['Bicluster']].append(line['Regulator'])
                if not line['Bicluster'] in biMutAndTf:
                    biMutAndTf[line['Bicluster']] = []
                biMutAndTf[line['Bicluster']].append({'mut':line['Mutation'].lstrip('X'),'tf':line['Regulator']})
print 'Done.\n'

# Load up TCGA causality TFs and miRNA
causalTfs_rep = {}
causalMiRNAs_rep = {}
biMutAndTf_rep = {}
biMutAndMiRNA_rep = {}
#with open('../../mesothelioma_TCGA/causal_v8/summaryCausality_CNV_3_1_2019_0.3_0.05_cleanedUp.csv','r') as inFile:
with open('../../mesothelioma_TCGA/causal_v9/summaryCausality_CNV_8_16_2019_0.3_0.05_cleanedUp.csv','r') as inFile:
    header = inFile.readline().strip().split(',')
    while 1:
        line = inFile.readline().strip().split(',')
        if not line or len(line)==1:
            break
        line = dict(zip(header,line))
        # Convert loci into hg38 and get rid of those that don't map
  TODO       if (line['Mutation'].find('p')>0 or line['Mutation'].find('q')>0) and line['Mutation'] in loci_mapping_rev and len(loci_mapping_rev[line['Mutation']])==1:
            line['Mutation'] = loci_mapping_rev[line['Mutation']]
        if line['Mutation'] in keepMut and float(line['leo.nb.AtoB'])>=leo_nb_AtoB_cutoff:
            if not line['Regulator'][0:5]=='MIMAT':
                if not line['Bicluster'] in causalTfs_rep:
                    causalTfs_rep[line['Bicluster']] = []
                if not line['Regulator'] in causalTfs_rep[line['Bicluster']]:
                    causalTfs_rep[line['Bicluster']].append(line['Regulator'])
                if not line['Bicluster'] in biMutAndTf_rep:
                    biMutAndTf_rep[line['Bicluster']] = []
                biMutAndTf_rep[line['Bicluster']].append({'mut':line['Mutation'].lstrip('X'),'tf':line['Regulator']})
            else:
                if not line['Bicluster'] in causalMiRNAs_rep:
                    causalMiRNAs_rep[line['Bicluster']] = []
                if not line['Regulator'] in causalMiRNAs_rep[line['Bicluster']]:
                    causalMiRNAs_rep[line['Bicluster']].append(line['Regulator'])
                if not line['Bicluster'] in biMutAndMiRNA_rep:
                    biMutAndMiRNA_rep[line['Bicluster']] = []
                biMutAndMiRNA_rep[line['Bicluster']].append({'mut':line['Mutation'],'miRNA':line['Regulator']})
print 'Done.\n'

# Read in correlaion information
cor_miR = {}
with open('miRNA_cor_mesoTCGA/eig_miR_cor.csv','r') as c1File:
    header = [i.strip('"') for i in c1File.readline().strip().split(',')]
    header.pop(0)
    while 1:
        inLine = c1File.readline()
        if not inLine:
            break
        splitUp = inLine.strip().split(',')
        bic1 = splitUp.pop(0).strip('"')
        cor_miR[bic1] = dict(zip(header,[float(i) for i in splitUp]))
pv_miR = {}
with open('miRNA_cor_mesoTCGA/eig_miR_pv.csv','r') as pvFile:
    header = [i.strip('"') for i in pvFile.readline().strip().split(',')]
    header.pop(0)
    while 1:
        inLine = pvFile.readline()
        if not inLine:
            break
        splitUp = inLine.strip().split(',')
        bic1 = splitUp.pop(0).strip('"')
        pv_miR[bic1] = dict(zip(header,[float(i) for i in splitUp]))


# Load up JACKS analysis data

''' area of focus: '''
# TODO ./data/CRISPER_Cas9_Screen!!!
essScore = pd.read_csv('MPM_gene_JACKS_results.txt', sep='\t', header=0, index_col=0)
pvals = pd.read_csv('MPM_gene_pval_JACKS_results.txt', sep='\t', header=0, index_col=0)
lfc = pd.read_csv('JACKS/MPM_logfoldchange_means.txt', sep='\t', header=0, index_col=0)
lfc[['M12T','M24']] = 2**lfc[['M12T','M24']]
lfc_mean = np.log2(lfc.groupby('gene').mean().loc[list(pvals.index)])

# Calculate 2-tailed p-values
pvals_2tailed = pvals
for gene in pvals.index:
    for sample in ['M12T','M24']:
        if essScore[sample].loc[gene]>0:
            pvals_2tailed[sample].loc[gene] = 2*(1-pvals_2tailed[sample].loc[gene])
        else:
            pvals_2tailed[sample].loc[gene] = 2*(pvals_2tailed[sample].loc[gene])

# Calculate FDR
M12T_FDR = statsmodels.stats.multitest.multipletests(pvals_2tailed['M12T'], alpha=0.05, method='fdr_bh')
M24_FDR = statsmodels.stats.multitest.multipletests(pvals_2tailed['M24'], alpha=0.05, method='fdr_bh')

# Hits M12T
upM12T = essScore.index[(essScore['M12T']>=1) & (M12T_FDR[2]<=0.05)]
downM12T = essScore.index[(essScore['M12T']<=-1) & (M12T_FDR[2]<=0.05)]
# len(upM12T) = 35
# len(downM12T) = 1660

# Hits M24
upM24 = essScore.index[(essScore['M24']>=1) & (M24_FDR[2]<=0.05)]
downM24 = essScore.index[(essScore['M24']<=-1) & (M24_FDR[2]<=0.05)]
# len(upM24) = 45
# len(downM24) = 1046

# Load common essential genes
# TODO ./data/CRISPER_Cas9_Screen/DepMap_20Q1/!!!
#comEss = list(pd.read_csv('mesoCRISPR/DepMap_20Q1/common_essentials.csv', header=0, sep=' ').index)
# TODO ./data/CRISPER_Cas9_Screen/DepMap_20Q1/!!! we ended up using this version of the file in the end
comEss = list(pd.read_csv('../DepMap_20Q1/Achilles_common_essentials.csv', header=0, sep=' ').index)
# len(comEss) = 2148

notComM12T = sorted(set(list(downM12T)).difference(comEss))
# len(notComM12T) = 398

notComM24 = sorted(set(list(downM24)).difference(comEss))
# len(notComM24) = 261
''' end of area of focus '''



notComBoth = sorted(set(notComM12T).intersection(notComM24))
# len(notComBoth) = 155
# Dump combined postProcessed file with clusters IDs
tfs = {}
correlatedTfs = {}
bestTfs = {}
correspondentTfs = {}
correspondentMiRNAs = {}
miRNAs = {}
cor_miRNAs = {}
goBP = {}
hallmarksBi = {}
siRNATFs = {}
miRNAmimics = {}
dumpBics = []
writeMe = ['bicluster,'+biclusterHeader.strip()+',Conditions']
for j in biclusters:
    # Process out tfs and miRNAs
    splitUp = dict(zip([i.strip('"') for i in biclusterHeader.strip().split(',')],biclusters[j].split(',')))
    if float(splitUp['Var. Exp. First PC'])>=0.3 and float(splitUp['Var. Exp. First PC Perm. P-Value'])<=0.05 and float(splitUp['mesoTCGA_pc1.perm.p'])<=0.05 and ((float(splitUp['OS.covAgeSex.p'])<=0.05 and float(splitUp['mesoTCGA_OS.age.sex.p'])<=0.05 and ((float(splitUp['OS.covAgeSex'])/abs(float(splitUp['OS.covAgeSex'])))==(float(splitUp['mesoTCGA_OS.age.sex'])/abs(float(splitUp['mesoTCGA_OS.age.sex']))))) or len([i for i in ['SustainedAngiogenesis', 'InsensitivityToAntigrowthSignals', 'EvadingApoptosis', 'LimitlessReplicativePotential', 'EvadingImmuneDetection', 'TissueInvasionAndMetastasis', 'SelfSufficiencyInGrowthSignals', 'TumorPromotingInflammation', 'ReprogrammingEnergyMetabolism', 'GenomeInstabilityAndMutation'] if not splitUp[i]=='NA' and float(splitUp[i])>=0.8])>0):
        # setting up datastructures
        tfs[j] = {'MEME':[], 'WEEDER':[], 'TFBS_DB':[]}
        correlatedTfs[j] = {'MEME':[], 'WEEDER':[], 'TFBS_DB':[]}
        bestTfs[j] = {'MEME':[], 'WEEDER':[], 'TFBS_DB':[]}
        miRNAs[j] = {'pita':[],'targetscan':[],'miRvestigator':[]}
        cor_miRNAs[j] = {'pita':[],'targetscan':[],'miRvestigator':[]}
        goBP[j] = []
        siRNATFs[j] = []
        miRNAmimics[j] = []
        hallmarksBi[j] = []
        writeMe.append(j+','+biclusters[j])
        # MEME1
        if (not splitUp['MEME Motif1 E-Value']=='NA') and float(splitUp['MEME Motif1 E-Value'])<=0.05 and (not splitUp['Up.MEME Motif1 Matches']=='NA'): # and len(splitUp[10].split(' '))<10:
            tfs[j]['MEME'] += splitUp['Up.MEME Motif1 Matches'].split(' ')
            if not splitUp['Up.MEME Motif1 Correlated Matches']=='NA':
                correlatedTfs[j]['MEME'] += [k.split(':')[0] for k in splitUp['Up.MEME Motif1 Correlated Matches'].split(' ')]
                for k in splitUp['Up.MEME Motif1 Correlated Matches'].split(' '):
                    corTF1 = k.split(':')
                    if float(corTF1[1])>0 and (float(splitUp['OS.covAgeSex.p'])<=0.05 and float(splitUp['mesoTCGA_OS.age.sex.p'])<=0.05 and ((float(splitUp['OS.covAgeSex'])/abs(float(splitUp['OS.covAgeSex'])))==(float(splitUp['mesoTCGA_OS.age.sex'])/abs(float(splitUp['mesoTCGA_OS.age.sex']))))) and float(splitUp['OS.covAgeSex'])>0:
                        siRNATFs[j].append({'tf':corTF1[0], 'tf.method':'MEME_motif1', 'tf.percTargets':'NA', 'tf_bic.R':corTF1[1], 'tf_bic.pv':corTF1[2], 'bic_surv.coef':splitUp['OS.covAgeSex'], 'bic_surv.pv':splitUp['OS.covAgeSex.p'], 'genes':splitUp['Genes']})
                bestTfs[j]['MEME'] += [splitUp['Up.MEME Motif1 Minimum Correlated'].split(':')[0]]
        # MEME2
        if (not splitUp['Up.MEME Motif2 E-Value']=='NA') and float(splitUp['Up.MEME Motif2 E-Value'])<=0.05 and (not splitUp['Up.MEME Motif2 Matches']=='NA'): #  and len(splitUp[14].split(' '))<10:
            tfs[j]['MEME'] += splitUp['Up.MEME Motif2 Matches'].split(' ')
            if not splitUp['Up.MEME Motif2 Correlated Matches']=='NA':
                correlatedTfs[j]['MEME'] += [k.split(':')[0] for k in splitUp['Up.MEME Motif2 Correlated Matches'].split(' ')]
                for k in splitUp['Up.MEME Motif2 Correlated Matches'].split(' '):
                    corTF1 = k.split(':')
                    if float(corTF1[1])>0 and (float(splitUp['OS.covAgeSex.p'])<=0.05 and float(splitUp['mesoTCGA_OS.age.sex.p'])<=0.05 and ((float(splitUp['OS.covAgeSex'])/abs(float(splitUp['OS.covAgeSex'])))==(float(splitUp['mesoTCGA_OS.age.sex'])/abs(float(splitUp['mesoTCGA_OS.age.sex']))))) and float(splitUp['OS.covAgeSex'])>0:
                        siRNATFs[j].append({'tf':corTF1[0], 'tf.method':'MEME_motif2', 'tf.percTargets':'NA', 'tf_bic.R':corTF1[1], 'tf_bic.pv':corTF1[2], 'bic_surv.coef':splitUp['OS.covAgeSex'], 'bic_surv.pv':splitUp['OS.covAgeSex.p'], 'genes':splitUp['Genes']})
                bestTfs[j]['MEME'] += [splitUp['Up.MEME Motif2 Minimum Correlated'].split(':')[0]]
        # Weeder
        if (not splitUp['Up.WEEDER Motif1 Matches']=='NA'): #  and len(splitUp[17].split(' '))<10:
            tfs[j]['WEEDER'] += splitUp['Up.WEEDER Motif1 Matches'].split(' ')
            if not splitUp['Up.WEEDER Motif1 Correlated Matches']=='NA':
                correlatedTfs[j]['WEEDER'] += [k.split(':')[0] for k in splitUp['Up.WEEDER Motif1 Correlated Matches'].split(' ')]
                for k in splitUp['Up.WEEDER Motif1 Correlated Matches'].split(' '):
                    corTF1 = k.split(':')
                    if float(corTF1[1])>0 and (float(splitUp['OS.covAgeSex.p'])<=0.05 and float(splitUp['mesoTCGA_OS.age.sex.p'])<=0.05 and ((float(splitUp['OS.covAgeSex'])/abs(float(splitUp['OS.covAgeSex'])))==(float(splitUp['mesoTCGA_OS.age.sex'])/abs(float(splitUp['mesoTCGA_OS.age.sex']))))) and float(splitUp['OS.covAgeSex'])>0:
                        siRNATFs[j].append({'tf':corTF1[0], 'tf.method':'WEEDER_motif1', 'tf.percTargets':'NA', 'tf_bic.R':corTF1[1], 'tf_bic.pv':corTF1[2], 'bic_surv.coef':splitUp['OS.covAgeSex'], 'bic_surv.pv':splitUp['OS.covAgeSex.p'], 'genes':splitUp['Genes']})
                bestTfs[j]['WEEDER'] += [splitUp['Up.WEEDER Motif1 Minimum Correlated'].split(':')[0]]
        # Weeder
        if (not splitUp['Up.WEEDER Motif2 Matches']=='NA'): #  and len(splitUp[20].split(' '))<10:
            tfs[j]['WEEDER'] += splitUp['Up.WEEDER Motif2 Matches'].split(' ')
            if not splitUp['Up.WEEDER Motif2 Correlated Matches']=='NA':
                correlatedTfs[j]['WEEDER'] += [k.split(':')[0] for k in splitUp['Up.WEEDER Motif2 Correlated Matches'].split(' ')]
                for k in splitUp['Up.WEEDER Motif2 Correlated Matches'].split(' '):
                    corTF1 = k.split(':')
                    if float(corTF1[1])>0 and (float(splitUp['OS.covAgeSex.p'])<=0.05 and float(splitUp['mesoTCGA_OS.age.sex.p'])<=0.05 and ((float(splitUp['OS.covAgeSex'])/abs(float(splitUp['OS.covAgeSex'])))==(float(splitUp['mesoTCGA_OS.age.sex'])/abs(float(splitUp['mesoTCGA_OS.age.sex']))))) and float(splitUp['OS.covAgeSex'])>0:
                        siRNATFs[j].append({'tf':corTF1[0], 'tf.method':'MEME_motif2', 'tf.percTargets':'NA', 'tf_bic.R':corTF1[1], 'tf_bic.pv':corTF1[2], 'bic_surv.coef':splitUp['OS.covAgeSex'], 'bic_surv.pv':splitUp['OS.covAgeSex.p'], 'genes':splitUp['Genes']})
                bestTfs[j]['WEEDER'] += [splitUp['Up.WEEDER Motif2 Minimum Correlated'].split(':')[0]]
        # TFBS_DB
        if (not splitUp['TFBS_DB.percTargets']=='NA') and float((splitUp['TFBS_DB.percTargets'].split(' '))[0])>=0.1 and float((splitUp['TFBS_DB.pValue'].split(' '))[0])<=0.05:
            tfs[j]['TFBS_DB'] += splitUp['TFBS_DB.TFs'].split(' ')
            if not splitUp['TFBS_DB.Correlated Matches']=='NA':
                correlatedTfs[j]['TFBS_DB'] += [k.split(':')[0] for k in splitUp['TFBS_DB.Correlated Matches'].split(' ')]
                for k in splitUp['TFBS_DB.Correlated Matches'].split(' '): # TODO HELLO !!!!! DUMPS OUT FINAL FILE!!!!!!!!!
                    corTF1 = k.split(':')
                    if float(corTF1[1])>0 and (float(splitUp['OS.covAgeSex.p'])<=0.05 and float(splitUp['mesoTCGA_OS.age.sex.p'])<=0.05 and ((float(splitUp['OS.covAgeSex'])/abs(float(splitUp['OS.covAgeSex'])))==(float(splitUp['mesoTCGA_OS.age.sex'])/abs(float(splitUp['mesoTCGA_OS.age.sex']))))) and float(splitUp['OS.covAgeSex'])>0:
                        siRNATFs[j].append({'tf':corTF1[0], 'tf.method':'TFBS_DB', 'tf.percTargets':splitUp['TFBS_DB.percTargets'].split(' ')[0], 'tf_bic.R':corTF1[1], 'tf_bic.pv':corTF1[2], 'bic_surv.coef':splitUp['OS.covAgeSex'], 'bic_surv.pv':splitUp['OS.covAgeSex.p'], 'genes':splitUp['Genes']})
                bestTfs[j]['TFBS_DB'] += [splitUp['TFBS_DB.Minimum Correlated'].split(':')[0]]
        # Correspondent TFs
        tmpCorrespondent = []
        if j in biMutAndTf:
            for tf1 in list(set(correlatedTfs[j]['MEME']+correlatedTfs[j]['WEEDER']+correlatedTfs[j]['TFBS_DB'])):
                for tf2 in biMutAndTf[j]:
                    if tf1==tf2['tf'] and not tf2['tf'] in tmpCorrespondent:
                        tmpCorrespondent.append(tf2['tf'])
        if j in biMutAndTf_rep:
            for tf1 in list(set(correlatedTfs[j]['MEME']+correlatedTfs[j]['WEEDER']+correlatedTfs[j]['TFBS_DB'])):
                for tf2 in biMutAndTf_rep[j]:
                    if tf1==tf2['tf'] and not tf2['tf'] in tmpCorrespondent:
                        tmpCorrespondent.append(tf2['tf'])
        #if (not splitUp['Correspondent.TFs']=='NA'):
        #    tmpCorrespondent += splitUp['Correspondent.TFs'].split(' ')
        if not len(tmpCorrespondent)==0:
            correspondentTfs[j] = sorted(list(set(tmpCorrespondent)))
            dumpBics.append(j)
        #correspondentMiRNAs[j] = []
        tmpCorrespondent = []
        # WEEDER:miRvestigator miRNA?
        if splitUp['3pUTR.WEEDER Motif1 Model']=='8mer':
            miRNAs[j]['miRvestigator'] += splitUp['3pUTR.WEEDER Motif1 Matches'].split(' ')
            for miR1 in splitUp['3pUTR.WEEDER Motif1 Matches'].split(' '):
                #print miR1
                if j in causalMiRNAs_rep:
                    if miR1 in causalMiRNAs_rep[j] and not miR1 in tmpCorrespondent:
                        #print j, miR1, causalMiRNAs_rep[j], 'WEEDER Motif1'
                        tmpCorrespondent.append(miR1)
                        dumpBics.append(j)
                if j in cor_miR and miR1 in cor_miR[j]:
                    if cor_miR[j][miR1]<0 and pv_miR[j][miR1]<=0.05:
                        #print j, miR1, cor_miR[j][miR1], pv_miR[j][miR1],'targetscan'
                        cor_miRNAs[j]['miRvestigator'].append(miR1)
                        dumpBics.append(j)
                        if (float(splitUp['OS.covAgeSex.p'])<=0.05 and float(splitUp['mesoTCGA_OS.age.sex.p'])<=0.05 and ((float(splitUp['OS.covAgeSex'])/abs(float(splitUp['OS.covAgeSex'])))==(float(splitUp['mesoTCGA_OS.age.sex'])/abs(float(splitUp['mesoTCGA_OS.age.sex']))))) and float(splitUp['OS.covAgeSex'])>0:
                            miRNAmimics[j].append({'miRNA':miR1, 'miRNA.method':'miRvestigator_Motif1', 'miRNA.percTargets':'NA', 'miR_bic.R':str(cor_miR[j][miR1]), 'miR_bic.pv':str(pv_miR[j][miR1]), 'bic_surv.coef':splitUp['OS.covAgeSex'], 'bic_surv.pv':splitUp['OS.covAgeSex.p'], 'genes':splitUp['Genes']})
        if splitUp['3pUTR.WEEDER Motif2 Model']=='8mer':
            miRNAs[j]['miRvestigator'] += splitUp['3pUTR.WEEDER Motif2 Matches'].split(' ')
            for miR1 in splitUp['3pUTR.WEEDER Motif2 Matches'].split(' '):
                if j in causalMiRNAs_rep:
                    if miR1 in causalMiRNAs_rep[j] and not miR1 in tmpCorrespondent:
                        #print j, miR1, causalMiRNAs_rep[j],'WEEDER Motif2'
                        tmpCorrespondent.append(miR1)
                        dumpBics.append(j)
                if j in cor_miR and miR1 in cor_miR[j]:
                    if cor_miR[j][miR1]<0 and pv_miR[j][miR1]<=0.05:
                        #print j, miR1, cor_miR[j][miR1], pv_miR[j][miR1],'targetscan'
                        cor_miRNAs[j]['miRvestigator'].append(miR1)
                        dumpBics.append(j)
                        if (float(splitUp['OS.covAgeSex.p'])<=0.05 and float(splitUp['mesoTCGA_OS.age.sex.p'])<=0.05 and ((float(splitUp['OS.covAgeSex'])/abs(float(splitUp['OS.covAgeSex'])))==(float(splitUp['mesoTCGA_OS.age.sex'])/abs(float(splitUp['mesoTCGA_OS.age.sex']))))) and float(splitUp['OS.covAgeSex'])>0:
                            miRNAmimics[j].append({'miRNA':miR1, 'miRNA.method':'miRvestigator_Motif2', 'miRNA.percTargets':'NA', 'miR_bic.R':str(cor_miR[j][miR1]), 'miR_bic.pv':str(pv_miR[j][miR1]), 'bic_surv.coef':splitUp['OS.covAgeSex'], 'bic_surv.pv':splitUp['OS.covAgeSex.p'], 'genes':splitUp['Genes']})
        # PITA miRNA?
        if (not splitUp['3pUTR_pita.percTargets']=='NA') and float((splitUp['3pUTR_pita.percTargets'].split(' '))[0])>=0.1 and float((splitUp['3pUTR_pita.pValue'].split(' '))[0])<=0.05: # and len(splitUp[36].split(' '))<10:
            miRNAs[j]['pita'] += splitUp['3pUTR_pita.miRNAs'].split(' ')
            for miR1 in splitUp['3pUTR_pita.miRNAs'].split(' '):
                if j in causalMiRNAs_rep:
                    if miR1 in causalMiRNAs_rep[j] and not miR1 in tmpCorrespondent:
                        #print j, miR1, causalMiRNAs_rep[j],'pita'
                        tmpCorrespondent.append(miR1)
                        dumpBics.append(j)
                if j in cor_miR and miR1 in cor_miR[j]:
                    if cor_miR[j][miR1]<0 and pv_miR[j][miR1]<=0.05:
                        #print j, miR1, cor_miR[j][miR1], pv_miR[j][miR1],'targetscan'
                        cor_miRNAs[j]['pita'].append(miR1)
                        dumpBics.append(j)
                        if (float(splitUp['OS.covAgeSex.p'])<=0.05 and float(splitUp['mesoTCGA_OS.age.sex.p'])<=0.05 and ((float(splitUp['OS.covAgeSex'])/abs(float(splitUp['OS.covAgeSex'])))==(float(splitUp['mesoTCGA_OS.age.sex'])/abs(float(splitUp['mesoTCGA_OS.age.sex']))))) and float(splitUp['OS.covAgeSex'])>0:
                            miRNAmimics[j].append({'miRNA':miR1, 'miRNA.method':'pita', 'miRNA.percTargets':splitUp['3pUTR_pita.percTargets'].split(' ')[0], 'miR_bic.R':str(cor_miR[j][miR1]), 'miR_bic.pv':str(pv_miR[j][miR1]), 'bic_surv.coef':splitUp['OS.covAgeSex'], 'bic_surv.pv':splitUp['OS.covAgeSex.p'], 'genes':splitUp['Genes']})
        # TargetScan miRNA?
        if (not splitUp['3pUTR_targetScan.percTargets']=='NA') and float((splitUp['3pUTR_targetScan.percTargets'].split(' '))[0])>=0.1 and float((splitUp['3pUTR_targetScan.pValue'].split(' '))[0])<=0.05: # and len(splitUp[39].split(' '))<10:
            miRNAs[j]['targetscan'] += splitUp['3pUTR_targetScan.miRNAs'].split(' ')
            for miR1 in splitUp['3pUTR_targetScan.miRNAs'].split(' '):
                if j in causalMiRNAs_rep:
                    if miR1 in causalMiRNAs_rep[j] and not miR1 in tmpCorrespondent:
                        #print j, miR1, causalMiRNAs_rep[j],'targetscan'
                        tmpCorrespondent.append(miR1)
                        dumpBics.append(j)
                if j in cor_miR and miR1 in cor_miR[j]:
                    if cor_miR[j][miR1]<0 and pv_miR[j][miR1]<=0.05:
                        #print j, miR1, cor_miR[j][miR1], pv_miR[j][miR1],'targetscan'
                        cor_miRNAs[j]['targetscan'].append(miR1)
                        dumpBics.append(j)
                        if (float(splitUp['OS.covAgeSex.p'])<=0.05 and float(splitUp['mesoTCGA_OS.age.sex.p'])<=0.05 and ((float(splitUp['OS.covAgeSex'])/abs(float(splitUp['OS.covAgeSex'])))==(float(splitUp['mesoTCGA_OS.age.sex'])/abs(float(splitUp['mesoTCGA_OS.age.sex']))))) and float(splitUp['OS.covAgeSex'])>0:
                            miRNAmimics[j].append({'miRNA':miR1, 'miRNA.method':'TargetScan', 'miRNA.percTargets':splitUp['3pUTR_targetScan.percTargets'].split(' ')[0], 'miR_bic.R':str(cor_miR[j][miR1]), 'miR_bic.pv':str(pv_miR[j][miR1]), 'bic_surv.coef':splitUp['OS.covAgeSex'], 'bic_surv.pv':splitUp['OS.covAgeSex.p'], 'genes':splitUp['Genes']})
        if not len(tmpCorrespondent)==0:
            correspondentMiRNAs[j] = list(set(tmpCorrespondent))
                    
        # GO Biological Processes
        if not splitUp['GO_Term_BP']=='NA':
            goBP[j] += splitUp['GO_Term_BP'].split(';')

        # Hallmarks of Cancer
        tmpHallmarks = getSignificantHallmarks(splitUp)
        for hallmark in tmpHallmarks:
            if not hallmark in hallmarksBi[j]:
                hallmarksBi[j].append(hallmark)
outFile = open('postProcessed_clustersOfBiclusters_CNA_CNVkit.csv','w')
outFile.write('\n'.join(writeMe))
outFile.close()


# Correspondent TFs TCGA
correspondentTfs_TCGA = {}
for bic1 in correlatedTfs:
    if bic1 in biMutAndTf_rep and bic1 in tfs:
        tmp = []
        for i in list(set(correlatedTfs[bic1]['MEME']+correlatedTfs[bic1]['WEEDER']+correlatedTfs[bic1]['TFBS_DB'])):
            for j in biMutAndTf_rep[bic1]:
                #print i, j
                if i==j['tf'] and not j in tmp:
                    #print bic1, i, j
                    tmp.append(j)
                #if i['tf']==j['tf']:
                #    print i,j
        if len(tmp)>0:
            correspondentTfs_TCGA[bic1] = tmp
#print 'correspondentTfs_TCGA',correspondentTfs_TCGA

# Dump out siRNATFs file
writeMe = ['Bicluster,TF_Entrez,TF_Symbol,TF_Method,TF_PercentTargets,OS.covAgeSex,OS.covAgeSex.p,TF_Bic.R,TF_Bic.p,Correspondent,Autoregulatory,Hallmarks,Bic_Genes']
for bic1 in siRNATFs:
    for tf1 in siRNATFs[bic1]:
        correspondent1 = 'No'
        if bic1 in correspondentTfs and tf1['tf'] in correspondentTfs[bic1]:
            correspondent1 = 'Yes'
        autoregulatory = 'No'
        if tf1['tf'] in tf1['genes'].split(' '):
            autoregulatory = 'Yes'
        if not tf1['tf.method']=='TFBS_DB':
            tf1['tf.percTargets'] = str(percentTargets(tf1['tf'], tf1['genes'].split(' '), tfbsDb, id2family, family2id, id2motif))
        writeMe.append(','.join([bic1, tf1['tf'], entrezId2gene[tf1['tf']], tf1['tf.method'], tf1['tf.percTargets'], tf1['bic_surv.coef'], tf1['bic_surv.pv'],  tf1['tf_bic.R'], tf1['tf_bic.pv'], correspondent1, autoregulatory, ' '.join(hallmarksBi[bic1]), tf1['genes']]))

with open('siRNATFs_10_17_2019.csv','w') as outFile:
    outFile.write('\n'.join(writeMe))


# Dump out siRNATFs file
writeMe = ['Bicluster,miRBase,miRNA,miR_Method,miR_PercentTargets,OS.covAgeSex,OS.covAgeSex.p,mIR_Bic.R,miR_Bic.p,Hallmarks,Bic_Genes']
for bic1 in miRNAmimics:
    for miR1 in miRNAmimics[bic1]:
        #correspondent1 = 'No'
        #if bic1 in correspondentTfs and miR1['tf'] in correspondentTfs[bic1]:
        #    correspondent1 = 'Yes'
        #if not tf1['tf.method']=='TFBS_DB':
        #    tf1['tf.percTargets'] = str(percentTargets(tf1['tf'], tf1['genes'].split(' '), tfbsDb, id2family, family2id, id2motif))
        writeMe.append(','.join([bic1, miR1['miRNA'], id2miR[miR1['miRNA']], miR1['miRNA.method'], miR1['miRNA.percTargets'], miR1['bic_surv.coef'], miR1['bic_surv.pv'],  miR1['miR_bic.R'], miR1['miR_bic.pv'], ' '.join(hallmarksBi[bic1]), miR1['genes']]))

with open('miRNAmimics_10_18_2019.csv','w') as outFile:
    outFile.write('\n'.join(writeMe))

# Causal TF in both Bueno and TCGA
tmp1 = [item for bic1 in biMutAndTf if bic1 in tfs for item in biMutAndTf[bic1]]
#missing = [i for i in tmp1 if not getMutName(i['mut'])[0] in entrezId2gene.keys()]
#tmp2 = entrezId2gene.keys()
#for j in tmp1:
#    if not (getMutName(j['mut'])[0] in tmp2):
#        print j
#        missing.append(getMutName(j['mut'])[0])
#print 'missing',missing


replicationCausality_mutReg = set([entrezId2gene[getMutName(i['mut'])[0]]+'_'+i['mut'].split('_')[1]+'|'+entrezId2gene[i['tf']] for i in tmp1])
tmp2 = [item for bic1 in biMutAndTf_rep if bic1 in tfs for item in biMutAndTf_rep[bic1]]
replicationCausality_mutReg_rep = set([entrezId2gene[getMutName(i['mut'])[0]]+'_'+getMutName(i['mut'])[1]+'|'+entrezId2gene[i['tf']] for i in tmp2])
common_mutReg2 = list(replicationCausality_mutReg.intersection(replicationCausality_mutReg_rep))
print 'common_mutReg2',common_mutReg2,'\n'

# Causal TF in both Bueno and TCGA & 
tmp3 = [item for bic1 in biMutAndTf if bic1 in tfs for item in biMutAndTf[bic1] if bic1 in correspondentTfs and item['tf'] in correspondentTfs[bic1]]
replicationCausality_mutReg_cor = set([entrezId2gene[getMutName(i['mut'])[0]]+'_'+getMutName(i['mut'])[1]+'|'+entrezId2gene[i['tf']] for i in tmp3])
tmp4 = [item for bic1 in biMutAndTf_rep if bic1 in tfs for item in biMutAndTf_rep[bic1] if bic1 in correspondentTfs_TCGA and item['tf'] in [item2['tf'] for item2 in correspondentTfs_TCGA[bic1]]]
replicationCausality_mutReg_rep_cor = set([entrezId2gene[getMutName(i['mut'])[0]]+'_'+getMutName(i['mut'])[1]+'|'+entrezId2gene[i['tf']] for i in tmp4])
common_mutReg3 = list(replicationCausality_mutReg_cor.intersection(replicationCausality_mutReg_rep_cor))
print 'common_mutReg3',common_mutReg3,'\n'

# Causality replication (mutReg) also correspondent for both
common_mutReg4 = list(set(list(replicationCausality_mutReg_cor.intersection(replicationCausality_mutReg_rep))+list(replicationCausality_mutReg.intersection(replicationCausality_mutReg_rep_cor))))
print 'common_mutReg4',common_mutReg4,'\n'

# Causality replication (full path)
replicatingCausality = {}
for bic1 in biMutAndTf:
    if bic1 in tfs and bic1 in biMutAndTf_rep:
        tmp = []
        for i in biMutAndTf[bic1]:
            for j in biMutAndTf_rep[bic1]:
                #print i, j
                if i['mut']==j['mut'] and i['tf']==j['tf']:
                    #print bic1, i, j
                    tmp.append(i)
                #if i['tf']==j['tf']:
                #    print i,j
        if len(tmp)>0:
            replicatingCausality[bic1] = tmp
#print 'replicatingCausality',replicatingCausality,'\n'


# Causality replication (full path) and correspondent
replicatingCausality2 = {}
for bic1 in biMutAndTf:
    if bic1 in tfs and bic1 in biMutAndTf_rep:
        tmp = []
        for i in biMutAndTf[bic1]:
            for j in biMutAndTf_rep[bic1]:
                #print i, j
                if i['mut']==j['mut'] and i['tf']==j['tf'] and ((bic1 in correspondentTfs and i['tf'] in correspondentTfs[bic1]) or (bic1 in correspondentTfs_TCGA and i['tf'] in correspondentTfs_TCGA[bic1])):
                    #print bic1, i, j
                    tmp.append(i)
                #if i['tf']==j['tf']:
                #    print i,j
        if len(tmp)>0:
            replicatingCausality2[bic1] = tmp
print 'replicatingCausality',replicatingCausality2,'\n'


# Dump out TFs to test
outFile = open('tf_regulators_CNA_CNVkit_8_13_2019.csv','w')
outFile.write('cluster,MEME,WEEDER,TFBS_DB,Number Predicted\n')
outFile.write('\n'.join([i+','+';'.join(list(set(tfs[i]['MEME'])))+','+';'.join(list(set(tfs[i]['WEEDER'])))+','+';'.join(list(set(tfs[i]['TFBS_DB'])))+','+str(len(list(set(list(set(tfs[i]['MEME']))+list(set(tfs[i]['WEEDER']))+list(set(tfs[i]['TFBS_DB'])))))) for i in tfs]))
outFile.close()

# Dump out TFs to test
outFile = open('correlated_tf_regulators_CNA_CNVkit_8_13_2019.csv','w')
outFile.write('cluster,MEME,WEEDER,TFBS_DB,Number Predicted\n')
outFile.write('\n'.join([i+','+';'.join(list(set(correlatedTfs[i]['MEME'])))+','+';'.join(list(set(correlatedTfs[i]['WEEDER'])))+','+';'.join(list(set(correlatedTfs[i]['TFBS_DB'])))+','+str(len(list(set(list(set(correlatedTfs[i]['MEME']))+list(set(correlatedTfs[i]['WEEDER']))+list(set(correlatedTfs[i]['TFBS_DB'])))))) for i in correlatedTfs]))
outFile.close()

# Dump out TFs to test
outFile = open('best_tf_regulators_CNA_CNVkit_8_13_2019.csv','w')
outFile.write('cluster,MEME,WEEDER,TFBS_DB,Number Predicted\n')
outFile.write('\n'.join([i+','+';'.join(list(set(bestTfs[i]['MEME'])))+','+';'.join(list(set(bestTfs[i]['WEEDER'])))+','+';'.join(list(set(bestTfs[i]['TFBS_DB'])))+','+str(len(list(set(list(set(bestTfs[i]['MEME']))+list(set(bestTfs[i]['WEEDER']))+list(set(bestTfs[i]['TFBS_DB'])))))) for i in bestTfs]))
outFile.close()

# Dump out miRNAs to test
outFile = open('miRNA_regulators_CNA_CNVkit_8_13_2019.csv','w')
outFile.write('Bicluster,miRvestigator,PITA,TargetScan\n')
outFile.write('\n'.join([i+','+' '.join(list(set(miRNAs[i]['miRvestigator'])))+','+' '.join(list(set(miRNAs[i]['pita'])))+','+' '.join(list(set(miRNAs[i]['targetscan']))) for i in miRNAs]))
outFile.close()

# Dump out miRNAs to test
outFile = open('correlated_miRNA_regulators_CNA_CNVkit_8_13_2019.csv','w')
outFile.write('Bicluster,miRvestigator,PITA,TargetScan\n')
outFile.write('\n'.join([i+','+' '.join(list(set(cor_miRNAs[i]['miRvestigator'])))+','+' '.join(list(set(cor_miRNAs[i]['pita'])))+','+' '.join(list(set(cor_miRNAs[i]['targetscan']))) for i in miRNAs]))
outFile.close()

# Dump out GO BPs to test
outFile = open('GO_BPs_metaClusters_CNA_CNVkit_8_13_2019.csv','w')
outFile.write('\n'.join([i+','+';'.join(list(set(goBP[i])))+','+str(len(list(set(goBP[i])))) for i in goBP]))
outFile.close()


# SIF file output
pams = []
LoFs = []
Acts = []
CNAdels = []
bics = []
tfs = []
miRNAs = []
miRNACausalFlows = []
miRNACorrelated = []
sifWriteMe = []
attWriteMe = []
cleanNamesWriteMe = []
lociMapping = {}
for bic1 in list(set(dumpBics)):
    if bic1 in correspondentTfs:
        for tf1 in correspondentTfs[bic1]:
            tf2 = entrezId2gene[tf1]
            if not tf2 in tfs:
                tfs.append(tf2)
            if not bic1 in bics:
                bics.append(bic1)
            sifWriteMe.append(tf2+' r2b '+bic1)
            attWriteMe.append(tf2+' TF')
            attWriteMe.append(bic1+' Bicluster')
            cleanNamesWriteMe.append(tf2+' '+tf2)
            if bic1 in biMutAndTf:
                for pair1 in biMutAndTf[bic1]:
                    if pair1['tf']==tf1:
                        tmp = getMutName(pair1['mut'])
                        #print(tmp)
                        if tmp[1]=='PAM' and not tmp[0] in pams:
                            pams.append(tmp[0])
                        if tmp[1]=='CNAdel' and not tmp[0] in CNAdels:
                            CNAdels.append(tmp[0])
                        if tmp[1]=='LoF' and not tmp[0] in LoFs:
                            LoFs.append(tmp[0])
                        if tmp[1]=='Act' and not tmp[0] in Acts:
                            Acts.append(tmp[0])
                        #print(tmp)
                        type1 = tmp[1]
                        if len(tmp)>1:
                            if tmp[0] in entrezId2gene:
                                mut1 = entrezId2gene[tmp[0]]+'_'+tmp[1]
                                loc1 = mapMut2Loci(pair1['mut'], loci_Bueno)
                                if loc1:
                                    if not loc1 in lociMapping:
                                        lociMapping[loc1] = []
                                    if not mut1 in lociMapping[loc1]:
                                        lociMapping[loc1].append(mut1)
                                cleanNamesWriteMe.append(mut1+' '+entrezId2gene[tmp[0]])
                            elif len(tmp)==3:
                                mut1 = pair1['mut']
                                type1 = tmp[2]
                                cleanNamesWriteMe.append(mut1+' '+tmp[0]+'_'+tmp[1])
                            else:
                                print tmp[0]
                                mut1 = ''
                        else:
                            if pair1['mut'] in entrezId2gene:
                                mut1 = entrezId2gene[pair1['mut']]
                            else:
                                print pair1['mut']
                                mut1=''
                        sifWriteMe.append(mut1+' g2r '+tf2)
                        attWriteMe.append(mut1+' '+type1)            
    if bic1 in correspondentMiRNAs:
        for miR1 in correspondentMiRNAs[bic1]:
            miR2 = id2miR[miR1]
            if not miR2 in miRNAs:
                miRNAs.append(miR2)
            if not bic1 in bics:
                bics.append(bic1)
            sifWriteMe.append(miR2+' r2b '+bic1)
            attWriteMe.append(miR2+' miRNA')
            attWriteMe.append(bic1+' Bicluster')
            cleanNamesWriteMe.append(miR2+' '+miR2)
            if bic1 in biMutAndMiRNA_rep:
                for pair1 in biMutAndMiRNA_rep[bic1]:
                    if pair1['miRNA']==miR1:
                        tmp = getMutName(pair1['mut'])
                        #print(tmp)
                        if tmp[1]=='PAM' and not tmp[0] in pams:
                            pams.append(tmp[0])
                        if tmp[1]=='CNAdel' and not tmp[0] in CNAdels:
                            CNAdels.append(tmp[0])
                        if tmp[1]=='LoF' and not tmp[0] in LoFs:
                            LoFs.append(tmp[0])
                        if tmp[1]=='Act' and not tmp[0] in Acts:
                            Acts.append(tmp[0])
                        type1 = tmp[1]
                        if len(tmp)>1:
                            if tmp[0] in entrezId2gene:
                                mut1 = entrezId2gene[tmp[0]]+'_'+tmp[1]
                                cleanNamesWriteMe.append(mut1+' '+entrezId2gene[tmp[0]])
                            elif len(tmp)==3:
                                mut1 = pair1['mut']
                                cleanNamesWriteMe.append(mut1+' '+tmp[0]+'_'+tmp[1])
                                type1 = tmp[2]
                            else:
                                print tmp[0]
                                mut1 = ''
                        else:
                            if pair1['mut'] in entrezId2gene:
                                mut1 = entrezId2gene[pair1['mut']]
                            else:
                                print pair1['mut']
                                mut1=''
                        sifWriteMe.append(mut1+' g2r '+miR2)
                        attWriteMe.append(mut1+' '+type1)
                        miRNACausalFlows.append([mut1,miR2,bic1])
    if bic1 in cor_miRNAs:
        for set1 in cor_miRNAs[bic1]:
            for miR1 in cor_miRNAs[bic1][set1]:
                miR2 = id2miR[miR1]
                if not miR2 in miRNAs:
                    miRNAs.append(miR2)
                if not bic1 in bics:
                    bics.append(bic1)
                sifWriteMe.append(miR2+' r2b '+bic1)
                attWriteMe.append(miR2+' miRNA')
                attWriteMe.append(bic1+' Bicluster')
                cleanNamesWriteMe.append(miR2+' '+miR2)
                miRNACorrelated.append([miR1,miR2,bic1])
    if (bic1 in correspondentTfs) or (bic1 in cor_miRNAs):
        if bic1 in hallmarksBi:
            for h1 in hallmarksBi[bic1]:
                sifWriteMe.append(bic1+' b2h '+h1)
                attWriteMe.append(bic1+' Bicluster')
                attWriteMe.append(h1+' Hallmark')
                cleanNamesWriteMe.append(h1+' '+h1)
# Write them out
with open('causalAndMechanistic_network_CNA_CNVkit_8_13_2019.sif','w') as outFile:
    outFile.write('\n'.join(list(set(sifWriteMe))))
with open('causalAndMechanistic_network_CNA_CNVkit_8_13_2019.txt','w') as outFile:
    outFile.write('\n'.join(list(set(attWriteMe))))
with open('causalAndMechanistic_network_CNA_CNVkit_cleanNames_8_13_2019.txt','w') as outFile:
    outFile.write('\n'.join(list(set(cleanNamesWriteMe))))

# Write out for Chuong
with open('corMirs_7_20_2019.csv','w') as outFile:
    outFile.write('miRNA,Bicluster,R,p-value,Hallmarks')
    for pair1 in miRNACorrelated:
        tmp1 = ''
        tmp2 = ''
        if pair1[2] in cor_miR and pair1[0] in cor_miR[pair1[2]]:
            tmp1 = cor_miR[pair1[2]][pair1[0]]
            tmp2 = pv_miR[pair1[2]][pair1[0]]
        outFile.write('\n'+pair1[1]+','+pair1[2]+','+str(tmp1)+','+str(tmp2)+','+' '.join(hallmarksBi[pair1[2]]))
      



print('Mutations:')
print('\tPAMs:',len(pams))
print('\tCNAdels:',len(CNAdels))
print('\tLoFs:',len(LoFs))
print('\tActs:',len(Acts))
print('Regulators:')
print('\tTFs:',len(tfs))
print('\tmiRNAs:',len(miRNAs))
print('Biclusters:',len(bics))

"""
######################
### Plot causality ###
######################
# Load expression
gexp1 = pd.read_csv('../mesothelioma_norm.txt', header=0, index_col=0, sep=' ')

# Load mutations
somMut1 = pd.read_csv('../BUENO_MESO_finalMutFile_deep_filtered_mmf_0.05.csv', header=0, index_col=0)

# Load phenotypes
pheno1 = pd.read_csv('../phenotypes_meso_noFilter.csv', header=0, index_col=0)

# Load bicluster eigengenes
be1 = pd.read_csv('eigengenes/biclusterEigengenes_pita.csv', header=0, index_col=0)
be1.columns = [i.lstrip('X') for i in be1.columns]
be1.index = ['pita_'+str(i) for i in be1.index]
be2 = pd.read_csv('eigengenes/biclusterEigengenes_targetscan.csv', header=0, index_col=0)
be2.columns = [i.lstrip('X') for i in be2.columns]
be2.index = ['targetscan_'+str(i) for i in be2.index]
be3 = pd.read_csv('eigengenes/biclusterEigengenes_tfbs_db.csv', header=0, index_col=0)
be3.columns = [i.lstrip('X') for i in be3.columns]
be3.index = ['tfbs_db_'+str(i) for i in be3.index]

# Combined for plotting
dAll = pd.concat([gexp1,somMut1,be1,be2,be3,pheno1.T], sort=True).T

pp = PdfPages('causality_plots_Genentech.pdf')
for bic1 in list(set(dumpBics)):
    if bic1 in correspondentTfs:
        for tf1 in correspondentTfs[bic1]:
            tf2 = entrezId2gene[tf1]
            if not tf2 in tfs:
                tfs.append(tf2)
            if not bic1 in bics:
                bics.append(bic1)
            if bic1 in biMutAndTf:
                for pair1 in biMutAndTf[bic1]:
                    if pair1['tf']==tf1:
                        tmp = getMutName(pair1['mut'])
                        #print(tmp)
                        if tmp[1]=='PAM' and not tmp[0] in pams:
                            pams.append(tmp[0])
                        if tmp[1]=='CNAdel' and not tmp[0] in CNAdels:
                            CNAdels.append(tmp[0])
                        if tmp[1]=='LoF' and not tmp[0] in LoFs:
                            LoFs.append(tmp[0])
                        if tmp[1]=='Act' and not tmp[0] in Acts:
                            Acts.append(tmp[0])
                        #print(tmp)
                        type1 = tmp[1]
                        if len(tmp)>1:
                            if tmp[0] in entrezId2gene:
                                mut1 = entrezId2gene[tmp[0]]+'_'+tmp[1]
                                loc1 = mapMut2Loci(pair1['mut'], loci_Bueno)
                                if loc1:
                                    if not loc1 in lociMapping:
                                        lociMapping[loc1] = []
                                    if not mut1 in lociMapping[loc1]:
                                        lociMapping[loc1].append(mut1)
                            elif len(tmp)==3:
                                mut1 = pair1['mut']
                                type1 = tmp[2]
                            else:
                                print tmp[0]
                                mut1 = ''
                        else:
                            if pair1['mut'] in entrezId2gene:
                                mut1 = entrezId2gene[pair1['mut']]
                            else:
                                print pair1['mut']
                                mut1=''
                        # Plot
                        if pair1['mut'] in dAll.columns:
                            print(pair1['mut']+' -> '+tf2+' -> '+bic1)
                            fig, ax = plt.subplots(ncols = 2, nrows = 2)
                            plt.subplots_adjust(wspace=0.5,hspace=0.5)
                            # Top-left: mutation vs. regulator
                            tt1 = stats.ttest_ind(dAll[int(tf1)][dAll[pair1['mut']]==1], dAll[int(tf1)][dAll[pair1['mut']]==0])
                            ax1 = sns.boxplot( y = list(dAll[int(tf1)]), x = list(dAll[pair1['mut']]), ax = ax[0,0]) #, color = ['#99A3A4','#B03A2E'])
                            ax1.set(xlabel = mut1, ylabel = tf2, title = "T = %.2f" % tt1[0] + '; '+"pv = %.2E" % tt1[1], xticklabels=['WT','Mutated'])
                            # Top-right: Boxplot regulator by mutation status
                            tt2 = stats.ttest_ind(dAll[bic1][dAll[pair1['mut']]==1], dAll[bic1][dAll[pair1['mut']]==0])
                            ax2 =sns.boxplot( y = list(dAll[bic1]), x = list(dAll[pair1['mut']]), ax = ax[0,1])
                            ax2.set(xlabel = mut1, ylabel = bic1, title = "T = %.2f" % tt2[0] + '; '+"pv = %.2E" % tt2[1], xticklabels=['WT','Mutated'])
                            # Bottom-left: Scatterplot bicluster by regulator
                            r3 = stats.pearsonr(list(dAll.dropna()[int(tf1)]), list(dAll.dropna()[bic1]))
                            ax3 = sns.regplot( y = list(dAll[bic1]), x = list(dAll[int(tf1)]), ax = ax[1,0]) # , hue = 'histology_WHO'
                            #sns.jointplot(x = int(tf1), y=bic1, data = dAll, kind="reg", ax=ax[1,0])
                            ax3.set(xlabel = tf2, ylabel = bic1, title = "R = %.2f" % r3[0] + '; '+"pv = %.2E" % r3[1])
                            # Bottom-right: Resdisual plot
                            dTmp = pd.DataFrame({ 'bic1':dAll[bic1],'tf1':dAll[int(tf1)]}).dropna()
                            lm1 = LinearRegression().fit(dTmp[['tf1']], dTmp['bic1'])
                            resid1 = dTmp['bic1']-lm1.predict(dTmp[['tf1']])
                            tt4 = stats.ttest_ind(resid1[dAll[pair1['mut']]==1], resid1[dAll[pair1['mut']]==0])
                            ax4 = sns.boxplot( y = list(resid1[dAll.dropna()[pair1['mut']].index]), x = list(dAll.dropna()[pair1['mut']]), ax = ax[1,1])
                            ax4.set(xlabel = mut1, ylabel = bic1, title = "T = %.2f" % tt4[0] + '; '+"pv = %.2E" % tt4[1], xticklabels=['WT','Mutated'])
                            pp.savefig()
                            plt.close()
pp.close()
"""
"""    if bic1 in correspondentMiRNAs:
        for miR1 in correspondentMiRNAs[bic1]:
            miR2 = id2miR[miR1]
            if not miR2 in miRNAs:
                miRNAs.append(miR2)
            if not bic1 in bics:
                bics.append(bic1)
            sifWriteMe.append(miR2+' r2b '+bic1)
            attWriteMe.append(miR2+' miRNA')
            attWriteMe.append(bic1+' Bicluster')
            cleanNamesWriteMe.append(miR2+' '+miR2)
            if bic1 in biMutAndMiRNA_rep:
                for pair1 in biMutAndMiRNA_rep[bic1]:
                    if pair1['miRNA']==miR1:
                        tmp = getMutName(pair1['mut'])
                        #print(tmp)
                        if tmp[1]=='PAM' and not tmp[0] in pams:
                            pams.append(tmp[0])
                        if tmp[1]=='CNAdel' and not tmp[0] in CNAdels:
                            CNAdels.append(tmp[0])
                        if tmp[1]=='LoF' and not tmp[0] in LoFs:
                            LoFs.append(tmp[0])
                        if tmp[1]=='Act' and not tmp[0] in Acts:
                            Acts.append(tmp[0])
                        type1 = tmp[1]
                        if len(tmp)>1:
                            if tmp[0] in entrezId2gene:
                                mut1 = entrezId2gene[tmp[0]]+'_'+tmp[1]
                                cleanNamesWriteMe.append(mut1+' '+entrezId2gene[tmp[0]])
                            elif len(tmp)==3:
                                mut1 = pair1['mut']
                                cleanNamesWriteMe.append(mut1+' '+tmp[0]+'_'+tmp[1])
                                type1 = tmp[2]
                            else:
                                print tmp[0]
                                mut1 = ''
                        else:
                            if pair1['mut'] in entrezId2gene:
                                mut1 = entrezId2gene[pair1['mut']]
                            else:
                                print pair1['mut']
                                mut1=''
                        sifWriteMe.append(mut1+' g2r '+miR2)
                        attWriteMe.append(mut1+' '+type1)
    if bic1 in cor_miRNAs:
        for set1 in cor_miRNAs[bic1]:
            for miR1 in cor_miRNAs[bic1][set1]:
                miR1 = id2miR[miR1]
                if not miR1 in miRNAs:
                    miRNAs.append(miR1)
                if not bic1 in bics:
                    bics.append(bic1)
                sifWriteMe.append(miR1+' r2b '+bic1)
                attWriteMe.append(miR1+' miRNA')
                attWriteMe.append(bic1+' Bicluster')
                cleanNamesWriteMe.append(miR1+' '+miR1)
    if (bic1 in correspondentTfs) or (bic1 in cor_miRNAs):
        if bic1 in hallmarksBi:
            for h1 in hallmarksBi[bic1]:
                sifWriteMe.append(bic1+' b2h '+h1)
                attWriteMe.append(bic1+' Bicluster')
                attWriteMe.append(h1+' Hallmark')
                cleanNamesWriteMe.append(h1+' '+h1)
 """