# coding: utf-8
from sqlalchemy import BINARY, Column, Float, ForeignKey, Integer, String, Text
from sqlalchemy.orm import relationship
from flask import Flask
from flask_sqlalchemy import SQLAlchemy

###########################
## Open database session ##
###########################
app = Flask(__name__)
#app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+mysqlconnector://root:root@mpm_web_db_1:3306/mpm_web'
app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+mysqlconnector://root:root@localhost:32000/mpm_web'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = 'FALSE'

db = SQLAlchemy(app)

# used
class BicGene(db.Model):
    __tablename__ = 'bic_gene'

    id = db.Column(db.Integer, primary_key=True)
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    gene_id = db.Column(db.ForeignKey('gene.id'), index=True)

    bicluster = db.relationship('Bicluster', primaryjoin='BicGene.bicluster_id == Bicluster.id', backref='bic_genes')
    gene = db.relationship('Gene', primaryjoin='BicGene.gene_id == Gene.id', backref='bic_genes')

# used
class BicGo(db.Model):
    __tablename__ = 'bic_go'

    id = db.Column(db.Integer, primary_key=True)
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    go_bp_id = db.Column(db.ForeignKey('go_bp.id'), index=True)

    bicluster = db.relationship('Bicluster', primaryjoin='BicGo.bicluster_id == Bicluster.id', backref='bic_goes')
    go_bp = db.relationship('GoBp', primaryjoin='BicGo.go_bp_id == GoBp.id', backref='bic_goes')

# used
class BicHal(db.Model):
    __tablename__ = 'bic_hal'

    id = db.Column(db.Integer, primary_key=True)
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    hallmark_id = db.Column(db.ForeignKey('hallmark.id'), index=True)

    bicluster = db.relationship('Bicluster', primaryjoin='BicHal.bicluster_id == Bicluster.id', backref='bic_hals')
    hallmark = db.relationship('Hallmark', primaryjoin='BicHal.hallmark_id == Hallmark.id', backref='bic_hals')

# used
class BicPat(db.Model):
    __tablename__ = 'bic_pat'

    id = db.Column(db.Integer, primary_key=True)
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    patient_id = db.Column(db.ForeignKey('patient.id'), index=True)

    bicluster = db.relationship('Bicluster', primaryjoin='BicPat.bicluster_id == Bicluster.id', backref='bic_pats')
    patient = db.relationship('Patient', primaryjoin='BicPat.patient_id == Patient.id', backref='bic_pats')

# used
class Bicluster(db.Model):
    __tablename__ = 'bicluster'

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(20), index=True)
    var_exp_fpc = db.Column(db.Float)
    var_exp_fpc_p_value = db.Column(db.Float)
    survival = db.Column(db.Float)
    survival_p_value = db.Column(db.Float)

# used
class CausalFlow(db.Model):
    __tablename__ = 'causal_flow'

    id = db.Column(db.Integer, primary_key=True)
    somatic_mutation_id = db.Column(db.ForeignKey('somatic_mutation.id'), index=True) # links through to mutation
    regulator_id = db.Column(db.Integer) # should this link to a table? miRNA table id or TFRegualtor table id. union types?
    regulator_type = db.Column(db.String(10)) # what is this? miRNA or Transcription Factor (table name!!!). should be "tf" or "mirna"
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    leo_nb_atob = db.Column(db.Float)
    mlogp_m_atob = db.Column(db.Float)

    bicluster = db.relationship('Bicluster', primaryjoin='CausalFlow.bicluster_id == Bicluster.id', backref='causal_flows')
    somatic_mutation = db.relationship('SomaticMutation', primaryjoin='CausalFlow.somatic_mutation_id == SomaticMutation.id', backref='causal_flows')

# used
class ExpDataset(db.Model):
    __tablename__ = 'exp_dataset'

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(20))
    type = db.Column(db.String(10))

# used
class Gene(db.Model):
    __tablename__ = 'gene'

    id = db.Column(db.Integer, primary_key=True)
    symbol = db.Column(db.String(100), index=True)
    entrez = db.Column(db.Integer)

# used
class GenePrior(db.Model):
    __tablename__ = 'gene_prior'

    id = db.Column(db.Integer, primary_key=True)
    gene_id = db.Column(db.ForeignKey('gene.id'), index=True)
    source = db.Column(db.String(11), index=True)
    description = db.Column(db.String(100), index=True)
    pmid = db.Column(db.Integer)

# used
class GeneExpression(db.Model):
    __tablename__ = 'gene_expression'

    id = db.Column(db.Integer, primary_key=True)
    exp_dataset_id = db.Column(db.ForeignKey('exp_dataset.id'), index=True)
    patient_id = db.Column(db.ForeignKey('patient.id'), index=True)
    gene_id = db.Column(db.ForeignKey('gene.id'), index=True)
    value = db.Column(db.Float)

    exp_dataset = db.relationship('ExpDataset', primaryjoin='GeneExpression.exp_dataset_id == ExpDataset.id', backref='gene_expressions')
    gene = db.relationship('Gene', primaryjoin='GeneExpression.gene_id == Gene.id', backref='gene_expressions')
    patient = db.relationship('Patient', primaryjoin='GeneExpression.patient_id == Patient.id', backref='gene_expressions')

# used
class GoBp(db.Model):
    __tablename__ = 'go_bp'

    id = db.Column(db.Integer, primary_key=True)
    go_id = db.Column(db.String(15), index=True)
    name = db.Column(db.Text)
    description = db.Column(db.Text)

# used
class GoGene(db.Model):
    __tablename__ = 'go_gene'

    id = db.Column(db.Integer, primary_key=True)
    go_bp_id = db.Column(db.ForeignKey('go_bp.id'), index=True)
    gene_id = db.Column(db.ForeignKey('gene.id'), index=True)

    gene = db.relationship('Gene', primaryjoin='GoGene.gene_id == Gene.id', backref='go_genes')
    go_bp = db.relationship('GoBp', primaryjoin='GoGene.go_bp_id == GoBp.id', backref='go_genes')

# used
class Hallmark(db.Model):
    __tablename__ = 'hallmark'

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100), index=True)

# used
class Mirna(db.Model):
    __tablename__ = 'mirna'

    id = db.Column(db.Integer, primary_key=True)
    mature_seq_id = db.Column(db.String(20)) # MIMAT[0-9]+
    name = db.Column(db.String(50), index=True)

# used
class MirnaPrior(db.Model):
    __tablename__ = 'mirna_prior'

    id = db.Column(db.Integer, primary_key=True)
    mirna_id = db.Column(db.ForeignKey('mirna.id'), index=True)
    source = db.Column(db.String(11), index=True)
    description = db.Column(db.String(100))
    pmid = db.Column(db.Integer)

# not used
class MirnaExpression(db.Model):
    __tablename__ = 'mirna_expression'

    id = db.Column(db.Integer, primary_key=True)
    exp_dataset_id = db.Column(db.ForeignKey('exp_dataset.id'), index=True)
    patient_id = db.Column(db.ForeignKey('patient.id'), index=True)
    mirna_id = db.Column(db.ForeignKey('mirna.id'), index=True) #   MIMAT[0-9]+
    value = db.Column(db.Float)

    exp_dataset = db.relationship('ExpDataset', primaryjoin='MirnaExpression.exp_dataset_id == ExpDataset.id', backref='mirna_expressions')
    mirna = db.relationship('Mirna', primaryjoin='MirnaExpression.mirna_id == Mirna.id', backref='mirna_expressions')
    patient = db.relationship('Patient', primaryjoin='MirnaExpression.patient_id == Patient.id', backref='mirna_expressions')

# used
class MirnaRegulator(db.Model):
    __tablename__ = 'mirna_regulator'

    id = db.Column(db.Integer, primary_key=True)
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    mirna_id = db.Column(db.ForeignKey('mirna.id'), index=True) #   MIMAT[0-9]+

    bicluster = db.relationship('Bicluster', primaryjoin='MirnaRegulator.bicluster_id == Bicluster.id', backref='mirna_regulators')
    mirna = db.relationship('Mirna', primaryjoin='MirnaRegulator.mirna_id == Mirna.id', backref='mirna_regulators')

# used
class MirnaTarget(db.Model):
    __tablename__ = 'mirna_targets'

    id = db.Column(db.Integer, primary_key=True)
    mirna_id = db.Column(db.ForeignKey('mirna.id'), index=True)
    gene_id = db.Column(db.ForeignKey('gene.id'), index=True)
    source = db.Column(db.String(11)) # either PITA or Target Scan, based on .json filename

    gene = db.relationship('Gene', primaryjoin='MirnaTarget.gene_id == Gene.id', backref='mirna_targets')
    mirna = db.relationship('Mirna', primaryjoin='MirnaTarget.mirna_id == Mirna.id', backref='mirna_targets')

# used
class Patient(db.Model):
    __tablename__ = 'patient'

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(50), index=True)

# used
class PhenoDatum(db.Model):
    __tablename__ = 'pheno_data'

    id = db.Column(db.Integer, primary_key=True)
    phenotype_id = db.Column(db.ForeignKey('phenotype.id'), index=True)
    patient_id = db.Column(db.ForeignKey('patient.id'), index=True)
    phenotype_value = db.Column(db.Float) # stores the value from a cell in the phenotype (like the halmark data type)
    phenotype_string = db.Column(db.String(50))
    # name = db.Column(db.String(50), nullable=False)
    # long_name = db.Column(db.String(50))

    patient = db.relationship('Patient', primaryjoin='PhenoDatum.patient_id == Patient.id', backref='pheno_data')
    phenotype = db.relationship('Phenotype', primaryjoin='PhenoDatum.phenotype_id == Phenotype.id', backref='pheno_data')

# used
class Phenotype(db.Model):
    __tablename__ = 'phenotype'

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(50), nullable=False) # column names (used for a relationship in PhenoDatum)
    long_name = db.Column(db.String(50)) # same as name

# used
class Replication(db.Model):
    __tablename__ = 'replication'

    id = db.Column(db.Integer, primary_key=True)
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    study = db.Column(db.String(50)) # "MESO TCGA"
    var_exp_fpc = db.Column(db.Float) # mesoTCGA_var.exp
    var_exp_fpc_p_value = db.Column(db.Float) # mesoTCGA_pc1.perm.p
    survival = db.Column(db.Float) # mesoTCGA_OS.age.sex
    survival_p_value = db.Column(db.Float) # mesoTCGA_OS.age.sex.p

    bicluster = db.relationship('Bicluster', primaryjoin='Replication.bicluster_id == Bicluster.id', backref='replications')

# used
class SomaticMutation(db.Model):
    __tablename__ = 'somatic_mutation'

    id = db.Column(db.Integer, primary_key=True)
    ext_id = db.Column(db.ForeignKey('gene.id'), index=True) # should be a gene id (link to gene table) (if the entrez id is the same as a gene we already have in the gene table, then link that to the gene table)
    mutation_type = db.Column(db.String(10)) # LoF, PAM, etc
    locus_id = db.Column(db.ForeignKey('locus.id'), index=True)

    locus = db.relationship('Locus', primaryjoin='SomaticMutation.locus_id == Locus.id', backref='somatic_mutations')

    mutation_name = ""

class BuenoDeepFilter(db.Model):
    __tablename__ = 'bueno_deep_filter'

    id = db.Column(db.Integer, primary_key=True)
    somatic_mutation_id = db.Column(db.ForeignKey('somatic_mutation.id'), index=True)
    patient_id = db.Column(db.ForeignKey('patient.id'), index=True)
    value = db.Column(db.Boolean)

    somatic_mutation = db.relationship('SomaticMutation', primaryjoin='BuenoDeepFilter.somatic_mutation_id == SomaticMutation.id', backref='bueno_deep_filters')
    patient = db.relationship('Patient', primaryjoin='BuenoDeepFilter.patient_id == Patient.id', backref='bueno_deep_filters')

# used
class TfRegulator(db.Model):
    __tablename__ = 'tf_regulator'

    id = db.Column(db.Integer, primary_key=True)
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    tf_id = db.Column(db.ForeignKey('gene.id'), index=True)
    r_value = db.Column(db.Float(precision='24,24'))
    p_value = db.Column(db.Float(precision='24,24'))
	# action = db.Column(db.String(15))

    bicluster = db.relationship('Bicluster', primaryjoin='TfRegulator.bicluster_id == Bicluster.id', backref='tf_regulators')
    gene = db.relationship('Gene', primaryjoin='TfRegulator.tf_id == Gene.id', backref='tf_regulators')

# used
class TfTarget(db.Model):
    __tablename__ = 'tf_targets'

    id = db.Column(db.Integer, primary_key=True)
    tf_id = db.Column(db.ForeignKey('gene.id'), index=True) # it's a gene ID
    gene_id = db.Column(db.ForeignKey('gene.id'), index=True) # from the JSON file

    gene = db.relationship('Gene', primaryjoin='TfTarget.gene_id == Gene.id', backref='gene_tf_targets')
    tf = db.relationship('Gene', primaryjoin='TfTarget.tf_id == Gene.id', backref='gene_tf_targets_0')

# used
class Locus(db.Model):
    __tablename__ = 'locus'

    id = db.Column(db.Integer, primary_key=True)
    locus_name = db.Column(db.String(10))
    mutation_type = db.Column(db.String(10))

# used
class Eigengene(db.Model):
    __tablename__ = 'eigengene'

    id = db.Column(db.Integer, primary_key=True)
    patient_id = db.Column(db.ForeignKey('patient.id'), index=True)
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    value = db.Column(db.Float)

    bicluster = db.relationship('Bicluster', primaryjoin='Eigengene.bicluster_id == Bicluster.id', backref='eigengene')
    patient = db.relationship('Patient', primaryjoin='Eigengene.patient_id == Patient.id', backref='eigengene')

class BiclusterPhenotypeSignificance(db.Model):
    __tablename__ = 'bicluster_phenotype_significance'

    id = db.Column(db.Integer, primary_key=True)
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    phenotype_id = db.Column(db.ForeignKey('phenotype.id'), index=True)
    p_value = db.Column(db.Float)
    r_value = db.Column(db.Float)

    bicluster = db.relationship('Bicluster', primaryjoin='BiclusterPhenotypeSignificance.bicluster_id == Bicluster.id', backref='bicluster_phenotype_significance')
    phenotype = db.relationship('Phenotype', primaryjoin='BiclusterPhenotypeSignificance.phenotype_id == Phenotype.id', backref='bicluster_phenotype_significance')

class CellLine(db.Model):
    __tablename__ = 'cell_line'

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(10))

class GeneJACKSResult(db.Model):
    __tablename__ = 'gene_jacks_result'

    id = db.Column(db.Integer, primary_key=True)
    gene_id = db.Column(db.ForeignKey('gene.id'), index=True)
    cell_line_id = db.Column(db.ForeignKey('cell_line.id'), index=True)
    score = db.Column(db.Float)
    std = db.Column(db.Float)
    p_value = db.Column(db.Float)

    gene = db.relationship('Gene', primaryjoin='GeneJACKSResult.gene_id == Gene.id', backref='gene_jacks_result')
    cell_line = db.relationship('CellLine', primaryjoin='GeneJACKSResult.cell_line_id == CellLine.id', backref='gene_jacks_result')

class GRNA(db.Model):
    __tablename__ = 'grna'

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(30))
    gene_id = db.Column(db.ForeignKey('gene.id'), index=True)

    gene = db.relationship('Gene', primaryjoin='GRNA.gene_id == Gene.id', backref='grna')

class GRNAJACKSResult(db.Model):
    __tablename__ = 'grna_jacks_result'

    id = db.Column(db.Integer, primary_key=True)
    grna_id = db.Column(db.ForeignKey('grna.id'), index=True)
    cell_line_id = db.Column(db.ForeignKey('cell_line.id'), index=True)
    mean = db.Column(db.Float)
    std = db.Column(db.Float)
    p_value = db.Column(db.Float)

    grna = db.relationship('GRNA', primaryjoin='GRNAJACKSResult.grna_id == GRNA.id', backref='grna_jacks_result')
    cell_line = db.relationship('CellLine', primaryjoin='GRNAJACKSResult.cell_line_id == CellLine.id', backref='grna_jacks_result')

class AchillesCommonEssential(db.Model):
    __tablename__ = 'achilles_common_essential'
    
    id = db.Column(db.Integer, primary_key=True)
    gene_id = db.Column(db.ForeignKey('gene.id'), index=True)

    gene = db.relationship('Gene', primaryjoin='AchillesCommonEssential.gene_id == Gene.id', backref='achilles_common_essential')

class GeneAchillesResult(db.Model):
    __tablename__ = 'achilles_results'

    id = db.Column(db.Integer, primary_key=True)
    gene_id = db.Column(db.ForeignKey('gene.id'), index=True)
    cell_line_id = db.Column(db.ForeignKey('cell_line.id'), index=True)
    score = db.Column(db.Float)
    p_value = db.Column(db.Float)

    gene = db.relationship('Gene', primaryjoin='GeneAchillesResult.gene_id == Gene.id', backref='achilles_results')
    cell_line = db.relationship('CellLine', primaryjoin='GeneAchillesResult.cell_line_id == CellLine.id', backref='achilles_results')

class GRNASequence(db.Model):
    __tablename__ = 'grna_sequence'

    id = db.Column(db.Integer, primary_key=True)
    sequence = db.Column(db.String(50))
    grna_id = db.Column(db.ForeignKey('grna.id'), index=True)

    grna = db.relationship('GRNA', primaryjoin='GRNASequence.grna_id == GRNA.id', backref='grna_sequence')