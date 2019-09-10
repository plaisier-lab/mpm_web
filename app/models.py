# coding: utf-8
from sqlalchemy import BINARY, Column, Float, ForeignKey, Integer, String, Text
from sqlalchemy.orm import relationship
from flask import Flask
from flask_sqlalchemy import SQLAlchemy

###########################
## Open database session ##
###########################
app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql+mysqlconnector://root:root@mpm_web_db_1:3306/mpm_web'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = 'FALSE'

db = SQLAlchemy(app)

class BicGene(db.Model):
    __tablename__ = 'bic_gene'

    id = db.Column(db.Integer, primary_key=True)
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    gene_id = db.Column(db.ForeignKey('gene.id'), index=True)

    bicluster = db.relationship('Bicluster', primaryjoin='BicGene.bicluster_id == Bicluster.id', backref='bic_genes')
    gene = db.relationship('Gene', primaryjoin='BicGene.gene_id == Gene.id', backref='bic_genes')


class BicGo(db.Model):
    __tablename__ = 'bic_go'

    id = db.Column(db.Integer, primary_key=True)
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    go_bp_id = db.Column(db.ForeignKey('go_bp.id'), index=True)

    bicluster = db.relationship('Bicluster', primaryjoin='BicGo.bicluster_id == Bicluster.id', backref='bic_goes')
    go_bp = db.relationship('GoBp', primaryjoin='BicGo.go_bp_id == GoBp.id', backref='bic_goes')


class BicHal(db.Model):
    __tablename__ = 'bic_hal'

    id = db.Column(db.Integer, primary_key=True)
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    hallmark_id = db.Column(db.ForeignKey('hallmark.id'), index=True)

    bicluster = db.relationship('Bicluster', primaryjoin='BicHal.bicluster_id == Bicluster.id', backref='bic_hals')
    hallmark = db.relationship('Hallmark', primaryjoin='BicHal.hallmark_id == Hallmark.id', backref='bic_hals')


class BicPat(db.Model):
    __tablename__ = 'bic_pat'

    id = db.Column(db.Integer, primary_key=True)
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    patient_id = db.Column(db.ForeignKey('patient.id'), index=True)

    bicluster = db.relationship('Bicluster', primaryjoin='BicPat.bicluster_id == Bicluster.id', backref='bic_pats')
    patient = db.relationship('Patient', primaryjoin='BicPat.patient_id == Patient.id', backref='bic_pats')


class Bicluster(db.Model):
    __tablename__ = 'bicluster'

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(20), index=True)
    var_exp_fpc = db.Column(db.Float)
    var_exp_fpc_p_value = db.Column(db.Float)
    survival = db.Column(db.Float)
    survival_p_value = db.Column(db.Float)


class CausalFlow(db.Model):
    __tablename__ = 'causal_flow'

    id = db.Column(db.Integer, primary_key=True)
    somatic_mutation_id = db.Column(db.ForeignKey('somatic_mutation.id'), index=True)
    regulator_id = db.Column(db.Integer)
    regulator_type = db.Column(db.String(10))
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    leo_nb_atob = db.Column(db.Float)
    mlogp_m_atob = db.Column(db.Float)

    bicluster = db.relationship('Bicluster', primaryjoin='CausalFlow.bicluster_id == Bicluster.id', backref='causal_flows')
    somatic_mutation = db.relationship('SomaticMutation', primaryjoin='CausalFlow.somatic_mutation_id == SomaticMutation.id', backref='causal_flows')


class ExpDataset(db.Model):
    __tablename__ = 'exp_dataset'

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(20))
    type = db.Column(db.String(10))


class Gene(db.Model):
    __tablename__ = 'gene'

    id = db.Column(db.Integer, primary_key=True)
    symbol = db.Column(db.String(100), index=True)
    entrez = db.Column(db.Integer)


class GenePrior(db.Model):
    __tablename__ = 'gene_prior'

    id = db.Column(db.Integer, primary_key=True)
    gene_id = db.Column(db.ForeignKey('gene.id'), index=True)
    source = db.Column(db.String(11), index=True)
    description = db.Column(db.String(100), index=True)
    pmid = db.Column(db.Integer)


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


class GoBp(db.Model):
    __tablename__ = 'go_bp'

    id = db.Column(db.Integer, primary_key=True)
    go_id = db.Column(db.String(15), index=True)
    name = db.Column(db.Text)
    description = db.Column(db.Text)


class GoGene(db.Model):
    __tablename__ = 'go_gene'

    id = db.Column(db.Integer, primary_key=True)
    go_bp_id = db.Column(db.ForeignKey('go_bp.id'), index=True)
    gene_id = db.Column(db.ForeignKey('gene.id'), index=True)

    gene = db.relationship('Gene', primaryjoin='GoGene.gene_id == Gene.id', backref='go_genes')
    go_bp = db.relationship('GoBp', primaryjoin='GoGene.go_bp_id == GoBp.id', backref='go_genes')


class Hallmark(db.Model):
    __tablename__ = 'hallmark'

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100), index=True)


class Mirna(db.Model):
    __tablename__ = 'mirna'

    id = db.Column(db.Integer, primary_key=True)
    mature_seq_id = db.Column(db.String(20))
    name = db.Column(db.String(50), index=True)


class MirnaPrior(db.Model):
    __tablename__ = 'mirna_prior'

    id = db.Column(db.Integer, primary_key=True)
    mirna_id = db.Column(db.ForeignKey('mirna.id'), index=True)
    source = db.Column(db.String(11), index=True)
    description = db.Column(db.String(100))
    pmid = db.Column(db.Integer)


class MirnaExpression(db.Model):
    __tablename__ = 'mirna_expression'

    id = db.Column(db.Integer, primary_key=True)
    exp_dataset_id = db.Column(db.ForeignKey('exp_dataset.id'), index=True)
    patient_id = db.Column(db.ForeignKey('patient.id'), index=True)
    mirna_id = db.Column(db.ForeignKey('mirna.id'), index=True)
    value = db.Column(db.Float)

    exp_dataset = db.relationship('ExpDataset', primaryjoin='MirnaExpression.exp_dataset_id == ExpDataset.id', backref='mirna_expressions')
    mirna = db.relationship('Mirna', primaryjoin='MirnaExpression.mirna_id == Mirna.id', backref='mirna_expressions')
    patient = db.relationship('Patient', primaryjoin='MirnaExpression.patient_id == Patient.id', backref='mirna_expressions')


class MirnaRegulator(db.Model):
    __tablename__ = 'mirna_regulator'

    id = db.Column(db.Integer, primary_key=True)
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    mirna_id = db.Column(db.ForeignKey('mirna.id'), index=True)

    bicluster = db.relationship('Bicluster', primaryjoin='MirnaRegulator.bicluster_id == Bicluster.id', backref='mirna_regulators')
    mirna = db.relationship('Mirna', primaryjoin='MirnaRegulator.mirna_id == Mirna.id', backref='mirna_regulators')


class MirnaTarget(db.Model):
    __tablename__ = 'mirna_targets'

    id = db.Column(db.Integer, primary_key=True)
    mirna_id = db.Column(db.ForeignKey('mirna.id'), index=True)
    gene_id = db.Column(db.ForeignKey('gene.id'), index=True)
    source = db.Column(db.String(11))

    gene = db.relationship('Gene', primaryjoin='MirnaTarget.gene_id == Gene.id', backref='mirna_targets')
    mirna = db.relationship('Mirna', primaryjoin='MirnaTarget.mirna_id == Mirna.id', backref='mirna_targets')


class Patient(db.Model):
    __tablename__ = 'patient'

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(50), index=True)


class PhenoDatum(db.Model):
    __tablename__ = 'pheno_data'

    id = db.Column(db.Integer, primary_key=True)
    phenotype_id = db.Column(db.ForeignKey('phenotype.id'), index=True)
    patient_id = db.Column(db.ForeignKey('patient.id'), index=True)
    name = db.Column(db.String(50), nullable=False)
    long_name = db.Column(db.String(50))

    patient = db.relationship('Patient', primaryjoin='PhenoDatum.patient_id == Patient.id', backref='pheno_data')
    phenotype = db.relationship('Phenotype', primaryjoin='PhenoDatum.phenotype_id == Phenotype.id', backref='pheno_data')


class Phenotype(db.Model):
    __tablename__ = 'phenotype'

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(50), nullable=False)
    long_name = db.Column(db.String(50))


class Replication(db.Model):
    __tablename__ = 'replication'

    id = db.Column(db.Integer, primary_key=True)
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    study = db.Column(db.String(50))
    var_exp_fpc = db.Column(db.Float)
    var_exp_fpc_p_value = db.Column(db.Float)
    survival = db.Column(db.Float)
    survival_p_value = db.Column(db.Float)

    bicluster = db.relationship('Bicluster', primaryjoin='Replication.bicluster_id == Bicluster.id', backref='replications')


class SomaticMutation(db.Model):
    __tablename__ = 'somatic_mutation'

    id = db.Column(db.Integer, primary_key=True)
    ext_id = db.Column(db.Integer, index=True)
    mutation_type = db.Column(db.String(10))


class TfRegulator(db.Model):
    __tablename__ = 'tf_regulator'

    id = db.Column(db.Integer, primary_key=True)
    bicluster_id = db.Column(db.ForeignKey('bicluster.id'), index=True)
    gene_id = db.Column(db.ForeignKey('gene.id'), index=True)
    action = db.Column(db.String(15))

    bicluster = db.relationship('Bicluster', primaryjoin='TfRegulator.bicluster_id == Bicluster.id', backref='tf_regulators')
    gene = db.relationship('Gene', primaryjoin='TfRegulator.gene_id == Gene.id', backref='tf_regulators')


class TfTarget(db.Model):
    __tablename__ = 'tf_targets'

    id = db.Column(db.Integer, primary_key=True)
    tf_id = db.Column(db.ForeignKey('gene.id'), index=True)
    gene_id = db.Column(db.ForeignKey('gene.id'), index=True)

    gene = db.relationship('Gene', primaryjoin='TfTarget.gene_id == Gene.id', backref='gene_tf_targets')
    tf = db.relationship('Gene', primaryjoin='TfTarget.tf_id == Gene.id', backref='gene_tf_targets_0')
