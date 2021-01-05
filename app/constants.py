NUM_PARTS = 5

HALLMARK_TO_ICON = {
	'Evading apoptosis': 'cellDeath.gif',
	'Evading immune detection': 'avoidImmuneDestruction.gif',
	'Genome instability and mutation': 'genomicInstability.gif',
	'Insensitivity to antigrowth signals': 'evadeGrowthSuppressors.gif',
	'Limitless replicative potential': 'immortality.gif',
  'Reprogramming energy metabolism': 'cellularEnergetics.gif',
	'Self sufficiency in growth signals': 'sustainedProliferativeSignalling.gif',
	'Sustained angiogenesis': 'angiogenesis.gif',
	'Tissue invasion and metastasis': 'invasion.gif',
	'Tumor promoting inflammation': 'promotingInflammation.gif'
}

HALLMARKS = [
	'Self sufficiency in growth signals',
	'Reprogramming energy metabolism',
	'Evading apoptosis',
	'Genome instability and mutation',
	'Sustained angiogenesis',
	'Tissue invasion and metastasis',
	'Tumor promoting inflammation',
	'Limitless replicative potential',
	'Evading immune detection',
	'Insensitivity to antigrowth signals',
]

ENRICHMENT_PHENOTYPES = [
	"Sarcomatoid", "Epithelioid", "Desmoplastic", "Biphasic",
]

SEX_PHENOTYPES = ['M', 'F']

PHENOTYPES_DICT = {
	"sex": SEX_PHENOTYPES,
	"histology_WHO": ENRICHMENT_PHENOTYPES,
}

GRAPH_COLOR_MAP = {
	'Sarcomatoid': '#2c6dd4',
	'Epithelioid': '#8a171a',
	'Desmoplastic': '#ed2024',
	'Biphasic': '#faa41a',
	'M': '#6abd45',
	'F': '#2c6dd4',
	'NA': 'grey',
}

GRAPH_COLOR_GRADIENTS = {
	'age_at_surgery': ('#333333', '#2cd46d'),
	'preop_treatment': ('#333333', '#2cd46d'),
	'B.cells.naive': ('#333333', '#2cd46d'),
	'B.cells.memory': ('#333333', '#2cd46d'),
	'Plasma.cells': ('#333333', '#2cd46d'),
	'T.cells.CD8': ('#333333', '#2cd46d'),
	'T.cells.CD4.naive': ('#333333', '#2cd46d'),
	'T.cells.CD4.memory.resting': ('#333333', '#2cd46d'),
	'T.cells.CD4.memory.activated': ('#333333', '#2cd46d'),
	'T.cells.follicular.helper': ('#333333', '#2cd46d'),
	'T.cells.regulatory..Tregs.': ('#333333', '#2cd46d'),
	'T.cells.gamma.delta': ('#333333', '#2cd46d'),
	'NK.cells.resting': ('#333333', '#2cd46d'),
	'NK.cells.activated': ('#333333', '#2cd46d'),
	'Monocytes': ('#333333', '#2cd46d'),
	'Macrophages.M0': ('#333333', '#2cd46d'),
	'Macrophages.M1': ('#333333', '#2cd46d'),
	'Macrophages.M2': ('#333333', '#2cd46d'),
	'Dendritic.cells.resting': ('#333333', '#2cd46d'),
	'Dendritic.cells.activated': ('#333333', '#2cd46d'),
	'Mast.cells.resting': ('#333333', '#2cd46d'),
	'Mast.cells.activated': ('#333333', '#2cd46d'),
	'Eosinophils': ('#333333', '#2cd46d'),
	'Neutrophils': ('#333333', '#2cd46d'),
}

SELECTABLE_PHENOTYPES = [
	("WHO Histology", "histology_WHO"),
	("Patient Sex", "sex"),
	("Age At Surgery", "age_at_surgery"),
	("Preop Treatment", "preop_treatment"),
	("B-Cells Naive", "B.cells.naive"),
	("B-Cells Memory", "B.cells.memory"),
	("Plasma Cells", "Plasma.cells"),
	("T-Cells CD8", "T.cells.CD8"),
	("T-Cells CD4 Naive", "T.cells.CD4.naive"),
	("T-Cells CD4 Memory Resting", "T.cells.CD4.memory.resting"),
	("T-Cells CD4 Memory Activated", "T.cells.CD4.memory.activated"),
	("T-Cells Follicular Helper", "T.cells.follicular.helper"),
	("T-Cells Regulatory TRegs", "T.cells.regulatory..Tregs."),
	("T-Cells Gamma Delta", "T.cells.gamma.delta"),
	("NK-Cells Resting", "NK.cells.resting"),
	("NK-Cells Activated", "NK.cells.activated"),
	("Monocytes", "Monocytes"),
	("Macrophages M0", "Macrophages.M0"),
	("Macrophages M1", "Macrophages.M1"),
	("Macrophages M2", "Macrophages.M2"),
	("Dendritic Cells Resting", "Dendritic.cells.resting"),
	("Dendritic Cells Activated", "Dendritic.cells.activated"),
	("Mast Cells Resting", "Mast.cells.resting"),
	("Mast Cells Activated", "Mast.cells.activated"),
	("Eosinophils", "Eosinophils"),
	("Neutrophils", "Neutrophils"),
]

# blacklist for phenotypes we don't want to check p-values for
SELECTABLE_PHENOTYPES_BLACKLIST = ["histology_WHO", "sex", "age_at_surgery", "preop_treatment"]

USE_PHENOTYPE_SCATTERPLOT = True

PHENOTYPE_INDEX_TO_UINAME = {item[1]: item[0] for item in SELECTABLE_PHENOTYPES}