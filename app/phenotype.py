"""
helper functions for displaying phenotype data to the user
"""

from constants import GRAPH_COLOR_GRADIENTS
from colour import Color
import numpy as np

def get_phenotype_min_max(cursor, phenotype_name, bueno_deep_filter=False):
	if bueno_deep_filter == False:
		cursor.execute("""
			SELECT phenotype_value FROM pheno_data pd
				JOIN patient p ON pd.patient_id=p.id
				JOIN phenotype pt ON pd.phenotype_id=pt.id
				WHERE pt.name=%s;
		""", [phenotype_name])
	else:
		cursor.execute("""
			SELECT phenotype_value FROM pheno_data pd
				JOIN patient p ON pd.patient_id=p.id
				JOIN phenotype pt ON pd.phenotype_id=pt.id
				WHERE pt.name=%s
				AND p.id IN (SELECT DISTINCT patient_id FROM bueno_deep_filter);
		""", [phenotype_name])
	data = cursor.fetchall()
	phenotype_data = [item[0] for item in data]

	if phenotype_data[0] == None:
		return None
	else:
		return (np.min(phenotype_data), np.max(phenotype_data))

def get_phenotype_color(phenotype_name, value, min_max):
	min, max = min_max
	color1 = Color(GRAPH_COLOR_GRADIENTS[phenotype_name][0])
	color2 = Color(GRAPH_COLOR_GRADIENTS[phenotype_name][1])
	percent = (value - min) / (max - min)
	hex = Color(rgb=(color1.red * (1 - percent) + color2.red * percent, color1.green * (1 - percent) + color2.green * percent, color1.blue * (1 - percent) + color2.blue * percent)).hex

	if len(hex) == 4:
		return hex + hex[1] + hex[1] + hex[1]
	else:
		return hex