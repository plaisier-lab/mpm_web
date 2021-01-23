import re

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