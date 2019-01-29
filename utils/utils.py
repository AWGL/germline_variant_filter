import yaml
import csv
import pandas as pd
import requests
import datetime

def parse_config(yaml_file):
	"""
	Parse the yaml config file

	"""

	config_dict = yaml.load(open(yaml_file))

	return config_dict


def parse_ped_file(ped_file):

	"""
	Parse a ped file into a dictionary with each sample as a key.

	"""

	ped_dict = {}
	with open(ped_file) as csvfile:
		spamreader = csv.reader(csvfile, delimiter='\t')
		for row in spamreader:
			ped_dict[row[1]] = {}
			ped_dict[row[1]]['familyID'] = row[0]
			ped_dict[row[1]]['paternalID'] = row[2]
			ped_dict[row[1]]['maternalID'] = row[3]
			ped_dict[row[1]]['sex'] = row[4]
			ped_dict[row[1]]['affected'] = row[5]


	return ped_dict


def parse_panel_app_dump(panel_app_dump):
	"""
	Parse the local panelapp dump
	"""
	panel_app_dict = {}
	with open(panel_app_dump) as csvfile:

		spamreader = csv.reader(csvfile, delimiter='\t')
	
		for row in spamreader:

			panel_app_dict[row[0]] = {}
			panel_app_dict[row[0]]['date'] = row[1]
			panel_app_dict[row[0]]['disease'] = row[2]
			panel_app_dict[row[0]]['inheritance'] = row[3]
			
	return panel_app_dict


def is_proband_in_trio(sample, ped_dict):
	"""
	Is the sample a proband in a trio?

	"""
	if ped_dict[sample]['paternalID'] == '0':

		return False

	if ped_dict[sample]['maternalID'] == '0': 

		return False

	return True


def select_variants_for_sample(df, sample):
	"""
	Checks whether the variant has not been filtered for that variant.
	
	checks if the sample is not reference for this variant.
	
	if both of these are True then we return True
	
	else we return False
	"""
	
	passes_filter = False
	is_variant = False
	
	genotypes = df[sample + '.GT']
	
	if '|' in genotypes:
		
		genotypes = genotypes.split('|')
	
	elif '/' in genotypes:
		
		genotypes = genotypes.split('/')
		
	else:
		
		raise Exception('There is a genotype without either | or / in.')
		
	if df['ALT'] in genotypes:
		
		is_variant = True
		
	if pd.isna(df[sample + '.FT']) == True or df[sample + '.FT'] == 'PASS':
		
		passes_filter = True
		
	if is_variant == True and passes_filter == True:
		  
		  return True
		  
	else:
		  
		  return False


def fix_column_names(columns):
	"""
	Change column names to valid python variable names by replacing '.' with '_'

	We need this for when we loop through the dataframe and put each VEP consequence on a seperate row.

	
	"""
	
	fixed_columns = []
	
	for column_name in columns:
		
		fixed_column_name = column_name.replace('.', '_').replace('-', '_')
		
		if fixed_column_name not in ['CHROM', 'POS', 'REF', 'ALT','ID','CSQ', 'QUAL', 'AC', 'TYPE','FILTER']:
			
			fixed_column_name = 'sample_' + fixed_column_name
		
		fixed_columns.append(fixed_column_name)
	
	return fixed_columns


def parse_csq(csq_string, csq_desc):
	"""
	Take the raw VEP string and parse it into a list of dictionaries containing the information.
	
	"""
	
	transcript_list = []
	
	# if we don't have a CSQ string then set everything to None
	if pd.isna(csq_string) == True:
		
		transcript_dict = {}
		
		for key in csq_desc:
			
			transcript_dict[key] = None
			
		transcript_list.append(transcript_dict)
		
		return transcript_list
	
	# Split out the consequence blocks
	transcripts = csq_string.split(',')
	
	for transcript in transcripts:
		
		transcript_dict = {}
		
		vep_data = transcript.split('|')
		
		for key, value in zip(csq_desc, vep_data):
			
			transcript_dict[key] = value
			
		transcript_list.append(transcript_dict)
		
	return transcript_list
		
		
def split_vep_transcripts(df, csq_desc, vep_fields, column_names):
	"""
	Takes a df and splits out the CSQ field so that each line in the dataframe is a \
	seperate transcript.
	
	df = dataframe
	csq_desc = the csq string from the vcf header
	vep_fields = the vep fields to extract
	column_names = the new columns for the new df
	
	"""
	
	new_df_list = []

	for variant in df.itertuples():
		
		csq = variant.CSQ
		
		# Parse csq string into a list of dictionaries
		transcript_list = parse_csq(csq, csq_desc)
		
		# Add the standard vcf columns

		vcf_row = []
		for data in variant:
			
			vcf_row.append(data)
		
		# Combine the standard columns and the genotype columns                
		for transcript in transcript_list:
			
			vep_row = []
			
			for vep_field in vep_fields:
				
				data = transcript[vep_field]
				
				if data == '.':
					
					data = None
				
				vep_row = vep_row + [data]
				
			final_row = vcf_row + vep_row
			
			# first item in list is the index so exclude
			new_df_list.append(final_row[1:])
			
	for vep_field in vep_fields:
		
		column_names = column_names + [vep_field]
		
	new_df = pd.DataFrame(new_df_list, columns = column_names)
	
	return new_df

def get_variant_key(df):
	"""
	Make a key for the variant.

	"""

	return f"{df['CHROM']}:{df['POS']}{df['REF']}>{df['ALT']}"

def get_worst_consequence(consequences, consequence_severity):
	"""
	Of all the transcripts that a variant falls in what is the worst consequence.

	"""
	
	worst_index = 9999
	
	for consequence in consequences:
		
		if consequence == None:
			
			break
		
		split_consequences = consequence.split('&')
		
		for c in split_consequences:
		
			index = consequence_severity.index(c)
		
			if index < worst_index:
			
				worst_index = index
	
	if worst_index == 9999:
		
		return None
	
	else:
			
		return consequence_severity[worst_index]
	
	
def consequence_filter(df, to_keep_consequences):
	"""
	Whether to filter the variant on consequence?
	
	"""
	
	if df['worst_consequence'] in to_keep_consequences:
		
		return False
	
	else:
		
		return True
	
	
def has_important_clinsig(df, clin_sig_words):
	"""
	Does the variant have a CLINSIG field with a keyword in?

	For example pathogenic.

	"""
	
	clin_sig = df['CLIN_SIG']
	
	if clin_sig == None:
		
		return None
	
	for word in clin_sig_words:
		
		if word.lower() in clin_sig.lower():
			
			return True
		
	return False


def fix_splice_ai(df, column):
	"""
	Split splice ai if there are more than one entry
	
	return the one related to out gene of interest
	
	"""
	
	spliceai_symbols = []
	
	if df[column] == '' or df[column] == None:
		
		return None
	
	if '&' in df[column]:
		
		spliceai_symbols = df['SpliceAI_SYMBOL'].split('&')
		
		for i, symbol in enumerate(spliceai_symbols):
			
			if symbol == df['SYMBOL']:
				
				return df[column].split('&')[i]
		
			
	return max(df[column].split('&'))

def has_affect_on_splicing(df, cutoff):
	"""
	Does a variant have an affect on splicing


	"""
	
	if df['fixed_SpliceAI_DS_AG'] >= cutoff or pd.isna(df['fixed_SpliceAI_DS_AG']):
		
		return True
	
	elif df['fixed_SpliceAI_DS_AL'] >= cutoff or pd.isna(df['fixed_SpliceAI_DS_AL']):
		
		return True
	
	elif df['fixed_SpliceAI_DS_DG'] >= cutoff or pd.isna(df['fixed_SpliceAI_DS_DG']):
		
		return True
	
	elif df['fixed_SpliceAI_DS_DL'] >= cutoff or pd.isna(df['fixed_SpliceAI_DS_DL']):
		
		return True
	
	else:
		
		return False


def any_has_splicing_affect(splicing):
	"""
	Do any of the transcripts have an affect on splicing?

	"""
	
	if any(splicing) == True:
		
		return True
	
	else:
		
		return False

def fix_ccrs(df):
	"""
	Fixes the ccrs annoatation so that when we have multiple annotations for a \
	variant we take the largest one.

	"""
	
	ccrs =  df['ccrs']
	
	if ccrs == '.' or ccrs == None:
		
		return None
	
	if '&' in ccrs:
		
		ccrs =ccrs.split('&')
	
		ccrs = [float(x) for x in ccrs]
	
		return max(ccrs)
	
	else:
		
		try:
		
			return float(ccrs)
		
		except:
			
			return None

def annotate_with_gnomad_scores(df, gnomad_scores_dict, value):
	"""
	Add the gnomad constraint score for that transcript.

	"""
	
	feature = df['Feature']
	
	if feature == None:
		
		return None
	
	feature = feature.split('.')[0]
	
	if feature in gnomad_scores_dict[value]:
		
		return gnomad_scores_dict[value][feature]
	
	else:
		
		return None

def get_panel_app_info(gene):
	"""
	Query panel app and get the data for a specific gene.

	"""
	
	url = f'https://panelapp.genomicsengland.co.uk/WebServices/search_genes/{gene}/?format=json&LevelOfConfidence=HighEvidence'
	
	try:
		panel_response = requests.get(url)

		panel_data = panel_response.json()
		
		if panel_data['meta']['numOfResults'] ==0:
		
			return None, None
	
		diseases = []
		modes_of_inheritance = []
	
		for gene_info in panel_data['results']:
				
			if float(gene_info['version']) > 1.0:
		   
				diseases.append(gene_info.get('SpecificDiseaseName', 'None'))
				
				modes_of_inheritance.append(gene_info.get('ModeOfInheritance', 'None'))
		
		return '|'.join(list(set(diseases))), '|'.join(set(modes_of_inheritance))
	
	except:
		
		return None, None

def apply_panel_app_data(df, panel_app_dict, panel_app_dump_max_time):
	"""
	Get the panel app data for each row in the dataframe.

	"""
	
	symbol = df['SYMBOL']
	
	if symbol in panel_app_dict:

		today = datetime.datetime.now()

		last_analyzed = panel_app_dict[symbol]['date']

		last_analyzed = last_analyzed.split(' ')[0].split('-')

		last_analyzed = datetime.datetime(int(last_analyzed[0]), int(last_analyzed[1]), int(last_analyzed[2]))

		difference = today - last_analyzed

		if difference.days < panel_app_dump_max_time:
		
			return panel_app_dict[symbol]['disease']
	
	diseases, inheritance = get_panel_app_info(symbol)
		
	date = datetime.datetime.now()
		
	panel_app_dict[symbol] = {'disease': diseases, 'inheritance': inheritance, 'date': str(date)}

	print (symbol)
		
	return panel_app_dict[symbol]['disease']


def parse_hpo_file(hpo_file):
	"""
	Parse the file containing the HPO terms and which gene they map to.

	"""
	
	hpo_dict = {}
	
	with open(hpo_file, 'r') as csvfile:
		
		spamreader = csv.reader(csvfile, delimiter='\t' )
		next(spamreader)
		i =1
		for row in spamreader:
			i = 1 +1
			if i >10:
				break
			
			if row[0] not in hpo_dict:
				
				hpo_dict[row[0]] = {}
				hpo_dict[row[0]][row[3]] = row[3]
				
			else:
				
				hpo_dict[row[0]][row[3]] = row[3]
				
	return hpo_dict


