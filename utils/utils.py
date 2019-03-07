import yaml
import csv
import pandas as pd
import requests
import datetime
import sqlite3

def parse_config(yaml_file):
	"""
	Parse the yaml config file containing the preferences.

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
	Parse the local panelapp dump. This stores the data retrieved from the PanelApp API locally so \
	we don't have to requery if the data was retrived within a certian time limit.
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


def select_variants_for_sample(df, sample, min_dp, min_gq, gt_depth_tag):
	"""
	Checks whether the variant has not been filtered for that variant.
	
	Checks if the sample is not reference for this variant.
	
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
		
	if df[sample + '.GQ'] >= min_gq and df[sample + '.' + gt_depth_tag] >= min_dp:
		
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
		
		if fixed_column_name not in ['CHROM', 'POS', 'REF', 'ALT','ID','CSQ', 'QUAL', 'AC', 'TYPE','FILTER', 'OLD_CLUMPED']:
			
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

		if len(csq_desc) != len(vep_data):

			raise Exception('VEP CSQ String Mismatch! Check the input CSQ description.')
		
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

	Takes a list consequence_severity which has all consequences in order of severity

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
	
	if df['WorstConsequence'] in to_keep_consequences:
		
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
	Split SpliceAi result if there are more than one entry.
	
	Return the one related to our gene of interest.

	Otherwise return the highest of the two.
	
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
	Does a variant have an affect on splicing?

	"""
	
	if df['SpliceAI_DS_AG'] >= cutoff or pd.isna(df['SpliceAI_DS_AG']):
		
		return True
	
	elif df['SpliceAI_DS_AL'] >= cutoff or pd.isna(df['SpliceAI_DS_AL']):
		
		return True
	
	elif df['SpliceAI_DS_DG'] >= cutoff or pd.isna(df['SpliceAI_DS_DG']):
		
		return True
	
	elif df['SpliceAI_DS_DL'] >= cutoff or pd.isna(df['SpliceAI_DS_DL']):
		
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

		try:
	
			return max(ccrs)

		except:

			return None
	
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
		panel_response = requests.get(url, timeout=5)

		panel_data = panel_response.json()
		
		if panel_data['meta']['numOfResults'] == 0:
		
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

def apply_panel_app_data_disease(df, panel_app_dict, panel_app_dump_max_time):
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
		
	return panel_app_dict[symbol]['disease']


def apply_panel_app_data_inheritance(df, panel_app_dict, panel_app_dump_max_time):
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
		
			return panel_app_dict[symbol]['inheritance']
	
	diseases, inheritance = get_panel_app_info(symbol)
		
	date = datetime.datetime.now()
		
	panel_app_dict[symbol] = {'disease': diseases, 'inheritance': inheritance, 'date': str(date)}
		
	return panel_app_dict[symbol]['inheritance']


def parse_hpo_file(hpo_file):
	"""
	Parse the file containing the HPO terms and which gene they map to.

	"""
	
	hpo_dict = {}
	
	with open(hpo_file, 'r') as csvfile:
		
		spamreader = csv.reader(csvfile, delimiter='\t' )
		next(spamreader)

		for row in spamreader:

			if row[0] not in hpo_dict:
				
				hpo_dict[row[0]] = {}
				hpo_dict[row[0]][row[3]] = row[3]
				
			else:
				
				hpo_dict[row[0]][row[3]] = row[3]
				
	return hpo_dict


def get_genotype(df, sample):

	"""
	Return whether the proband gentype is HET or HOM

	"""

	sample_gt = df[f'sample_{sample}_GT']
	alt = df['ALT']

	if '|' in sample_gt:
			
		sample_gt = sample_gt.split('|')
	
	elif '/' in sample_gt:
			
		sample_gt = sample_gt.split('/')
			
	else:
			
		raise Exception('There is a sample genotype without either | or / in.')

	if sample_gt.count(alt) == 1:

		return 'HET'

	elif sample_gt.count(alt) == 2:

		return 'HOM'

	else:

		return 'UNKNOWN'


def fix_exon(df):

	"""
	Return a fixed version of the exon for display in excel

	"""

	exon = df['EXON']

	if exon == None or exon == '':

		return ''

	else:

		return exon.replace('/', '|')



def fix_intron(df):

	"""
	Return a fixed version of the intron for display in excel

	"""

	intron = df['INTRON']

	if intron == None or intron == '':

		return ''

	else:

		return intron.replace('/', '|')



def get_hgvsc(df):
	"""
	Return formatted hgvsc.

	"""

	hgvsc = df['HGVSc']

	if hgvsc == None or hgvsc == '':

		return ''

	else:

		return hgvsc.split(':')[1]


def get_hgvsp(df):
	"""
	Return formatted hgvsp

	"""

	hgvsp = df['HGVSp']

	if hgvsp == None or hgvsp == '':

		return ''

	else:

		return hgvsp.split(':')[1]


def annotate_hpo(df, patient_hpos, hpo_dict):
	"""
	Add the count of the matching HPO terms.

	"""

	gene_id = df['Gene']

	counter = 0

	if gene_id in hpo_dict:
		
		for patient_hpo in set(patient_hpos):

			if patient_hpo in hpo_dict[gene_id]:

				counter = counter + 1

	return counter


def check_picks(df, pick_dict):
	"""
	Check every sample has a pick flag.

	If it doesn't then set all to 1.

	"""

	variant_id = df['VariantId']

	if pick_dict['PICK'][variant_id] != '1':

		return '1'

	else:

		return df['PICK']


def are_arguments_valid(args, config_dict):

	"""
	Check whether the arguments and config make sense.

	"""

	if args.local_panel_app_dump != None and args.panelapp == False:

		print ('Cannot select to use PanelApp dump and not select to add PanelApp Data.')
		return False

	if args.spliceai == True and 'SpliceAI' not in args.csq[0]:

		print ('Input file does not contain the required annotations for the spliceai option.')
		return False

	if args.smart_synonymous == True and (('SpliceAI' not in args.csq[0] and 'CLIN_SIG' not in args.csq) or args.spliceai ==False) :

		print ('Input file does not contain the required annotations for the smart-synonymous option or you have not selected to parse the spliceAI columns.')
		return False

	if args.add_ccrs == True and 'ccrs' not in args.csq[0]:

		print ('Input file does not contain the required annotations for the add-ccrs option.')
		return False

	if args.add_ccrs == True and ('CCR_percentile' not in  config_dict['final_fields_trio'] or 'CCR_percentile' not in config_dict['final_fields_single']):

		print ('CCR_percentile needs to be in the final fields variable in the config file if you select the add-ccrs option.')
		return False	

	if args.gnomad_constraint_scores == True and 'gnomad_gene_scores' not in config_dict:

		print ('You need to add the location of the gnomad_gene_scores in the config file to use the gnomad-constraint-scores option.')
		return False

	if args.patient_hpos != None and 'hpo_file' not in config_dict:

		print ('Cannot add patient HPO terms as we do not have a HPO gene map location specified in the config file.')
		return False

	if args.patient_hpos != None and ('HPOCount' not in config_dict['final_fields_trio'] or 'HPOCount' not in config_dict['final_fields_single']) :

		print ('Cannot add patient HPO terms and not put HPOCount in the final fields variable.')
		return False

	return True


def fix_gnomad(df, column):
	"""
	Fix gnomad columns when we have multiple hits.

	return largest of the two.

	"""
	
	gnomad =  df[column]
	
	# If we have blank or None return NOne
	if gnomad == '.' or gnomad == None:
		
		return None
	
	# If for whatever reason we get two AFs
	if '&' in gnomad:
		
		gnomad = gnomad.split('&')
	
		if '.' in gnomad:

			gnomad = [x for x in gnomad if x != '.']

		# If the list is empty after removing . 
		if len(gnomad) == 0:

			return None

		gnomad = [float(x) for x in gnomad]
		
		# Else return the largest (as this is the max_Af field)
		gnomad = max(gnomad)
	
	else:
		
		try:
		
			return float(gnomad)
		
		except:
			
			return None



def add_to_db(df, db_path):
	"""
	Add variants to in house database
	"""

	conn = sqlite3.connect(db_path)
	c = conn.cursor()

	# Create table if it does not exist
	c.execute('create table if not exists Variants (sample_id text, run_id text, variant text, gt text)')

	# Add each variant
	for row in df.itertuples():
		
		sample_id = row.SampleId
		run_id = row.RunId
		variant = row.VariantId
		gt = row.Genotype
		
		data = (sample_id, run_id, variant, gt)
		
		# Check whether the variant is already in for this sample
		c.execute("SELECT * FROM Variants WHERE sample_id = ? AND variant = ?", (data[0], data[2]))
		my_variant =  c.fetchall()

		# If not add the variant
		if len(my_variant) == 0:
		
			c.execute("INSERT INTO Variants VALUES (?,?,?,?)", data)
		
	conn.commit()
	conn.close()

	return None


def get_db_frequencies(db_path):
	"""
	Given a path to the in house database return a dict with the count of how \
	many samples we have seen that variant in

	"""
	
	conn = sqlite3.connect(db_path)
	
	c = conn.cursor()
	
	c.execute("SELECT * FROM Variants")
	
	all_variants = c.fetchall()
	
	var_df = pd.DataFrame(all_variants, columns=['sample_id', 'run_id', 'variant', 'gt'])
	
	conn.commit()
	conn.close()
	
	var_sample_dict = var_df.groupby('variant')['sample_id'].apply(lambda x: "|".join(x)).to_dict()

	return var_sample_dict
	

def add_in_house_count(df, freq_dict):
	"""
	Add how many times we have seen this variant in our own patients
	"""
	
	if df['VariantId'] in freq_dict:
		
		
		samples = freq_dict[df['VariantId']].split('|')
				
		if df['SampleId'] in samples:
			
			samples.remove(df['SampleId'])
						
		else:
			
			samples = samples
			
		return len(samples)
			
	
	else:
		
		return 0














