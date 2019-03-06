from utils.utils import *
from utils.inheritance_utils import *
import pandas as pd
import csv
import requests
import logging
import argparse

########################################################################################################################################################
# Set up Logger
########################################################################################################################################################

version = '0.0.1'

logger = logging.getLogger('germline_variant_filter')
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    '%(levelname)s\t%(asctime)s\t%(name)s\t%(message)s'
)
handler.setFormatter(formatter)
logger.addHandler(handler)

########################################################################################################################################################
# Parse Arguments
########################################################################################################################################################


parser = argparse.ArgumentParser(
	formatter_class=argparse.RawTextHelpFormatter,
	description='A program to filter variants from germline sequencing pipelines. See readme for details.')

	
parser.add_argument('--config', type=str, nargs=1, required=True,
					help='Filepath to YAML config file. See config/ directory for an example of what these look like.')

parser.add_argument('--ped', type=str, nargs=1, required=True,
					help='Filepath to PED file.')

parser.add_argument('--input', type=str, nargs=1, required=True,
					help='Filepath to the input CSV file. See readme for details on the required format.')

parser.add_argument('--panelapp', action='store_true',
					help='Whether to add PanelApp annotations. Default = False.')

parser.add_argument('--local-panel-app-dump', type=str, nargs=1,
					help='Filepath local panelapp store.')

parser.add_argument('--csq', type=str, nargs=1, required=True,
					help='The VEP CSQ string. For example Allele|Consequence|IMPACT|SYMBOL...')

parser.add_argument('--spliceai', action='store_true',
					help='Attempt to parse SpliceAI annotations added by VEP. Default = False.')

parser.add_argument('--smart-synonymous', action='store_true',
					help='Smart synonymous variant filtering. See readme for details. Default = False.')

parser.add_argument('--add-ccrs', action='store_true',
					help='Add the CCR annotations from VEP. Default = False.')

parser.add_argument('--gnomad-constraint-scores', action='store_true',
					help='Add the per gene gnomad constraint scores. Default = False.')

parser.add_argument('--patient-hpos', type=str, nargs=1,
					help='Filepath to file containing patient HPO terms. See examples/ directory for information on the format of this file.')

parser.add_argument('--worksheet', type=str, nargs=1, required =True,
					help='The worksheet ID.')

parser.add_argument('--results-dir', type=str, nargs=1, required =True,
					help='Where to put the results.')

parser.add_argument('--unaffected-parent-filter', action='store_true',
					help='Apply the unaffected parent filter function. EXPERIMENTAL')

args = parser.parse_args()

config = args.config[0]
ped_file = args.ped[0]
csv_file = args.input[0]
csq_desc = args.csq[0]
csq_desc = csq_desc.split('|')
parse_splice_ai = args.spliceai
smart_synonymous_filtering = args.smart_synonymous
add_ccrs = args.add_ccrs
add_gnomad_constaint_scores = args.gnomad_constraint_scores
add_panel_app_info = args.panelapp
worksheet = args.worksheet[0]
results_dir = args.results_dir[0]
unaffected_parent_filter = args.unaffected_parent_filter


if args.local_panel_app_dump != None:

	local_panel_app_dump = args.local_panel_app_dump[0]
	use_local_panel_app_dump = True

else:

	use_local_panel_app_dump = False

if args.patient_hpos != None:

	add_hpo = True
	patient_hpos = args.patient_hpos[0]

else:

	add_hpo = False

########################################################################################################################################################
# Parse Config files
########################################################################################################################################################

config_dict = parse_config(config)

if are_arguments_valid(args, config_dict) == False:

	raise Exception('Invalid command line options.')

ped_dict = parse_ped_file(ped_file)
vep_fields = config_dict['vep_fields']
samples = list(ped_dict.keys())
default_cutoff_gnomad_genomes = config_dict['default_cutoff_gnomad_genomes']
default_cutoff_gnomad_exomes = config_dict['default_cutoff_gnomad_exomes']
splice_ai_cutoff = config_dict['splice_ai_cutoff']
consequence_severity = config_dict['consequence_severity']
clin_sig_words = config_dict['clin_sig_words']
to_keep_consequences = config_dict['to_keep_consequences']

min_dp = config_dict['min_dp']
min_gq = config_dict['min_gq']

min_parental_depth_dn = config_dict['min_parental_depth_dn']
min_parental_gq_dn = config_dict['min_parental_gq_dn']
min_parental_depth_uid = config_dict['min_parental_depth_uid']
min_parental_gq_uid = config_dict['min_parental_gq_uid']

gnomad_gene_scores = config_dict['gnomad_gene_scores']
panel_app_dump_max_time = config_dict['panel_app_dump_max_time']

if add_hpo == True:

	hpo_file = config_dict['hpo_file']

other_gnomadg = config_dict['other_gnomadg']
other_gnomade = config_dict['other_gnomade']
other_ac = config_dict['other_ac']

compound_het_gnomadg = config_dict['compound_het_gnomadg']
compound_het_gnomade = config_dict['compound_het_gnomade']

dom_autosomal_gnomadg = config_dict['dom_autosomal_gnomadg']
dom_autosomal_gnomade = config_dict['dom_autosomal_gnomade']
dom_autosomal_ac = config_dict['dom_autosomal_ac']

dom_x_female_gnomadg = config_dict['dom_x_female_gnomadg']
dom_x_female_gnomade = config_dict['dom_x_female_gnomade']
dom_x_female_ac = config_dict['dom_x_female_ac']

x_linked_male_gnomadg = config_dict['x_linked_male_gnomadg']
x_linked_male_gnomade = config_dict['x_linked_male_gnomade']
x_linked_male_ac = config_dict['x_linked_male_ac']

y_linked_male_gnomadg = config_dict['y_linked_male_gnomadg']
y_linked_male_gnomade = config_dict['y_linked_male_gnomade']
y_linked_male_ac = config_dict['y_linked_male_ac']

reccessive_autosomal_gnomadg = config_dict['reccessive_autosomal_gnomadg']
reccessive_autosomal_gnomade = config_dict['reccessive_autosomal_gnomade']
reccessive_autosomal_ac = config_dict['reccessive_autosomal_ac']

reccessive_x_female_gnomadg = config_dict['reccessive_x_female_gnomadg']
reccessive_x_female_gnomade = config_dict['reccessive_x_female_gnomade']
reccessive_x_female_ac = config_dict['reccessive_x_female_ac']

mito_gnomadg = config_dict['mito_gnomadg']
mito_gnomade = config_dict['mito_gnomade']

de_novo_gnomadg = config_dict['de_novo_gnomadg']
de_novo_gnomade = config_dict['de_novo_gnomade']
de_novo_ac = config_dict['de_novo_ac']

upi_gnomadg = config_dict['upi_gnomadg']
upi_gnomade = config_dict['upi_gnomade']

final_fields_trio  =config_dict['final_fields_trio']
final_fields_single  =config_dict['final_fields_single']

gt_depth_tag = config_dict['gt_depth_tag']

########################################################################################################################################################
# Load misc data
########################################################################################################################################################

# Load Gnomad gene constraint scores
if add_gnomad_constaint_scores == True:

	logger.info('Parsing gnomad constraint scores.')
	gnomad_scores_df = pd.read_csv(gnomad_gene_scores, index_col=0)
	gnomad_scores_dict = gnomad_scores_df.to_dict()

# If we want to add panel app data
if add_panel_app_info == True:

	if use_local_panel_app_dump == True:

		logger.info('Reading Local Panel App data.')

		try:

			panel_app_dict = parse_panel_app_dump(local_panel_app_dump)
			
		except:

			logger.info('Could not local Panel App data.')
			panel_app_dict = {}
	else:

		panel_app_dict = {}

# If we want to annotate variants with HPO matches
if add_hpo == True:

	logger.info('Parsing HPO Gene Map File.')

	hpo_dict = parse_hpo_file(hpo_file)

	try:

		patient_hpos = pd.read_csv(patient_hpos, sep='\t').to_dict(orient='list')

	except:

		logger.warning('Could not read HPO gene map')

########################################################################################################################################################
# Initial Preprocessing of the Data
########################################################################################################################################################

logger.info('Parsing CSV into dataframe.')

# Parse CSV into dataframe
df = pd.read_csv(csv_file, sep='\t', dtype={'CHROM': object})

# Filter out variants that fail variant level QC
df = df[(df['FILTER'] == 'PASS') | (df['FILTER'] == '') | (pd.isna(df['FILTER']))]

# For each sample in the PED file create a column which specifies whether the variant is relevant for that sample
for sample in samples:
	df[sample + '_is_relevant'] = df.apply(select_variants_for_sample, args=(sample,min_dp, min_gq, gt_depth_tag), axis=1)

# Fix column names
df.columns = fix_column_names(df.columns)

# Parse CSQ data - putting each consequence block on its own line.
vep_df = split_vep_transcripts(df, csq_desc, vep_fields, list(df.columns))


########################################################################################################################################################
# Initial Frequency Filter
########################################################################################################################################################

logger.info('Filtering on default filtering settings.')

#Parse columns where we have two results e.g 0.001&0.3
vep_df['gnomADg_AF_POPMAX'] = vep_df.apply(fix_gnomad, axis=1, args=('gnomADg_AF_POPMAX',))
vep_df['gnomADe_AF_POPMAX'] = vep_df.apply(fix_gnomad, axis=1, args=('gnomADe_AF_POPMAX',))
vep_df['gnomADg_AF_POPMAX'] = pd.to_numeric(vep_df['gnomADg_AF_POPMAX'])
vep_df['gnomADe_AF_POPMAX'] = pd.to_numeric(vep_df['gnomADe_AF_POPMAX'])


vep_df.fillna(value = {'gnomADg_AF_POPMAX':0.0, 'gnomADe_AF_POPMAX':0.0}, inplace=True)

# Filter on gnomad genomes and exomes - if data is missing or we have less than largest cutoff  e.g. (1%)
vep_df = vep_df[((vep_df['gnomADg_AF_POPMAX'] <= default_cutoff_gnomad_genomes) | (pd.isna(vep_df['gnomADg_AF_POPMAX']) )) &
			   ((vep_df['gnomADe_AF_POPMAX'] <= default_cutoff_gnomad_exomes ) | (pd.isna(vep_df['gnomADe_AF_POPMAX'])))]

# Also create the variant key e.g.12:12345A>G
vep_df['VariantId'] = vep_df.apply(get_variant_key,axis=1)

########################################################################################################################################################
# Process SpliceAI columns if requested
########################################################################################################################################################

if parse_splice_ai == True:

	logger.info('Fixing SpliceAI columns.')

	# Apply fix for splice AI columns

	vep_df['SpliceAI_DS_AG'] = vep_df.apply(fix_splice_ai, axis=1, args=('SpliceAI_DS_AG',))
	vep_df['SpliceAI_DS_AL'] = vep_df.apply(fix_splice_ai, axis=1, args=('SpliceAI_DS_AL',))
	vep_df['SpliceAI_DS_DG'] = vep_df.apply(fix_splice_ai, axis=1, args=('SpliceAI_DS_DG',))
	vep_df['SpliceAI_DS_DL'] = vep_df.apply(fix_splice_ai, axis=1, args=('SpliceAI_DS_DL',))

	vep_df['SpliceAI_DS_AG'] = pd.to_numeric(vep_df['SpliceAI_DS_AG'])
	vep_df['SpliceAI_DS_AL'] = pd.to_numeric(vep_df['SpliceAI_DS_AL'])
	vep_df['SpliceAI_DS_DG'] = pd.to_numeric(vep_df['SpliceAI_DS_DG'])
	vep_df['SpliceAI_DS_DL'] = pd.to_numeric(vep_df['SpliceAI_DS_DL'])

	vep_df['has_affect_on_splicing'] = vep_df.apply(has_affect_on_splicing, axis=1, args=(splice_ai_cutoff,))
	vep_df['any_has_splicing_affect'] = vep_df.groupby('VariantId')['has_affect_on_splicing'].transform(any_has_splicing_affect)

########################################################################################################################################################
# Consequence Filtering
########################################################################################################################################################

logger.info('Filtering on Consequence.')

# Get worst consequence (in any feature)
vep_df['WorstConsequence'] = vep_df.groupby('VariantId')['Consequence'].transform(get_worst_consequence, consequence_severity)

# Has the variant got a relevant clinical consequence e.g. Pathogenic?
vep_df['has_important_clinsig'] = vep_df.apply(has_important_clinsig, axis=1, args=(clin_sig_words,))

# Apply consequence filter
vep_df['consequence_filter'] = vep_df.apply(consequence_filter, axis=1, args=(to_keep_consequences,))

if smart_synonymous_filtering == True:

	vep_df = vep_df[(vep_df['consequence_filter'] == False) | ((vep_df['WorstConsequence'] == 'synonymous_variant') & ((vep_df['has_important_clinsig'] == True) | (vep_df['any_has_splicing_affect'] == True))) ]

else:

	vep_df = vep_df[(vep_df['consequence_filter'] == False)]


########################################################################################################################################################
# Per Sample Processing
########################################################################################################################################################

# List of samples in which we have no variants
no_variants_samples = []

for sample in samples:

	logger.info(f'Processing sample: {sample}')

	proband_in_trio = is_proband_in_trio(sample, ped_dict)

	logger.info(f'{sample} is proband in trio:\t{proband_in_trio}')

	sample_sex = ped_dict[sample]['sex']

	logger.info(f'{sample} sex:\t{sample_sex}')

	if sample_sex == '0':

		logger.warning(f'{sample}: sex is unknown - downstream calculations will assume patient is Male. We reccomend rerunning program when sex is known.')
		sample_sex = 'Unknown'

	elif sample_sex =='1':

		sample_sex = 'Male'

	elif sample_sex == '2':

		sample_sex = 'Female'

	else:

		logger.warning(f'{sample}: sex is unknown - downstream calculations will assume patient is Male. We reccomend rerunning program when sex is known.')
		sample_sex = 'Unknown'

	# Get variants relevant to this sample
	sample_df = vep_df[vep_df[f'sample_{sample}_is_relevant'] == True]

	logger.info(f'{sample}: Found {sample_df.shape[0]} relevant variants in sample.')

	# If there are no variants for this sample
	if sample_df.shape[0] == 0:

		logger.warning(f'{sample}: sample has 0 variants - cannot carry out downstream filtering.')
		no_variants_samples.append(sample)
		continue

	# Create compound HET dict
	compound_het_dict = sample_df.groupby('Feature').count()['CHROM'].to_dict()
	compound_het_dict[None] = 0

	# Seperate workflows for trios and single samples - annotate each variant with inheritance pattern
	if proband_in_trio == True:

		mother = ped_dict[sample]['maternalID']
		father = ped_dict[sample]['paternalID']

		sample_df['Workflow'] = sample_df.apply(annotate_workflow_trio, axis=1, args=(sample, mother, father, sample_sex, compound_het_dict, min_parental_depth_dn, min_parental_gq_dn, min_parental_depth_uid, min_parental_gq_uid, gt_depth_tag))

	else:

		sample_df['Workflow'] = sample_df.apply(annotate_workflow_single, axis=1, args=(sample, sample_sex, compound_het_dict))


	# Filter on least restrictive workflow
	wf_restrictiveness = config_dict['wf_restrictiveness']

	workflows = list(sample_df['Workflow'].value_counts().index)

	master_sample_df = pd.DataFrame(columns =sample_df.columns)


	# For each workflow filter on the specific workflow settings
	for workflow in workflows:
		
		workflow_df = sample_df[sample_df['Workflow'] == workflow]
		
		# Get the least restrictive workflow if we have multiple annotations e.g DOMINANT_AUTOSOMAL|DENOVO_HC
		workflow_list = workflow.split('|')
		
		least_restrictive_idx = 9999999
		
		for wf in workflow_list:
			
			if wf_restrictiveness.index(wf) < least_restrictive_idx:
				
				least_restrictive_idx = wf_restrictiveness.index(wf)
				
		wf_to_use = wf_restrictiveness[least_restrictive_idx]
		
		if wf_to_use == 'OTHER':
			
			workflow_df = workflow_df[((workflow_df['gnomADg_AF_POPMAX'] < other_gnomadg) | (pd.isna(workflow_df['gnomADg_AF_POPMAX']) )) &
				   ((workflow_df['gnomADe_AF_POPMAX'] < other_gnomade ) | (pd.isna(workflow_df['gnomADe_AF_POPMAX']))) &
				   (workflow_df['AC'] < other_ac)
				   ]

		elif wf_to_use == 'COMPOUND_HET':

			workflow_df = workflow_df[((workflow_df['gnomADg_AF_POPMAX'] < compound_het_gnomadg) | (pd.isna(workflow_df['gnomADg_AF_POPMAX']) )) &
				   ((workflow_df['gnomADe_AF_POPMAX'] < compound_het_gnomade ) | (pd.isna(workflow_df['gnomADe_AF_POPMAX'])))
				   ]

		elif wf_to_use == 'UNIPARENTAL_ISODISOMY':

			workflow_df = workflow_df[((workflow_df['gnomADg_AF_POPMAX'] < upi_gnomadg) | (pd.isna(workflow_df['gnomADg_AF_POPMAX']) )) &
				   ((workflow_df['gnomADe_AF_POPMAX'] < upi_gnomade ) | (pd.isna(workflow_df['gnomADe_AF_POPMAX'])))
				   ]

		elif wf_to_use == 'MITOCHONDRIAL':

			workflow_df = workflow_df[((workflow_df['gnomADg_AF_POPMAX'] < mito_gnomadg) | (pd.isna(workflow_df['gnomADg_AF_POPMAX']) )) &
				   ((workflow_df['gnomADe_AF_POPMAX'] < mito_gnomade ) | (pd.isna(workflow_df['gnomADe_AF_POPMAX'])))
				   ]			   

		elif wf_to_use == 'RECCESSIVE_X_FEMALE':
			
			workflow_df = workflow_df[((workflow_df['gnomADg_AF_POPMAX'] < reccessive_x_female_gnomadg) | (pd.isna(workflow_df['gnomADg_AF_POPMAX']) )) &
				   ((workflow_df['gnomADe_AF_POPMAX'] < reccessive_x_female_gnomade ) | (pd.isna(workflow_df['gnomADe_AF_POPMAX']))) &
				   (workflow_df['AC'] < reccessive_x_female_ac)
				   ]
			
		elif wf_to_use == 'RECCESSIVE_AUTOSOMAL':
			
			workflow_df = workflow_df[((workflow_df['gnomADg_AF_POPMAX'] < reccessive_autosomal_gnomadg) | (pd.isna(workflow_df['gnomADg_AF_POPMAX']) )) &
				   ((workflow_df['gnomADe_AF_POPMAX'] < reccessive_autosomal_gnomade ) | (pd.isna(workflow_df['gnomADe_AF_POPMAX']))) &
				   (workflow_df['AC'] < reccessive_autosomal_ac)
				   ]      
		
		elif wf_to_use == 'X_LINKED_MALE':
			
			workflow_df = workflow_df[((workflow_df['gnomADg_AF_POPMAX'] < x_linked_male_gnomadg) | (pd.isna(workflow_df['gnomADg_AF_POPMAX']) )) &
				   ((workflow_df['gnomADe_AF_POPMAX'] < x_linked_male_gnomade ) | (pd.isna(workflow_df['gnomADe_AF_POPMAX']))) &
				   (workflow_df['AC'] < x_linked_male_ac)
				   ]
		elif wf_to_use == 'Y_LINKED_MALE':
			
			workflow_df = workflow_df[((workflow_df['gnomADg_AF_POPMAX'] < y_linked_male_gnomadg) | (pd.isna(workflow_df['gnomADg_AF_POPMAX']) )) &
				   ((workflow_df['gnomADe_AF_POPMAX'] < y_linked_male_gnomade ) | (pd.isna(workflow_df['gnomADe_AF_POPMAX']))) &
				   (workflow_df['AC'] < y_linked_male_ac)
				   ]

		elif wf_to_use == 'DOMINANT_AUTOSOMAL':
			
			workflow_df = workflow_df[((workflow_df['gnomADg_AF_POPMAX'] < dom_autosomal_gnomadg) | (pd.isna(workflow_df['gnomADg_AF_POPMAX']) )) &
				   ((workflow_df['gnomADe_AF_POPMAX'] < dom_autosomal_gnomade ) | (pd.isna(workflow_df['gnomADe_AF_POPMAX']))) &
				   (workflow_df['AC'] < dom_autosomal_ac)
				   ]

		elif wf_to_use == 'DOMINANT_X_FEMALE':
			
			workflow_df = workflow_df[((workflow_df['gnomADg_AF_POPMAX'] < dom_x_female_gnomadg) | (pd.isna(workflow_df['gnomADg_AF_POPMAX']) )) &
				   ((workflow_df['gnomADe_AF_POPMAX'] < dom_x_female_gnomade ) | (pd.isna(workflow_df['gnomADe_AF_POPMAX']))) &
				   (workflow_df['AC'] < dom_x_female_ac)
				   ]

		elif wf_to_use == 'DE_NOVO_HC' or wf_to_use == 'DE_NOVO_LC':

			workflow_df = workflow_df[((workflow_df['gnomADg_AF_POPMAX'] < de_novo_gnomadg) | (pd.isna(workflow_df['gnomADg_AF_POPMAX']) )) &
				   ((workflow_df['gnomADe_AF_POPMAX'] < de_novo_gnomade ) | (pd.isna(workflow_df['gnomADe_AF_POPMAX']))) &
				   (workflow_df['AC'] < de_novo_ac)
				   ]		
				
		else:
			
			raise Exception('Unknown workflow')
			
		master_sample_df = pd.concat([master_sample_df, workflow_df])

	# If there are no variants left for this sample
	if master_sample_df.shape[0] == 0:

		logger.warning(f'{sample}: sample has 0 variants - cannot carry out downstream filtering.')
		no_variants_samples.append(sample)
		continue 

	logger.info(f"{sample}: variants remaining after filtering: {master_sample_df.groupby('VariantId').count().shape[0]}")

	# If we want to add the constrained coding regions
	if add_ccrs == True:

		master_sample_df['CCR_percentile'] = master_sample_df.apply(fix_ccrs,axis=1)

	# If we want to add the gnomad per gene constraint scores
	if add_gnomad_constaint_scores == True:

		master_sample_df['pLI'] = master_sample_df.apply(annotate_with_gnomad_scores, axis=1, args=(gnomad_scores_dict,'pLI'))
		master_sample_df['oe_lof'] = master_sample_df.apply(annotate_with_gnomad_scores, axis=1, args=(gnomad_scores_dict, 'oe_lof'))
		master_sample_df['oe_lof_lower'] = master_sample_df.apply(annotate_with_gnomad_scores, axis=1, args=(gnomad_scores_dict, 'oe_lof_lower'))
		master_sample_df['oe_lof_upper'] = master_sample_df.apply(annotate_with_gnomad_scores, axis=1, args=(gnomad_scores_dict, 'oe_lof_upper'))

	# Add gene information from panel app
	if add_panel_app_info == True:

		master_sample_df['DiseaseName'] = master_sample_df.apply(apply_panel_app_data_disease, axis=1,args =(panel_app_dict,panel_app_dump_max_time,))
		master_sample_df['ModeOfInheritance'] = master_sample_df.apply(apply_panel_app_data_inheritance, axis=1,args =(panel_app_dict,panel_app_dump_max_time,))

	if add_hpo == True:

		if sample in patient_hpos:

			master_sample_df['HPOCount'] = master_sample_df.apply(annotate_hpo, axis=1,args=(patient_hpos[sample], hpo_dict))
			master_sample_df['HPOCountMax'] = master_sample_df.groupby('VariantId')['HPOCount'].transform('max')

		else:
			
			master_sample_df['HPOCount'] = 'NA'
			master_sample_df['HPOCountMax'] = 'NA'
			logger.warning(f'{sample}: Could not find any HPO terms for this sample in file.')


	# Check each variant has at least one pick flag
	pick_dict = master_sample_df[['VariantId','PICK']].groupby('VariantId').sum().to_dict()

	# Process CSV for saving to file.
	master_sample_df['SampleId'] = sample
	master_sample_df['RunId'] = worksheet
	master_sample_df['Genotype'] = master_sample_df.apply(get_genotype, axis=1, args=(sample,))
	master_sample_df['Proband'] = master_sample_df['sample_' + sample + '_GT']

	# Change column names if we have trio samples
	if proband_in_trio == True:

		master_sample_df['Father'] = master_sample_df['sample_' + father + '_GT']
		master_sample_df['Mother'] = master_sample_df['sample_' + mother + '_GT']

	master_sample_df['Gene'] = master_sample_df['SYMBOL']
	master_sample_df['Transcript'] = master_sample_df['Feature']
	master_sample_df['Exon'] = master_sample_df.apply(fix_exon, axis=1)
	master_sample_df['Intron'] = master_sample_df.apply(fix_intron, axis=1)
	master_sample_df['HGVSc'] = master_sample_df.apply(get_hgvsc, axis=1)
	master_sample_df['HGVSp'] = master_sample_df.apply(get_hgvsp, axis=1)

	if proband_in_trio == True and unaffected_parent_filter == True:

		logger.info('Applying unaffected parent filter - EXPERIENTIAL - if one or more of the parents are unaffected we will get errors!')

		master_sample_df = apply_unaffected_parent_filter(master_sample_df)

	# Check every variant has at least one PICK flag
	master_sample_df['Pick'] = master_sample_df.apply(check_picks, axis=1, args=(pick_dict,))

	if proband_in_trio == True:

		with open(f'{results_dir}/{sample}.csv', 'w') as f:
			f.write(f'#Variant Germline Filter Version {version}|Proband={sample}|father={father}|mother={mother}|proband_sex={sample_sex}\n')

		master_sample_df[final_fields_trio].to_csv(f'{results_dir}/{sample}.csv', sep='\t', float_format='%.6f', mode='a', index=False)

	else:

		with open(f'{results_dir}/{sample}.csv', 'w') as f:
			f.write(f'#Variant Germline Filter Version {version}|Proband={sample}|proband_sex={sample_sex}\n')

		master_sample_df[final_fields_single].to_csv(f'{results_dir}/{sample}.csv', sep='\t', float_format='%.6f', mode='a', index=False)

########################################################################################################################################################
# Clean up
########################################################################################################################################################

# Write panel app data to local dump if requested
if use_local_panel_app_dump:

	panel_app_df = pd.DataFrame(panel_app_dict)

	panel_app_df.transpose().to_csv(local_panel_app_dump, sep='\t')


# Write blank files for samples with no variants in.
for sample in no_variants_samples:

	with open(f'{results_dir}/{sample}.csv', 'w') as f:
			f.write(f'')









			   