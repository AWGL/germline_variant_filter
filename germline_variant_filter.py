from utils.utils import *
from utils.inheritance_utils import *
import pandas as pd
import csv
import requests


########################################################################################################################################################
# Parse Arguments
########################################################################################################################################################

config = 'config/config.yaml'
ped_file = 'test/pedigree.ped'
csv_file = '../180126/180126_D00501_0177_BH2LC7BCX2_filtered_annotated_roi_norm_vep.csv'
local_panel_app_dump = '../panelapp.csv'


csq_desc = 'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|REFSEQ_MATCH|GIVEN_REF|USED_REF|BAM_EDIT|SOURCE|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|gnomADg|gnomADg_AF_AFR|gnomADg_AF_AMR|gnomADg_AF_ASJ|gnomADg_AF_EAS|gnomADg_AF_FIN|gnomADg_AF_NFE|gnomADg_AF_OTH|gnomADg_AF_POPMAX|gnomADe|gnomADe_AF_POPMAX|ccrs|SpliceAI|SpliceAI_DS_AG|SpliceAI_DS_AL|SpliceAI_DS_DG|SpliceAI_DS_DL|SpliceAI_SYMBOL'
csq_desc = csq_desc.split('|')

parse_splice_ai = True
smart_synonymous_filtering = True
add_ccrs = True
add_gnomad_constaint_scores = True
add_panel_app_info = True
use_local_panel_app_dump = True

########################################################################################################################################################
# Parse Config files
########################################################################################################################################################

config_dict = parse_config(config)
ped_dict = parse_ped_file(ped_file)

vep_fields = config_dict['vep_fields']
samples = list(ped_dict.keys())[0:1]

default_cutoff_gnomad_genomes = config_dict['default_cutoff_gnomad_genomes']
default_cutoff_gnomad_exomes = config_dict['default_cutoff_gnomad_exomes']
splice_ai_cutoff = config_dict['splice_ai_cutoff']
consequence_severity = config_dict['consequence_severity']
clin_sig_words = config_dict['clin_sig_words']
to_keep_consequences = config_dict['to_keep_consequences']

min_parental_depth_dn = config_dict['min_parental_depth_dn']
min_parental_gq_dn = config_dict['min_parental_gq_dn']
min_parental_depth_uid = config_dict['min_parental_depth_uid']
min_parental_gq_uid = config_dict['min_parental_gq_uid']

gnomad_gene_scores = config_dict['gnomad_gene_scores']
panel_app_dump_max_time = config_dict['panel_app_dump_max_time']

########################################################################################################################################################
# Load miscallenous data
########################################################################################################################################################


if add_gnomad_constaint_scores == True:

	gnomad_scores_df = pd.read_csv(gnomad_gene_scores, index_col=0)
	gnomad_scores_dict = gnomad_scores_df.to_dict()


if add_panel_app_info == True:

	if use_local_panel_app_dump == True:

		try:

			panel_app_dict = parse_panel_app_dump(local_panel_app_dump)
			
		except:

			panel_app_dict = {}
	else:

		panel_app_dict = {}

########################################################################################################################################################
# Main Program
########################################################################################################################################################

# Parse CSV into dataframe
df = pd.read_csv(csv_file, sep='\t', dtype={'CHROM': object})

# Filter out variants that fail variant level QC
df = df[df['FILTER'] == 'PASS']

# For each sample in the PED file create a column which specifies whether the variant is relevant for that sample
for sample in samples:
	df[sample + '_is_relevant'] = df.apply(select_variants_for_sample, args=(sample,), axis=1)

# Fix column names
df.columns = fix_column_names(df.columns)

# Parse CSQ data - putting each consequence block on its own line.
vep_df = split_vep_transcripts(df, csq_desc, vep_fields, list(df.columns))


########################################################################################################################################################
# Initial Frequency Filter
########################################################################################################################################################

# Convert columns to numerical
vep_df['gnomADg_AF_POPMAX'] = pd.to_numeric(vep_df['gnomADg_AF_POPMAX'])
vep_df['gnomADe_AF_POPMAX'] = pd.to_numeric(vep_df['gnomADe_AF_POPMAX'])

# Filter on gnomad genomes and exomes - if data is missing or we have less than largest cutoff  e.g. (1%)
vep_df = vep_df[((vep_df['gnomADg_AF_POPMAX'] <= default_cutoff_gnomad_genomes) | (pd.isna(vep_df['gnomADg_AF_POPMAX']) )) &
			   ((vep_df['gnomADe_AF_POPMAX'] <= default_cutoff_gnomad_exomes ) | (pd.isna(vep_df['gnomADe_AF_POPMAX'])))]


# Also create the variant key e.g.12:12345A>G
vep_df['variant'] = vep_df.apply(get_variant_key,axis=1)

########################################################################################################################################################
# Process SpliceAI columns if requested
########################################################################################################################################################

if parse_splice_ai == True:

	# Apply fix for splice AI columns

	vep_df['fixed_SpliceAI_DS_AG'] = vep_df.apply(fix_splice_ai, axis=1, args=('SpliceAI_DS_AG',))
	vep_df['fixed_SpliceAI_DS_AL'] = vep_df.apply(fix_splice_ai, axis=1, args=('SpliceAI_DS_AL',))
	vep_df['fixed_SpliceAI_DS_DG'] = vep_df.apply(fix_splice_ai, axis=1, args=('SpliceAI_DS_DG',))
	vep_df['fixed_SpliceAI_DS_DL'] = vep_df.apply(fix_splice_ai, axis=1, args=('SpliceAI_DS_DL',))

	vep_df['fixed_SpliceAI_DS_AG'] = pd.to_numeric(vep_df['fixed_SpliceAI_DS_AG'])
	vep_df['fixed_SpliceAI_DS_AL'] = pd.to_numeric(vep_df['fixed_SpliceAI_DS_AL'])
	vep_df['fixed_SpliceAI_DS_DG'] = pd.to_numeric(vep_df['fixed_SpliceAI_DS_DG'])
	vep_df['fixed_SpliceAI_DS_DL'] = pd.to_numeric(vep_df['fixed_SpliceAI_DS_DL'])

	vep_df['has_affect_on_splicing'] = vep_df.apply(has_affect_on_splicing, axis=1, args=(splice_ai_cutoff,))
	vep_df['any_has_splicing_affect'] = vep_df.groupby('variant')['has_affect_on_splicing'].transform(any_has_splicing_affect)

########################################################################################################################################################
# Consequence Filtering
########################################################################################################################################################

# Get worst consequence (in any feature)
vep_df['worst_consequence'] = vep_df.groupby('variant')['Consequence'].transform(get_worst_consequence, consequence_severity)

# Has the variant got a relevant clinical consequence e.g. Pathogenic?
vep_df['has_important_clinsig'] = vep_df.apply(has_important_clinsig, axis=1, args=(clin_sig_words,))

# Apply consequence filter
vep_df['consequence_filter'] = vep_df.apply(consequence_filter, axis=1, args=(to_keep_consequences,))

if smart_synonymous_filtering == True:

	# TODO: Add check whether we have the data for this!
	vep_df = vep_df[(vep_df['consequence_filter'] == False) | ((vep_df['worst_consequence'] == 'synonymous_variant') & ((vep_df['has_important_clinsig'] == True) | (vep_df['has_affect_on_splicing'] == True))) ]

else:

	vep_df = vep_df[(vep_df['consequence_filter'] == False)]


########################################################################################################################################################
# Per Sample Processing
########################################################################################################################################################

for sample in samples:

	proband_in_trio = is_proband_in_trio(sample, ped_dict)

	print (f'{sample} is proband in trio:\t{proband_in_trio}')

	sample_sex = ped_dict[sample]['sex']

	print (f'{sample} sex:\t{sample_sex}')

	if sample_sex == '0':

		print (f'{sample} WARNING: sex is unknown - downstream calculations will assume patient is Male. We reccomend rerunning program when sex is known.')

	elif sample_sex =='1':

		sample_sex = 'Male'

	elif sample_sex == '1':

		sample_sex = 'Female'

	else:

		sample_sex = '0'

	# Get variants relevant to this sample
	sample_df = vep_df[vep_df[f'sample_{sample}_is_relevant'] == True]

	print(f'{sample} INFO: We have found {sample_df.shape[0]} variants in sample {sample}.')

	if sample_df.shape[0] == 0:

		print (f'{sample} WARNING: sample has 0 variants - cannot carry out downstream filtering.')
		continue


	# Create compound HET dict
	compound_het_dict = sample_df.groupby('Feature').count()['CHROM'].to_dict()
	compound_het_dict[None] = 0

	# Seperate workflows for trios and single samples
	if is_proband_in_trio == True:

		mother = ped_dict[sample]['maternalID']
		father = ped_dict[sample]['paternalID']

		sample_df['workflow'] = sample_df.apply(annotate_workflow_trio, axis=1, args=(sample, mother, father, sample_sex, compound_het_dict, min_parental_depth_dn, min_parental_gq_dn, min_parental_depth_uid, min_parental_gq_uid))

	else:

		sample_df['workflow'] = sample_df.apply(annotate_workflow_single, axis=1, args=(sample, sample_sex, compound_het_dict))


	# Filter on least restrictive workflow

	wf_restrictiveness = ['OTHER', 'UNIPARENTAL_ISODISOMY', 'COMPOUND_HET', 'RECCESSIVE_SEX', 'RECESSIVE_AUTOSOMAL', 'DOMINANT_SEX', 'DOMINANT_AUTOSOMAL', 'DENOVO_HC', 'DENOVO_LC' ]

	workflows = list(sample_df['workflow'].value_counts().index)

	master_sample_df = pd.DataFrame(columns =sample_df.columns)

	for workflow in workflows:
		
		workflow_df = sample_df[sample_df['workflow'] == workflow]
		
		workflow_list = workflow.split('|')
		
		least_restrictive_idx = 99999
		
		for wf in workflow_list:
			
			if wf_restrictiveness.index(wf) < least_restrictive_idx:
				
				least_restrictive_idx = wf_restrictiveness.index(wf)
				
		wf_to_use = wf_restrictiveness[least_restrictive_idx]
		
		if wf_to_use == 'OTHER' or wf_to_use == 'COMPOUND_HET' or wf_to_use == 'UNIPARENTAL_ISODISOMY':
			
			workflow_df = workflow_df[((workflow_df['gnomADg_AF_POPMAX'] <0.01) | (pd.isna(workflow_df['gnomADg_AF_POPMAX']) )) &
				   ((workflow_df['gnomADe_AF_POPMAX'] <0.01 ) | (pd.isna(workflow_df['gnomADe_AF_POPMAX'])))
				   ]
			
		elif wf_to_use == 'RECCESSIVE_SEX':
			
			workflow_df = workflow_df[((workflow_df['gnomADg_AF_POPMAX'] <0.01) | (pd.isna(workflow_df['gnomADg_AF_POPMAX']) )) &
				   ((workflow_df['gnomADe_AF_POPMAX'] <0.01 ) | (pd.isna(workflow_df['gnomADe_AF_POPMAX']))) &
				   (workflow_df['AC'] < 8)
				   ]
			
		elif wf_to_use == 'RECESSIVE_AUTOSOMAL':
			
			  workflow_df = workflow_df[((workflow_df['gnomADg_AF_POPMAX'] <0.01) | (pd.isna(workflow_df['gnomADg_AF_POPMAX']) )) &
				   ((workflow_df['gnomADe_AF_POPMAX'] <0.01 ) | (pd.isna(workflow_df['gnomADe_AF_POPMAX']))) &
				   (workflow_df['AC'] < 4)
				   ]      
		
		elif wf_to_use == 'DOMINANT_SEX':
			
			  workflow_df = workflow_df[((workflow_df['gnomADg_AF_POPMAX'] <0.0075) | (pd.isna(workflow_df['gnomADg_AF_POPMAX']) )) &
				   ((workflow_df['gnomADe_AF_POPMAX'] <0.001 ) | (pd.isna(workflow_df['gnomADe_AF_POPMAX']))) &
				   (workflow_df['AC'] < 6)
				   ]
		elif wf_to_use == 'DOMINANT_AUTOSOMAL' or wf_to_use == 'DENOVO_HC' or wf_to_use == 'DENOVO_LC':
			
			  workflow_df = workflow_df[((workflow_df['gnomADg_AF_POPMAX'] <0.0075) | (pd.isna(workflow_df['gnomADg_AF_POPMAX']) )) &
				   ((workflow_df['gnomADe_AF_POPMAX'] <0.001 ) | (pd.isna(workflow_df['gnomADe_AF_POPMAX']))) &
				   (workflow_df['AC'] < 4)
				   ]
				
		else:
			
			raise Exception('Unknown workflow')
			
		master_sample_df = pd.concat([master_sample_df, workflow_df])

	print (f"{sample} INFO: Variants remaining after filtering: {master_sample_df.groupby('variant').count().shape[0]}")

	# If we want to add the constrained coding regions then do that.
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

		master_sample_df['disease'] = master_sample_df.apply(apply_panel_app_data, axis=1,args =(panel_app_dict,panel_app_dump_max_time,))


	print (master_sample_df[['variant', 'SYMBOL']].head())


if use_local_panel_app_dump:

	panel_app_df = pd.DataFrame(panel_app_dict)

	panel_app_df.transpose().to_csv(local_panel_app_dump, sep='\t')

















			   