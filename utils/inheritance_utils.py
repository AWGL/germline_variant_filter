"""
Functions for working out the inheritiance of variants.

"""

def is_compound_het(chrom, ref, alt, sample_gt, sample_sex,transcript, compound_het_dict):
	"""
	Takes a transcript and returns True if the transcript has more than one variant in.

	"""

	count = compound_het_dict[transcript]
	
	if (chrom not in ['X', 'Y', 'MT'] or chrom =='X' and sample_sex =='Female') and count >1 and sample_gt.count(alt) ==1:

		return True

	else:

		return False


def is_dominant_autosomal(chrom, ref, alt, sample_gt):
	"""
	Returns whether the variant matches a dominant/HET pattern and is not on a autosome

	For example the gentype if ref/alt or alt/ref

	Note: Does not look at parents to check for mendellian violation.


	"""

	if sample_gt.count(alt) ==1 and chrom not in ['X', 'Y', 'MT']:

		return True

	else:

		return False


def is_dominant_x_female(chrom, ref, alt,sample_gt, sample_sex):
	"""
	Does the variant match a dominant/HET pattern and on the X chromosome in a female sample?
	"""
	
	if sample_gt.count(alt) ==1 and chrom == 'X' and sample_sex == 'Female':
		
		return True
	
	else:
		
		return False  


def is_x_linked_male(chrom, ref, alt,sample_gt, sample_sex):
	"""
	Is the variant on the X chromsome and the sample male
	"""
	
	if alt in sample_gt and chrom == 'X' and sample_sex == 'Male':
		
		return True
	
	else:
		
		return False  


def is_y_linked_male(chrom, ref, alt,sample_gt, sample_sex):
	"""
	Is the variant on the Y chromsome and the sample male?
	"""
	
	if alt in sample_gt and chrom == 'Y' and sample_sex == 'Male':
		
		return True
	
	else:
		
		return False


def is_reccessive_autosomal(chrom, ref, alt, sample_gt):
	"""
	Is the variant homozygous alt and on an autosome?
	
	"""
	
	if sample_gt.count(alt) == 2 and chrom not in ['X', 'Y', 'MT']:
		
		return True
	
	else:
		
		return False


def is_reccessive_x_female(chrom, ref, alt, sample_gt, sample_sex):
	"""
	Is the variant hom alt in a female on the x chromsome
  
	"""
	
	if sample_gt.count(alt) == 2 and chrom == 'X' and sample_sex == 'Female':
		
		return True
	
	else:
		
		return False


def is_mitochrondial(chrom, ref, alt, sample_gt):
	"""
	Is the variant on the mito chromsome?
	
	"""
	
	if alt in sample_gt and chrom == 'MT':
		
		return True
	
	else:
		
		return False


def is_high_confidence_de_novo(chrom, ref, alt, sample_gt, mother_gt, father_gt, mother_dp, father_dp, mother_gq, father_gq, min_parental_depth, min_parental_gq):
	"""
	Returns True if the variant is in the proband and not in the mother and father.

	Parent gentypes must also match min_parental_depth and min_parental_gq requirements.

	genotype arguments such as sample_gt are a list containing both genotypes e.g ['A', 'G']

	"""

	if (alt in sample_gt and alt not in mother_gt and alt not in father_gt) and (mother_dp >= min_parental_depth and father_dp >= min_parental_depth and mother_gq >= min_parental_gq and father_gq >= min_parental_gq):
				
		return True
		
	else:
				
		return False
		
		
		
def is_low_confidence_de_novo(chrom, ref, alt, sample_gt, mother_gt, father_gt):
	"""
	Returns True if the sample could possibly be a de novo.

	Does not take into account parental genotype quality.

	A call in the proband but missing calls in parents will return True

	"""

	if (alt in sample_gt and alt not in mother_gt and alt not in father_gt):
			
		return True
	
	else:
			
		return False
		
		
def is_uniparental_isodisomy(chrom, ref, alt, sample_gt, mother_gt, father_gt, sample_sex, min_parental_depth, min_parental_gq, father_dp, mother_dp, father_gq, mother_gq):
	"""
	Returns True if the variant matches a uniparental isodisomy pattern of inheritance.

	"""
			
	if chrom not in ['X', 'Y', 'MT']:
			
		if sample_gt.count(alt) == 2 and (mother_gt.count(alt) ==1 and father_gt.count(alt) == 0) and mother_dp >= min_parental_depth and father_dp >= min_parental_depth and father_gq >= min_parental_gq and mother_gq >= min_parental_gq:
					
			return True
			
		if sample_gt.count(alt) == 2 and (father_gt.count(alt) ==1 and mother_gt.count(alt) == 0) and mother_dp >= min_parental_depth and father_dp >= min_parental_depth and father_gq >= min_parental_gq and mother_gq >= min_parental_gq:
					
			return True
			
	elif chrom == 'X' and sample_sex == 'Female':
			
		if sample_gt.count(alt) == 2 and (mother_gt.count(alt) ==1 and father_gt.count(alt) == 0) and mother_dp >= min_parental_depth and father_dp >= min_parental_depth and father_gq >= min_parental_gq and mother_gq >= min_parental_gq:
					
			return True
			
		if sample_gt.count(alt) == 2 and (mother_gt.count(alt) == 0 and father_gt.count(alt) ==2) and mother_dp >= min_parental_depth and father_dp >= min_parental_depth and father_gq >= min_parental_gq and mother_gq >= min_parental_gq:
					
			return True
			
	return False


def annotate_workflow_trio(df, sample, mother, father, sample_sex, compound_het_dict, min_parental_depth_dn, min_parental_gq_dn, min_parental_depth_uid, min_parental_gq_uid, gt_depth_tag):
	"""
	For samples which are a trio add the workflow annotation to the dataframe.

	"""

	workflows = []

	chrom = df['CHROM']
	ref = df['REF']
	alt = df['ALT']
	sample_gt = df['sample_' + sample + '_GT']
	mother_gt = df['sample_' + mother + '_GT']
	father_gt = df['sample_' + father + '_GT']
	mother_dp = df['sample_' + mother + '_' + gt_depth_tag]
	father_dp = df['sample_' + father + '_' + gt_depth_tag]
	mother_gq = df['sample_' + mother + '_GQ']
	father_gq = df['sample_' + father + '_GQ']
	transcript = df['Feature']

	if '|' in sample_gt:

		sample_gt = sample_gt.split('|')

	elif '/' in sample_gt:

		sample_gt = sample_gt.split('/')

	else:

		raise Exception('There is a sample genotype without either | or / in.')


	if '|' in mother_gt:

		mother_gt = mother_gt.split('|')

	elif '/' in mother_gt:

		mother_gt = mother_gt.split('/')

	else:

		raise Exception('There is a maternal genotype without either | or / in.')


	if '|' in father_gt:

		father_gt = father_gt.split('|')

	elif '/' in father_gt:

		father_gt = father_gt.split('/')

	else:

		raise Exception('There is a paternal genotype without either | or / in.')

	if is_compound_het(chrom, ref, alt, sample_gt, sample_sex, transcript, compound_het_dict) == True:

		workflows.append('COMPOUND_HET')

	if is_dominant_autosomal(chrom, ref, alt, sample_gt) == True:

		workflows.append('DOMINANT_AUTOSOMAL')

	if is_dominant_x_female(chrom, ref, alt, sample_gt, sample_sex) == True:

		workflows.append('DOMINANT_X_FEMALE')

	if is_x_linked_male(chrom, ref, alt, sample_gt, sample_sex) == True:

		workflows.append('X_LINKED_MALE')
		
	if is_y_linked_male(chrom, ref, alt, sample_gt, sample_sex) == True:
		
		workflows.append('Y_LINKED_MALE')

	if is_reccessive_autosomal(chrom, ref, alt, sample_gt) == True:

		workflows.append('RECCESSIVE_AUTOSOMAL')
		
	if is_reccessive_x_female(chrom, ref, alt, sample_gt, sample_sex) == True:
		
		workflows.append('RECCESSIVE_X_FEMALE')
		
	if is_mitochrondial(chrom, ref, alt, sample_gt) ==True:
		
		workflows.append('MITOCHONDRIAL')

	if is_high_confidence_de_novo(chrom, ref, alt, sample_gt, mother_gt, father_gt,mother_dp,father_dp,mother_gq,father_gq, min_parental_depth_dn, min_parental_gq_dn) == True:

		workflows.append('DENOVO_HC')

	elif is_low_confidence_de_novo(chrom, ref, alt, sample_gt, mother_gt, father_gt) == True:

		workflows.append('DENOVO_LC')

	if is_uniparental_isodisomy(chrom, ref, alt, sample_gt, mother_gt, father_gt, sample_sex, min_parental_depth_uid, min_parental_gq_uid, father_dp, mother_dp, father_gq, mother_gq) == True:

		workflows.append('UNIPARENTAL_ISODISOMY')


	# If the variant does not fit in any then add other 
	if len(workflows) == 0 or (len(workflows) ==1 and ('DENOVO_HC' in workflows or 'DENOVO_LC' in workflows)): 

		workflows.append('OTHER')

	return '|'.join(workflows)


def annotate_workflow_single(df,sample, sample_sex, compound_het_dict):
	"""
	Annotate workflow for single samples not in a trio.

	"""

	workflows = []

	chrom = df['CHROM']
	ref = df['REF']
	alt = df['ALT']
	sample_gt = df['sample_' + sample + '_GT']
	transcript = df['Feature']


	if '|' in sample_gt:

		sample_gt = sample_gt.split('|')

	elif '/' in sample_gt:

		sample_gt = sample_gt.split('/')

	else:

		raise Exception('There is a sample genotype without either | or / in.')

	if is_compound_het(chrom, ref, alt, sample_gt, sample_sex, transcript, compound_het_dict) == True:

		workflows.append('COMPOUND_HET')

	if is_dominant_autosomal(chrom, ref, alt, sample_gt) == True:

		workflows.append('DOMINANT_AUTOSOMAL')

	if is_dominant_x_female(chrom, ref, alt, sample_gt, sample_sex) == True:

		workflows.append('DOMINANT_X_FEMALE')

	if is_x_linked_male(chrom, ref, alt, sample_gt, sample_sex) == True:

		workflows.append('X_LINKED_MALE')
		
	if is_y_linked_male(chrom, ref, alt, sample_gt, sample_sex) == True:
		
		workflows.append('Y_LINKED_MALE')

	if is_reccessive_autosomal(chrom, ref, alt, sample_gt) == True:

		workflows.append('RECCESSIVE_AUTOSOMAL')
		
	if is_reccessive_x_female(chrom, ref, alt, sample_gt, sample_sex) == True:
		
		workflows.append('RECCESSIVE_X_FEMALE')
		
	if is_mitochrondial(chrom, ref, alt, sample_gt) == True:
		
		workflows.append('MITOCHONDRIAL')

	if len(workflows) == 0:

		workflows.append('OTHER')

	return '|'.join(workflows)

