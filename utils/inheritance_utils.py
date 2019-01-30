def is_compound_het(transcript, compound_het_dict):
		"""
		Takes a transcript and returns True if the transcript has more than one variant in.
		
		"""
		
		count = compound_het_dict[transcript]
		
		if count > 1:
				
				return True
		
		else:
				
				return False

def is_dominant_autosomal(chrom, ref, alt, sample_gt, sample_sex):
		"""
		For a single sample i.e. not trio returns True if the variant genotype is het.
		
		For example the gentype if ref/alt or alt/ref or 
		
		
		"""
		
		if (chrom == 'X' or chrom == 'Y') and sample_sex == 'Male':
				
				return False
		
		else:
						
				if ref in sample_gt and alt in sample_gt:
						
						return True
				
				else:
						
						return False
				
def is_dominant_sex(chrom, ref, alt, sample_gt, sample_sex):
		"""
		If the variant is in a male and on the x chromsome return true

		"""
		
		if (chrom == 'X' or chrom == 'Y') and sample_sex == 'Male':
				
				return True
		
		else:
				
				return False
		
def is_recessive_autosomal(chrom, ref, alt, sample_gt, sample_sex):
		"""
		Does the genotype have two alt allelles
		
		false if male and on x chrom
		
		"""
		
		if (chrom == 'X' or chrom == 'Y') and sample_sex == 'Male':
				
				return False
		
		else:
						
				if sample_gt[0] == alt and sample_gt[1] == alt:
						
						return True
				
				else:
						
						return False
				
def is_recessive_sex(chrom, ref, alt, sample_gt, sample_sex):
		"""
		Does the genotype have two alt allelles and is on x chromsome and female
		
		"""
		
		if chrom == 'X' and sample_sex == 'Female':
						
				if sample_gt[0] == alt and sample_gt[1] == alt:
						
						return True
				
				else:
						
						return False
				
		else:
				
				return False
		
def is_high_confidence_de_novo(chrom,
															 ref,
															 alt,
															 sample_gt,
															 mother_gt,
															 father_gt,
															mother_dp,
															father_dp,
															mother_gq,
															father_gq,
															min_parental_depth,
															min_parental_gq):
				
				
		if (alt in sample_gt and alt not in mother_gt and alt not in father_gt) and (mother_dp >= min_parental_depth and father_dp >= min_parental_depth and mother_gq >=min_parental_gq and father_gq >=min_parental_gq):
				
				return True
		
		else:
				
				return False
		
		
		
def is_low_confidence_de_novo(chrom,
															 ref,
															 alt,
															 sample_gt,
															 mother_gt,
															 father_gt,
															mother_dp,
															father_dp):

				
				
		if (alt in sample_gt and alt not in mother_gt and alt not in father_gt) and (mother_dp >= 1 and father_dp >= 1):
				
				return True
		
		else:
				
				return False
		
		
def is_uniparental_isodisomy(chrom, ref, alt, sample_gt, mother_gt,
														 father_gt,
														 sample_sex,
														 min_parental_depth,
														 min_parental_gq,
														 father_dp,
														 mother_dp,
														 father_gq,
														 mother_gq):
				
		if chrom != 'X':
				
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


def annotate_workflow_trio(df,sample, mother, father, sample_sex, compound_het_dict, min_parental_depth_dn, min_parental_gq_dn, min_parental_depth_uid, min_parental_gq_uid ):
		
		workflows = []
		
		chrom = df['CHROM']
		ref = df['REF']
		alt = df['ALT']
		sample_gt = df['sample_' + sample + '_GT']
		mother_gt = df['sample_' + mother + '_GT']
		father_gt = df['sample_' + father + '_GT']
		mother_dp = df['sample_' + mother + '_DP']
		father_dp = df['sample_' + father + '_DP']
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
		
		if is_compound_het(transcript, compound_het_dict) == True:
				
				workflows.append('COMPOUND_HET')
				
		if is_dominant_autosomal(chrom, ref, alt, sample_gt, sample_sex) == True:
				
				workflows.append('DOMINANT_AUTOSOMAL')
				
		if is_dominant_sex(chrom, ref, alt, sample_gt, sample_sex) == True:
				
				workflows.append('DOMINANT_SEX')
				
		if is_recessive_autosomal(chrom, ref, alt, sample_gt, sample_sex) == True:
				
				workflows.append('RECESSIVE_AUTOSOMAL')
				
		if is_recessive_sex(chrom, ref, alt, sample_gt, sample_sex) == True:
				
				workflows.append('RECCESSIVE_SEX')
				
		if is_high_confidence_de_novo(chrom, ref, alt, sample_gt, mother_gt, father_gt,mother_dp,father_dp,mother_gq,father_gq, min_parental_depth_dn, min_parental_gq_dn) == True:
				
				workflows.append('DENOVO_HC')
				
		if is_low_confidence_de_novo(chrom, ref, alt, sample_gt,mother_gt, father_gt, mother_dp, father_dp) == True:
				
				workflows.append('DENOVO_LC')
				
		if is_uniparental_isodisomy(chrom, ref,alt, sample_gt, mother_gt, father_gt, sample_sex, min_parental_depth_uid, min_parental_gq_uid, father_dp, mother_dp, father_gq, mother_gq):
				
				workflows.append('UNIPARENTAL_ISODISOMY')
				
		if len(workflows) == 0:
				
				workflows.append('OTHER')

		if 'DENOVO_HC' in workflows and 'DENOVO_LC' in workflows:

			workflows.remove('DENOVO_LC')
				
		return '|'.join(workflows)

def annotate_workflow_single(df,sample, sample_sex, compound_het_dict):
		
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
		
		if is_compound_het(transcript, compound_het_dict) == True:
				
				workflows.append('COMPOUND_HET')
				
		if is_dominant_autosomal(chrom, ref, alt, sample_gt, sample_sex) == True:
				
				workflows.append('DOMINANT_AUTOSOMAL')
				
		if is_dominant_sex(chrom, ref, alt, sample_gt, sample_sex) == True:
				
				workflows.append('DOMINANT_SEX')
				
		if is_recessive_autosomal(chrom, ref, alt, sample_gt, sample_sex) == True:
				
				workflows.append('RECESSIVE_AUTOSOMAL')
				
		if is_recessive_sex(chrom, ref, alt, sample_gt, sample_sex) == True:
				
				workflows.append('RECCESSIVE_SEX')
				
		if len(workflows) == 0:
				
				workflows.append('OTHER')
				
		return '|'.join(workflows)