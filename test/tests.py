from utils.utils import *
from utils.inheritance_utils import *
import unittest
import pandas as pd

class WorkFlowTrioTest(unittest.TestCase):
	"""
	Test the workflow assignment for trios

	"""


	def test_trio_male(self):

		df = pd.read_csv('test/test_data_male_trio.csv', sep='\t')

		compound_het_dict = df.groupby('Feature').count()['CHROM'].to_dict()

		sample = 'proband'
		mother = 'mother'
		father = 'father'
		sample_sex = 'Male'
		min_parental_depth_dn =10
		min_parental_gq_dn = 20
		min_parental_depth_uid = 10
		min_parental_gq_uid =10

		df['Workflow'] = df.apply(annotate_workflow_trio, axis=1, args=(sample, mother, father, sample_sex, compound_het_dict, min_parental_depth_dn, min_parental_gq_dn, min_parental_depth_uid, min_parental_gq_uid, 'DP'))

		for row in df.itertuples():

			self.assertEqual(row.Workflow, row.Comment)
	
		
	def test_trio_female(self):

		df = pd.read_csv('test/test_data_female_trio.csv', sep='\t')

		compound_het_dict = df.groupby('Feature').count()['CHROM'].to_dict()

		sample = 'proband'
		mother = 'mother'
		father = 'father'
		sample_sex = 'Female'
		min_parental_depth_dn =10
		min_parental_gq_dn = 20
		min_parental_depth_uid = 10
		min_parental_gq_uid =10

		df['Workflow'] = df.apply(annotate_workflow_trio, axis=1, args=(sample, mother, father, sample_sex, compound_het_dict, min_parental_depth_dn, min_parental_gq_dn, min_parental_depth_uid, min_parental_gq_uid, 'DP'))

		for row in df.itertuples():

			self.assertEqual(row.Workflow, row.Comment)

class WorkFlowSingleTest(unittest.TestCase):

	"""
	Test the workflow assignment for single samples.

	"""

	def test_single_male(self):

		df = pd.read_csv('test/test_data_male_single.csv', sep='\t')

		compound_het_dict = df.groupby('Feature').count()['CHROM'].to_dict()

		sample = 'proband'
		mother = 'mother'
		father = 'father'
		sample_sex = 'Male'
		min_parental_depth_dn =10
		min_parental_gq_dn = 20
		min_parental_depth_uid = 10
		min_parental_gq_uid =10

		df['Workflow'] = df.apply(annotate_workflow_single, axis=1, args=(sample, sample_sex, compound_het_dict))

		for row in df.itertuples():

			self.assertEqual(row.Workflow, row.Comment)


	def test_single_female(self):

		df = pd.read_csv('test/test_data_female_single.csv', sep='\t')

		compound_het_dict = df.groupby('Feature').count()['CHROM'].to_dict()

		sample = 'proband'
		mother = 'mother'
		father = 'father'
		sample_sex = 'Female'
		min_parental_depth_dn =10
		min_parental_gq_dn = 20
		min_parental_depth_uid = 10
		min_parental_gq_uid =10

		df['Workflow'] = df.apply(annotate_workflow_single, axis=1, args=(sample, sample_sex, compound_het_dict))

		for row in df.itertuples():

			self.assertEqual(row.Workflow, row.Comment)





if __name__ == '__main__':
    unittest.main()

