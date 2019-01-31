from utils.utils import *
from utils.inheritance_utils import *
import unittest


class InteritancePatternTest(unittest.TestCase):


	def test_is_compound_het(self):

		ch_dict = {'1': 1, '2': 1, '3': 0, '4':3, '5': 2} 

		self.assertEqual(is_compound_het('1', ch_dict), False)
		self.assertEqual(is_compound_het('2', ch_dict), False)
		self.assertEqual(is_compound_het('3', ch_dict), False)
		self.assertEqual(is_compound_het('4', ch_dict), True)
		self.assertEqual(is_compound_het('5', ch_dict), True)

	def test_is_dominant_autosomal(self):

		# Typical dominant variant
		variant_chrom = '1'
		variant_ref = 'A'
		variant_alt = 'G'
		variant_gt = ['A', 'G']
		variant_sex = 'Male'

		self.assertEqual(is_dominant_autosomal(variant_chrom, variant_ref, variant_alt, variant_gt, variant_sex), True)

		# Typical dominant variant
		variant_chrom = '1'
		variant_ref = 'A'
		variant_alt = 'G'
		variant_gt = ['A', 'G']
		variant_sex = 'Female'

		self.assertEqual(is_dominant_autosomal(variant_chrom, variant_ref, variant_alt, variant_gt, variant_sex), True)

		# Typical dominant variant
		variant_chrom = '1'
		variant_ref = 'A'
		variant_alt = 'G'
		variant_gt = ['G', 'A']
		variant_sex = 'Male'

		self.assertEqual(is_dominant_autosomal(variant_chrom, variant_ref, variant_alt, variant_gt, variant_sex), True)

		# Typical reccessive variant
		variant_chrom = '1'
		variant_ref = 'A'
		variant_alt = 'G'
		variant_gt = ['G', 'G']
		variant_sex = 'Male'

		self.assertEqual(is_dominant_autosomal(variant_chrom, variant_ref, variant_alt, variant_gt, variant_sex), False)

		# On X chrom on a male
		variant_chrom = 'X'
		variant_ref = 'A'
		variant_alt = 'G'
		variant_gt = ['A', 'G']
		variant_sex = 'Male'

		self.assertEqual(is_dominant_autosomal(variant_chrom, variant_ref, variant_alt, variant_gt, variant_sex), False)

		# On X chrom on a female
		variant_chrom = 'X'
		variant_ref = 'A'
		variant_alt = 'G'
		variant_gt = ['A', 'G']
		variant_sex = 'Female'

		self.assertEqual(is_dominant_autosomal(variant_chrom, variant_ref, variant_alt, variant_gt, variant_sex), True)


		# Reccessive X chrom on a female
		variant_chrom = 'X'
		variant_ref = 'A'
		variant_alt = 'G'
		variant_gt = ['G', 'G']
		variant_sex = 'Female'

		self.assertEqual(is_dominant_autosomal(variant_chrom, variant_ref, variant_alt, variant_gt, variant_sex), False)

		# HET on X chrom on a male
		variant_chrom = 'X'
		variant_ref = 'A'
		variant_alt = 'G'
		variant_gt = ['A', 'G']
		variant_sex = 'Male'

		self.assertEqual(is_dominant_autosomal(variant_chrom, variant_ref, variant_alt, variant_gt, variant_sex), False)


		# HET on Y chrom on a female
		variant_chrom = 'Y'
		variant_ref = 'A'
		variant_alt = 'G'
		variant_gt = ['A', 'G']
		variant_sex = 'Female'

		self.assertEqual(is_dominant_autosomal(variant_chrom, variant_ref, variant_alt, variant_gt, variant_sex), False)















if __name__ == '__main__':
    unittest.main()