from fuc import pyvcf
import numpy as np

vcf_input = input('Enter the input vcf file name: ')

vf = pyvcf.VcfFrame.from_file(vcf_input)
just_gt = vf.strip('GT')
i = just_gt.df.apply(lambda row: len(set(row[9:])) != 1, axis=1)
no_the_same_variant = just_gt.df[i]

no_the_same_variant.astype("string")

# replace missing values signed as '.' with 'NaN'
ntsv_nan = no_the_same_variant.replace('.', np.nan)

# remove records with >=6 samples with missing values ('.' replaced with NaN)
drop_dots = ntsv_nan.dropna(thresh=7, subset=['P1-25', 'P1-93', 'P1-88', 'P1-6', 'P1-89', 'P1-26', 'P1-12', 'P1-92',
                                              'P1-22', 'P1-90', 'P1-28', 'P1-95'])

# choosing rows that are >= 3 1/1 AND >=3 0/0 (solution from CAPSLOCK)
max_rows = 2

cond_11x3_true = (drop_dots.iloc[:, 8:] == "1/1").sum(axis=1) >= 3
cond_00x3_true = (drop_dots.iloc[:, 8:] == "0/0").sum(axis=1) >= 3

# both triplets
cond_both_true = (cond_11x3_true & cond_00x3_true)
# identify consecutive false
cond_max_2_rows = (cond_11x3_true & cond_00x3_true).cumsum()
# identify where consecutive false is > max_rows
cond_max_2_rows = cond_max_2_rows.groupby(cond_max_2_rows).transform('size') <= max_rows+1

"""
# at least one triplet
cond_1_true_and_max_2_gaps = (cond_11x3_true & ~cond_00x3_true) | (~cond_11x3_true & cond_00x3_true)
# pick only if consecutive false <= max_rows and at least one triplet
cond_1_true_and_max_2_gaps = (cond_1_true_and_max_2_gaps & cond_max_2_rows)
# filter out where neither condition is respected
result = drop_dots[cond_both_true|cond_1_true_and_max_2_gaps]
"""

result = drop_dots[cond_both_true | cond_max_2_rows]

print(result.to_string())

# replace missing values signed as 'NaN' with '.' for saving df to VcfFrame and vcf file
final_df = result.replace(np.nan, '.')

# create new VcfFrame
to_vcf = pyvcf.VcfFrame(['##fileformat=VCFv4.3'], final_df)

# save VcfFrame to vcf file
vcf_output = input('Enter the output vcf file name: ')
chr4_test_homo_regions = to_vcf.to_file(vcf_output)
