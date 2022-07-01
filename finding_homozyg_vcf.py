from fuc import pyvcf
import numpy as np

vf = pyvcf.VcfFrame.from_file('chr4_testdata1.vcf')
just_gt = vf.strip('GT')
i = just_gt.df.apply(lambda row: len(set(row[9:])) != 1, axis=1)
no_the_same_variant = just_gt.df[i]

no_the_same_variant.astype("string")

# replace missing values signed as '.' with 'NaN'
ntsv_nan = no_the_same_variant.replace('.', np.nan)

print(ntsv_nan.to_string())

# remove records with >=6 samples with missing values ('.' replaced with NaN)
drop_dots = ntsv_nan.dropna(thresh=7, subset=['P1-25', 'P1-93', 'P1-88', 'P1-6', 'P1-89', 'P1-26', 'P1-12', 'P1-92',
                                              'P1-22', 'P1-90', 'P1-28', 'P1-95'])

# choosing rows that are >= 3 1/1 AND >=3 0/0 (solution from CAPSLOCK)
cond_11x3_true = (drop_dots.iloc[:, 8:] == "1/1").sum(axis=1) >= 3
cond_00x3_true = (drop_dots.iloc[:, 8:] == "0/0").sum(axis=1) >= 3
result = drop_dots[cond_11x3_true & cond_00x3_true]

print(result.to_string())

# replace missing values signed as 'NaN' with '.' for saving df to VcfFrame and vcf file
final_df = result.replace(np.nan, '.')

print(final_df.to_string())

# create new VcfFrame
to_vcf = pyvcf.VcfFrame(['##fileformat=VCFv4.3'], final_df)

# save VcfFrame to vcf file
chr4_test_homo_regions = to_vcf.to_file('1_chr4_test_homo_regions.vcf')
