from fuc import pyvcf
import pandas as pd
import time

start_time = time.time()

# open vcf file as pandas df
# vf1 = pyvcf.VcfFrame.from_file('P1_chr1_3_homozyg_reg.vcf')
vf1 = pyvcf.VcfFrame.from_file('P4_chr1_3_homozyg_reg.vcf')
df1 = vf1.df

# vf2 = pyvcf.VcfFrame.from_file('P1_chr4_6_homozyg_reg.vcf')
vf2 = pyvcf.VcfFrame.from_file('P4_chr4_6_homozyg_reg.vcf')
df2 = vf2.df

# vf3 = pyvcf.VcfFrame.from_file('P1_chr7_9_homozyg_reg.vcf')
vf3 = pyvcf.VcfFrame.from_file('P4_chr7_9_homozyg_reg.vcf')
df3 = vf3.df

homozyg_reg = pd.concat([df1, df2, df3], ignore_index=True)
print(homozyg_reg.to_string(max_rows=50))

to_vcf = pyvcf.VcfFrame(['##fileformat=VCFv4.3'], homozyg_reg)
homo_regions = to_vcf.to_file('P4_EL10_chr1_9_homozyg_reg.vcf')
