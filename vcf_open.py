from fuc import pyvcf

# open vcf file with fuc.pyvcf
vf = pyvcf.VcfFrame.from_file('P1_no_Un_NW_MT.vcf')

# remove useless columns, leave GT
just_gt = vf.strip('GT')

# remove records with the same genotypes for all samples
i = just_gt.df.apply(lambda r: len(set(r[9:])) != 1, axis=1)
no_the_same_variant = just_gt.df[i]

# create new VcfFrame
to_vcf = pyvcf.VcfFrame([''], no_the_same_variant)

# save VcfFrame to vcf file
P1_no_the_same_variant = to_vcf.to_file('P1_no_the_same_variant.vcf')

print(type(no_the_same_variant))



