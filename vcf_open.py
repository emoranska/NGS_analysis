from fuc import pyvcf
import numpy as np

import pandas as pd

vf = pyvcf.VcfFrame.from_file('P1_short_test.vcf')

just_gt = vf.strip('GT')

i = just_gt.df.apply(lambda r: len(set(r[9:])) != 1, axis=1)

no_the_same_variant = just_gt.df[i]

to_vcf = pyvcf.VcfFrame([''], no_the_same_variant)

P1_no_the_same_variant = to_vcf.to_file('P1_no_the_same_variant')

print(just_gt.df, no_the_same_variant, type(just_gt), type(no_the_same_variant))



