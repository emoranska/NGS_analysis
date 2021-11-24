import pandas as pd

# reading csv files as pandas DataFrame
data1 = pd.read_csv('P1_no_Un_NW_MT.vcf.default.csv', sep='\t')
data2 = pd.read_csv('P1_no_Un_NW_MT.vcf.sample.GT.csv', sep='\t')
#print(data1, data2)

# reset index in the second DataFrame (data2) which columns will be added to data1
data2 = data2.reset_index(drop=True)
#print(data2)

# concatenate data1 and data2 - add data2 on the right of data1
just_gt = pd.concat([data1, data2], axis=1)

# displaying result
#print(just_gt.head())

# remove records with the same genotypes for all samples
i = just_gt.apply(lambda r: len(set(r[9:])) != 1, axis=1)
no_the_same_variant = just_gt[i]

# save results to csv file
# no_the_same_variant_csv = no_the_same_variant.to_csv('P1_no_the_same_variant.csv', index=False)

print(no_the_same_variant.head())
