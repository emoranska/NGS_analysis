import pandas as pd

# reading csv files as pandas DataFrame
# data1 = pd.read_csv('P1_no_Un_NW_MT.vcf.default.csv', sep='\t')
# data2 = pd.read_csv('P1_no_Un_NW_MT.vcf.sample.GT.csv', sep='\t')
# print(data1, data2)


def read_csv(csv):
    data = pd.read_csv(csv, sep='\t')
    return data


# reset index in the second DataFrame (data2_reset) which columns will be added to data1
# concatenate data1 and data2_reset - add data2_reset on the right of data1
def concat_samples_gt(csv_1, csv_2):
    data1 = read_csv(csv_1)
    data2_reset = read_csv(csv_2).reset_index(drop=True)
    just_gt = pd.concat([data1, data2_reset], axis=1)
    return just_gt

# displaying result
# print(just_gt.head())


# remove records with the same genotypes for all samples
def remove_the_same_gt(csv_1, csv_2):
    data = concat_samples_gt(csv_1, csv_2)
    i = data.apply(lambda r: len(set(r[9:])) != 1, axis=1)
    no_the_same_variant = data[i]
    return no_the_same_variant


# save results to csv file
# no_the_same_variant_csv = no_the_same_variant.to_csv('P1_no_the_same_variant.csv', index=False)
def to_csv(pd_data, filename):
    csv_out = pd_data.to_csv(filename, index=False)
    return csv_out


# print(no_the_same_variant.head())
