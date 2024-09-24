import pandas as pd

genes = pd.read_csv('../files/genes.csv', sep=',').rename(columns={'Gene': 'GeneID'})
genes['gene_length_kb'] = (genes['End'] - genes['Start']) / 1000
print(genes.to_string(max_rows=30))

raw_counts = pd.read_csv('../files/kallisto_P1.csv', sep='\t')
print(raw_counts.to_string(max_rows=30))


def tpm_calc():
    genes_and_raw = genes.merge(raw_counts, how='outer')
    print(genes_and_raw.to_string(max_rows=30))

    # calculate RPK values
    cols = genes_and_raw.filter(like='P1-').columns
    # cols = genes_and_raw.filter(like='P4-').columns

    genes_and_raw[cols] = genes_and_raw[cols].div(genes_and_raw['gene_length_kb'], axis=0)
    print(genes_and_raw.to_string(max_rows=30))

    rpk = genes_and_raw[cols]

    # calculate scaling factor
    scaling_factor = genes_and_raw[cols].sum() / 1000000
    scaling_factor = scaling_factor.reset_index()
    scaling_factor.columns = ['sample', 'factor']
    print(scaling_factor.to_string(max_rows=30))

    tpm = rpk / scaling_factor.set_index('sample')['factor']
    tpm['GeneID'] = genes_and_raw['GeneID']
    print(tpm)

    tpm_sum = tpm.sum()
    print(tpm_sum)

    # tpm.to_csv('../files/P4_kallisto_tpm.csv', sep='\t', index=False)
    return tpm


# tpm_calc()
