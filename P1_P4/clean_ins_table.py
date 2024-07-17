import pandas as pd

# cleaning .csv file from duplicated MITE family names for visual intend
df = pd.read_csv('../files/P1_MITEs_ins_all_count.csv', sep='\t')

s1 = df.duplicated(subset=['family'])
s2 = df.duplicated(subset=['family', 'ins'])

df['family'] = df['family'].where(~s1, '')
df['ins'] = df['ins'].where(~s2, '')

print(df)
df.to_csv('../files/P1_MITEs_ins_all_count_cleaned.csv', sep='\t', index=False)
