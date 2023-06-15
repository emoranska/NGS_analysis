from fuc import pyvcf
import pandas as pd
import time

start_time = time.time()

# vf = pyvcf.VcfFrame.from_file('P1_chr1_100_homozyg_reg.vcf')
vf = pyvcf.VcfFrame.from_file('P1_chr1_3_homozyg_reg.vcf')
print(type(vf))
df = vf.df
print(type(df))

# print(df.iloc[0,11], type(df.iloc[0, 11]))

# df['1/1'] = df.apply(lambda x: x[9:] == '1/1', axis=1)
# i = df.apply(lambda x: 1 if df.loc[x, 'P1-25'] == '1/1' else 0, axis=1)
# df['1/1'] = df[i]

# df['1/1'] = ['P1-25' if x == '1/1' else 0 for x in df['P1-25']]
# print(df)

df['pos_back'] = df['POS'].shift(1)
df['pos_next'] = df['POS'].shift(-1)

cols = df.filter(like='P1-').columns
sep = ','
df['one_one'] = pd.Series(df[cols].eq('1/1').dot((cols + sep))).str.rstrip(sep)
df['zero_zero'] = pd.Series(df[cols].eq('0/0').dot((cols + sep))).str.rstrip(sep)

# print(df.to_string())
#

df['one_one_back'] = df['one_one'].shift(1).str.split(',')
df['one_one_next'] = df['one_one'].shift(-1).str.split(',')
df['zero_zero_back'] = df['zero_zero'].shift(1).str.split(',')
df['zero_zero_next'] = df['zero_zero'].shift(-1).str.split(',')

df['one_one'] = df['one_one'].str.split(',')
df['zero_zero'] = df['zero_zero'].str.split(',')

# print(df.loc[0, 'one_one'], type(df.loc[0, 'one_one']))
print(df.loc[1, 'one_one_back'], type(df.loc[0, 'one_one_back']))
# print(df.loc[1, 'zero_zero_next'], type(df.loc[1, 'zero_zero_next']))
#
# for x in df.loc[1, 'one_one_back']:
#     print(x, type(x))

# print(df.to_string(max_rows=20))

# df['one_one'] = df.apply(lambda row: set(row['one_one']), axis=1)
# df['one_one_back'] = df.apply(lambda row: set(row['one_one_back']), axis=1)
# df['check'] = df.apply(lambda row: row['one_one'] in row['one_one_back'], axis=1)

# df['start_end'] = df.query("pos - pos_back > 10000")
# df['start_end'] = ['start' if x.pos - x.pos_back > 10000 else 0 for x in df]

# for x in df.itertuples():
#     # if x.POS > x.pos_back:
#     if x.one_one in x.one_one_back:
#         df['start_end'] = 'start'
#     else:
#         df['start_end'] = 0

df_part = df.drop(['ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'P1-25', 'P1-93', 'P1-88', 'P1-6',
                   'P1-89', 'P1-26', 'P1-12', 'P1-92', 'P1-22', 'P1-90', 'P1-28', 'P1-95'], axis=1)

print(df_part.to_string(max_rows=20))

# df_part['common_one_next'] = [set(a).intersection(b) for a, b in zip(df_part.one_one, df_part.one_one_back)]

# FOR FIND STARTS OF REGIONS OF HOMOZYGOSITY - METHOD 1
df_part.fillna(0)

notna_only = df_part[['one_one', 'zero_zero', 'pos_back', 'pos_next', 'one_one_back', 'zero_zero_back', 'one_one_next',
                      'zero_zero_next']].notna().all(axis=1)
df_part.loc[notna_only, 'one_one'] = [sorted(a) for a in df_part.loc[notna_only, 'one_one']]
df_part.loc[notna_only, 'zero_zero'] = [sorted(a) for a in df_part.loc[notna_only, 'zero_zero']]

df_part.loc[notna_only, 'inter1'] = [list(sorted(set(a).intersection(b))) for a, b in
                                     zip(df_part.loc[notna_only, 'one_one'], df_part.loc[notna_only, 'one_one_next'])]

df_part.loc[notna_only, 'inter2'] = [list(sorted(set(a).intersection(b))) for a, b in
                            zip(df_part.loc[notna_only, 'one_one'], df_part.loc[notna_only, 'zero_zero_next'])]

df_part.loc[notna_only, 'inter3'] = [list(sorted(set(a).intersection(b))) for a, b in
                            zip(df_part.loc[notna_only, 'zero_zero'], df_part.loc[notna_only, 'one_one_next'])]

df_part.loc[notna_only, 'inter4'] = [list(sorted(set(a).intersection(b))) for a, b in
                            zip(df_part.loc[notna_only, 'zero_zero'], df_part.loc[notna_only, 'zero_zero_next'])]

df_part.loc[notna_only, 'inter5'] = [list(sorted(set(a).intersection(b))) for a, b in
                                     zip(df_part.loc[notna_only, 'one_one'], df_part.loc[notna_only, 'one_one_back'])]

df_part.loc[notna_only, 'inter6'] = [list(sorted(set(a).intersection(b))) for a, b in
                            zip(df_part.loc[notna_only, 'one_one'], df_part.loc[notna_only, 'zero_zero_back'])]

df_part.loc[notna_only, 'inter7'] = [list(sorted(set(a).intersection(b))) for a, b in
                            zip(df_part.loc[notna_only, 'zero_zero'], df_part.loc[notna_only, 'one_one_back'])]

df_part.loc[notna_only, 'inter8'] = [list(sorted(set(a).intersection(b))) for a, b in
                            zip(df_part.loc[notna_only, 'zero_zero'], df_part.loc[notna_only, 'zero_zero_back'])]

# df_part.loc[notna_only, 'start'] = ['start' if a - b > 10000 and (len(c) >= 3 or len(d) >= 3) and
#                                                  (len(e) >= 3 or len(f) >= 3) else 0 for a, b, c, d, e, f
#                                         in zip(df_part.loc[notna_only, 'POS'], df_part.loc[notna_only, 'pos_back'],
#                                                df_part.loc[notna_only, 'inter1'], df_part.loc[notna_only, 'inter2'],
#                                                df_part.loc[notna_only, 'inter3'],
#                                                df_part.loc[notna_only, 'inter4'])]

# df_part.loc[notna_only, 'start'] =  ['start' if ((a-b) > 10000 and (c-a) < 10000) and (len(d) >= 3 or len(e) >= 3)
#                                     and (len(f) >= 3 or len(g) >= 3) and ((len(h) < 3 or len(i) < 3) or
#                                     (len(j) < 3 or len(k) < 3)) else 0 for a, b, c, d, e, f, g, h, i, j, k
#                                     in zip(df_part.loc[notna_only, 'POS'], df_part.loc[notna_only, 'pos_back'],
#                                     df_part.loc[notna_only, 'pos_next'], df_part.loc[notna_only, 'inter1'],
#                                     df_part.loc[notna_only, 'inter2'], df_part.loc[notna_only, 'inter3'],
#                                     df_part.loc[notna_only, 'inter4'], df_part.loc[notna_only, 'inter5'],
#                                     df_part.loc[notna_only, 'inter6'], df_part.loc[notna_only, 'inter7'],
#                                     df_part.loc[notna_only, 'inter8'])]

df_part.loc[notna_only, 'start'] = ['start' if (c-a) < 10000 and (len(d) >= 3 or len(e) >= 3)
                                    and (len(f) >= 3 or len(g) >= 3) and ((len(h) < 3 or len(i) < 3) or
                                    (len(j) < 3 or len(k) < 3)) else 0 for a, b, c, d, e, f, g, h, i, j, k
                                    in zip(df_part.loc[notna_only, 'POS'], df_part.loc[notna_only, 'pos_back'],
                                    df_part.loc[notna_only, 'pos_next'], df_part.loc[notna_only, 'inter1'],
                                    df_part.loc[notna_only, 'inter2'], df_part.loc[notna_only, 'inter3'],
                                    df_part.loc[notna_only, 'inter4'], df_part.loc[notna_only, 'inter5'],
                                    df_part.loc[notna_only, 'inter6'], df_part.loc[notna_only, 'inter7'],
                                    df_part.loc[notna_only, 'inter8'])]

# df_part.loc[notna_only, 'check'] = [set(a) <= set(b) for a, b in
#                             zip(df_part.loc[notna_only, 'one_one'], df_part.loc[notna_only, 'one_one_back'])]

# FOR FIND ENDS OF REGIONS OF HOMOZYGOSITY - METHOD 1
# df_part.loc[notna_only, 'end'] = ['end' if b - a > 10000 and (len(c) >= 3 or len(d) >= 3) and
#                                                  (len(e) >= 3 or len(f) >= 3) else 0 for a, b, c, d, e, f
#                                         in zip(df_part.loc[notna_only, 'POS'], df_part.loc[notna_only, 'pos_next'],
#                                                df_part.loc[notna_only, 'inter5'], df_part.loc[notna_only, 'inter6'],
#                                                df_part.loc[notna_only, 'inter7'],
#                                                df_part.loc[notna_only, 'inter8'])]

# df_part.loc[notna_only, 'end'] = ['end' if ((c-a) > 10000 and (a-b) < 10000) and (len(h) >= 3 or len(i) >= 3)
#                             and (len(j) >= 3 or len(k) >= 3) and ((len(d) < 3 or len(e) < 3) or len(f) < 3 or
#                                     len(g) < 3) else 0 for a, b, c, d, e, f, g, h, i, j, k
#                                     in zip(df_part.loc[notna_only, 'POS'], df_part.loc[notna_only, 'pos_back'],
#                             df_part.loc[notna_only, 'pos_next'], df_part.loc[notna_only, 'inter1'],
#                                     df_part.loc[notna_only, 'inter2'], df_part.loc[notna_only, 'inter3'],
#                                     df_part.loc[notna_only, 'inter4'], df_part.loc[notna_only, 'inter5'],
#                                     df_part.loc[notna_only, 'inter6'], df_part.loc[notna_only, 'inter7'],
#                                     df_part.loc[notna_only, 'inter8'])]

df_part.loc[notna_only, 'end'] = ['end' if (a-b) < 10000 and (len(h) >= 3 or len(i) >= 3) and (len(j) >= 3 or len(k) >= 3) and
                                           ((len(d) < 3 or len(e) < 3) or len(f) < 3 or len(g) < 3) else 0
                                  for a, b, c, d, e, f, g, h, i, j, k
                                  in zip(df_part.loc[notna_only, 'POS'], df_part.loc[notna_only, 'pos_back'],
                                  df_part.loc[notna_only, 'pos_next'], df_part.loc[notna_only, 'inter1'],
                                  df_part.loc[notna_only, 'inter2'], df_part.loc[notna_only, 'inter3'],
                                  df_part.loc[notna_only, 'inter4'], df_part.loc[notna_only, 'inter5'],
                                  df_part.loc[notna_only, 'inter6'], df_part.loc[notna_only, 'inter7'],
                                  df_part.loc[notna_only, 'inter8'])]

print(df_part.to_string(max_rows=320))

# df_part['start'] = df_part[df_part['start'] == 'start']
#
# df_start_end = df_part['chr', 'pos', 'start', 'end']

df_start_end = (df_part.loc[(df_part['start'] == 'start') | (df_part['end'] == 'end')]).\
    drop(['pos_back', 'pos_next', 'one_one_back', 'one_one_next', 'zero_zero_back', 'zero_zero_next', 'inter1',
          'inter2', 'inter3', 'inter4', 'inter5', 'inter6', 'inter7', 'inter8'], axis=1)

print(df_start_end.to_string(max_rows=150))

a = ['start', 'end']
df_bins = df_start_end[~(df_start_end['start'].isin(a) & (df_start_end['end'].isin(a)))]

print(df_bins.to_string(max_rows=150))
#
# rslt_df = dataframe.loc[(dataframe['Age'] == 21) &
#               dataframe['Stream'].isin(options)]

# df_start_end.to_csv('P1_bins_test_chr1_3.csv', sep='\t', index=False)

print("--- %s seconds ---" % (time.time() - start_time))
