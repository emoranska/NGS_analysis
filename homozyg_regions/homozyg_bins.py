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

# print(df_part.to_string(max_rows=20))

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

# print(df_part.to_string(max_rows=320))

# df_part['start'] = df_part[df_part['start'] == 'start']
#
# df_start_end = df_part['chr', 'pos', 'start', 'end']

df_start_end = (df_part.loc[(df_part['start'] == 'start') | (df_part['end'] == 'end')]).\
    drop(['pos_back', 'pos_next', 'one_one_back', 'one_one_next', 'zero_zero_back', 'zero_zero_next', 'inter1',
          'inter2', 'inter3', 'inter4', 'inter5', 'inter6', 'inter7', 'inter8'], axis=1)

# print(df_start_end.to_string(max_rows=150))

a = ['start', 'end']
df_bins = df_start_end[~(df_start_end['start'].isin(a) & (df_start_end['end'].isin(a)))]

# print(df_bins.to_string(max_rows=200))

df_bins['one_and_zero'] = (df_bins['one_one'] + df_bins['zero_zero']).apply(sorted)
df_bins['one_zero_next'] = df_bins['one_and_zero'].shift(-1)
df_bins['one_zero_back'] = df_bins['one_and_zero'].shift(1)

# remove rows with the same samples in 1/1 and 0/0 inside the region of homozygosity
df_filter = df_bins.loc[~((df_bins['one_and_zero'] == df_bins['one_zero_next']) &
                          (df_bins['one_and_zero'] == df_bins['one_zero_back']))]

# print(df_filter.to_string(max_rows=50))

df_filter['one_one_next'] = df_filter['one_one'].shift(-1)
df_filter['zero_zero_next'] = df_filter['zero_zero'].shift(-1)
df_filter['one_one_back'] = df_filter['one_one'].shift(1)
df_filter['zero_zero_back'] = df_filter['zero_zero'].shift(1)
# df_filter['end_next'] = df_filter['end'].shift(-1)
# df_filter['pos_next'] = df_filter['POS'].shift(-1)
df_filter.fillna(0)

not_na = df_filter[['one_one_next', 'one_one_back', 'zero_zero_back', 'zero_zero_next']].notna().all(axis=1)

df_filter.loc[not_na, 'inter1'] = [list(sorted(set(a).intersection(b))) for a, b in
                                     zip(df_filter.loc[not_na, 'one_one'], df_filter.loc[not_na, 'one_one_next'])]

df_filter.loc[not_na, 'inter2'] = [list(sorted(set(a).intersection(b))) for a, b in
                                     zip(df_filter.loc[not_na, 'one_one'], df_filter.loc[not_na, 'zero_zero_next'])]

df_filter.loc[not_na, 'inter3'] = [list(sorted(set(a).intersection(b))) for a, b in
                                     zip(df_filter.loc[not_na, 'zero_zero'], df_filter.loc[not_na, 'one_one_next'])]

df_filter.loc[not_na, 'inter4'] = [list(sorted(set(a).intersection(b))) for a, b in
                                     zip(df_filter.loc[not_na, 'zero_zero'], df_filter.loc[not_na, 'zero_zero_next'])]

df_filter.loc[not_na, 'inter5'] = [list(sorted(set(a).intersection(b))) for a, b in
                                     zip(df_filter.loc[not_na, 'one_one'], df_filter.loc[not_na, 'one_one_back'])]

df_filter.loc[not_na, 'inter6'] = [list(sorted(set(a).intersection(b))) for a, b in
                                     zip(df_filter.loc[not_na, 'one_one'], df_filter.loc[not_na, 'zero_zero_back'])]

df_filter.loc[not_na, 'inter7'] = [list(sorted(set(a).intersection(b))) for a, b in
                                     zip(df_filter.loc[not_na, 'zero_zero'], df_filter.loc[not_na, 'one_one_back'])]

df_filter.loc[not_na, 'inter8'] = [list(sorted(set(a).intersection(b))) for a, b in
                                     zip(df_filter.loc[not_na, 'zero_zero'], df_filter.loc[not_na, 'zero_zero_back'])]


df_filter.loc[not_na, 'start2'] = ['start' if ((len(a) >= 3 or len(b) >= 3) and (len(c) >= 3 or len(d) >= 3)) and
                                              ((0 < len(e) < 3 or 0 < len(f) < 3) or (0 < len(g) < 3 or 0 < len(h) < 3)
                                               or ((len(e) >= 3 or len(f) >= 3) and (len(g) == 0 and len(h) == 0))) else
                                   '0' for a, b, c, d, e, f, g, h in zip(df_filter.loc[not_na, 'inter1'],
                                                                       df_filter.loc[not_na, 'inter2'],
                                                                       df_filter.loc[not_na, 'inter3'],
                                                                       df_filter.loc[not_na, 'inter4'],
                                                                       df_filter.loc[not_na, 'inter5'],
                                                                       df_filter.loc[not_na, 'inter6'],
                                                                       df_filter.loc[not_na, 'inter7'],
                                                                       df_filter.loc[not_na, 'inter8'])]

df_filter.loc[not_na, 'end2'] = ['end' if ((len(e) >= 3 or len(f) >= 3) and (len(g) >= 3 or len(h) >= 3)) and
                                          ((0 < len(a) < 3 or 0 < len(b) < 3) or (0 < len(c) < 3 or 0 < len(d) < 3) or
                                           ((len(a) >= 3 or len(b) >= 3) and (len(c) == 0 and len(d) == 0)) or
                                           ((len(a) >= 3 and len(b) == 0) and (len(c) >= 3 and len(d) == 0))) else '0' for
                                 a, b, c, d, e, f, g, h in zip(df_filter.loc[not_na, 'inter1'],
                                                               df_filter.loc[not_na, 'inter2'],
                                                               df_filter.loc[not_na, 'inter3'],
                                                               df_filter.loc[not_na, 'inter4'],
                                                               df_filter.loc[not_na, 'inter5'],
                                                               df_filter.loc[not_na, 'inter6'],
                                                               df_filter.loc[not_na, 'inter7'],
                                                               df_filter.loc[not_na, 'inter8'])]

df_filter_drop = df_filter.drop(['one_and_zero', 'one_zero_next', 'one_zero_back', 'one_one_back', 'one_one_next',
                                 'zero_zero_back', 'zero_zero_next', 'inter1', 'inter2', 'inter3', 'inter4', 'inter5',
                                 'inter6', 'inter7', 'inter8'], axis=1)

# print(df_filter_drop.to_string(max_rows=200))

b = ['0']
df_2_filter = (df_filter_drop.loc[~(df_filter_drop['start2'].isin(b) & df_filter_drop['end2'].isin(b))]).reset_index(drop=True)
df_2_filter.loc[0, 'start2'] = 'start'
df_2_filter.iloc[-1, 7] = 'end'

# df_2_filter = df_filter.loc[~((df_filter['inter1'] == df_filter['inter5']) & (df_filter['inter4'] == df_filter['inter8']))]

print(df_2_filter.to_string(max_rows=300))

df_2_filter['end2_ok'] = df_2_filter['end2'].shift(-1)
df_2_filter['pos_end2ok'] = df_2_filter['POS'].shift(-1)
df_2_filter['one_one_end2ok'] = df_2_filter['one_one'].shift(-1)
df_2_filter['zero_zero_end2ok'] = df_2_filter['zero_zero'].shift(-1)

df_2_filter_drop = df_2_filter.drop(['start', 'end'], axis=1)

print(df_2_filter_drop.to_string(max_rows=200))

df_final = (df_2_filter_drop.loc[(df_2_filter_drop['start2'].isin(a) & df_2_filter_drop['end2_ok'].isin(a))]).\
    reset_index(drop=True).drop(['end2'], axis=1)
print(df_final.to_string(max_rows=100))

df_final_2 = df_final.drop(['start2', 'end2_ok'], axis=1).rename(columns={'POS':'start', 'pos_end2ok':'end',
                                                                  'one_one': '1/1_start', 'zero_zero': '0/0_start',
                                                                  'one_one_end2ok': '1/1_end', 'zero_zero_end2ok': '0/0_end'}).\
    reindex(columns=['CHROM', 'start', 'end', '1/1_start', '0/0_start', '1/1_end', '0/0_end'])

print(df_final_2.to_string(max_rows=100))

# do dokończenia
# df_3_filter = df_filter.loc[(len(df_bins['inter1']) >= 3 or len(df_bins['inter2']) >= 3) and (len(df_bins['inter3']) >= 3 or len(df_bins['inter4']) >= 3)]

# df_filter = df_bins.loc[~(((df_bins['one_one'] == df_bins['one_one_next']) & (df_bins['zero_zero'] == df_bins['zero_zero_next'])
#                         & (df_bins['one_one'] == df_bins['one_one_back']) & (df_bins['zero_zero'] == df_bins['zero_zero_back'])) |
#                           ((df_bins['one_one'] == df_bins['zero_zero_next']) & (
#                                       df_bins['zero_zero'] == df_bins['one_one_next'])
#                            & (df_bins['one_one'] == df_bins['one_one_back']) & (
#                                        df_bins['zero_zero'] == df_bins['zero_zero_back']))
#                           )]


#
# rslt_df = dataframe.loc[(dataframe['Age'] == 21) &
#               dataframe['Stream'].isin(options)]

# df_2_filter.to_csv('P1_start_end_test_chr1_3.csv', sep='\t', index=False)

print("--- %s seconds ---" % (time.time() - start_time))
