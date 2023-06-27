from fuc import pyvcf
import pandas as pd
import time

start_time = time.time()

# open vcf file as pandas df
vf = pyvcf.VcfFrame.from_file('P1_chr1_3_homozyg_reg.vcf')
df = vf.df

# create columns with POS back and next
df['pos_back'] = df['POS'].shift(1)
df['pos_next'] = df['POS'].shift(-1)

# create columns with strings of samples with 1/1 and 0/0 respectively
cols = df.filter(like='P1-').columns
sep = ','
df['one_one'] = pd.Series(df[cols].eq('1/1').dot((cols + sep))).str.rstrip(sep)
df['zero_zero'] = pd.Series(df[cols].eq('0/0').dot((cols + sep))).str.rstrip(sep)
# print(df.to_string())

# create columns with lists of samples for 1/1, 0/0, next and back respectively
df['one_one_back'] = df['one_one'].shift(1).str.split(',')
df['one_one_next'] = df['one_one'].shift(-1).str.split(',')
df['zero_zero_back'] = df['zero_zero'].shift(1).str.split(',')
df['zero_zero_next'] = df['zero_zero'].shift(-1).str.split(',')
df['one_one'] = df['one_one'].str.split(',')
df['zero_zero'] = df['zero_zero'].str.split(',')

print(df.loc[1, 'one_one_back'], type(df.loc[0, 'one_one_back']))
# print(df.to_string(max_rows=20))

# remove useless columns
df_part = df.drop(['ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'P1-25', 'P1-93', 'P1-88', 'P1-6',
                   'P1-89', 'P1-26', 'P1-12', 'P1-92', 'P1-22', 'P1-90', 'P1-28', 'P1-95'], axis=1)
# print(df_part.to_string(max_rows=20))

df_part.fillna(0)
notna_only = df_part[['one_one', 'zero_zero', 'pos_back', 'pos_next', 'one_one_back', 'zero_zero_back', 'one_one_next',
                      'zero_zero_next']].notna().all(axis=1)
df_part.loc[notna_only, 'one_one'] = [sorted(a) for a in df_part.loc[notna_only, 'one_one']]
df_part.loc[notna_only, 'zero_zero'] = [sorted(a) for a in df_part.loc[notna_only, 'zero_zero']]

# create columns with intersections of current, back and next samples for 1/1 and 0/0 to be able to filter df
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

# find start points based of distance between next and current POS and length of lists with intersections of samples
# for 1/1 and 0/0
df_part.loc[notna_only, 'start'] = ['start' if (c-a) < 10000 and (len(d) >= 3 or len(e) >= 3)
                                    and (len(f) >= 3 or len(g) >= 3) and ((len(h) < 3 or len(i) < 3) or
                                    (len(j) < 3 or len(k) < 3)) else 0 for a, b, c, d, e, f, g, h, i, j, k
                                    in zip(df_part.loc[notna_only, 'POS'], df_part.loc[notna_only, 'pos_back'],
                                    df_part.loc[notna_only, 'pos_next'], df_part.loc[notna_only, 'inter1'],
                                    df_part.loc[notna_only, 'inter2'], df_part.loc[notna_only, 'inter3'],
                                    df_part.loc[notna_only, 'inter4'], df_part.loc[notna_only, 'inter5'],
                                    df_part.loc[notna_only, 'inter6'], df_part.loc[notna_only, 'inter7'],
                                    df_part.loc[notna_only, 'inter8'])]

# find ends based of distance between current and back POS and length of lists with intersections of samples for 1/1
# and 0/0
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

# filter records containting 'start' OR 'end' and drop useless columns
df_start_end = (df_part.loc[(df_part['start'] == 'start') | (df_part['end'] == 'end')]).\
    drop(['pos_back', 'pos_next', 'one_one_back', 'one_one_next', 'zero_zero_back', 'zero_zero_next', 'inter1',
          'inter2', 'inter3', 'inter4', 'inter5', 'inter6', 'inter7', 'inter8'], axis=1)
# print(df_start_end.to_string(max_rows=150))

# remove records with 'start' and 'end' simultaneously
a = ['start', 'end']
df_bins = df_start_end[~(df_start_end['start'].isin(a) & (df_start_end['end'].isin(a)))]

# create column with one common list of samples for 1/1 and 0/0 and then, with next and back list of samples
df_bins['one_and_zero'] = (df_bins['one_one'] + df_bins['zero_zero']).apply(sorted)
df_bins['one_zero_next'] = df_bins['one_and_zero'].shift(-1)
df_bins['one_zero_back'] = df_bins['one_and_zero'].shift(1)

# remove rows with the same samples for 1/1 and 0/0 inside the region of homozygosity (between start and end points)
df_filter = df_bins.loc[~((df_bins['one_and_zero'] == df_bins['one_zero_next']) &
                          (df_bins['one_and_zero'] == df_bins['one_zero_back']))]

# create columns with next or back samples for 1/1 and 0/0
df_filter['one_one_next'] = df_filter['one_one'].shift(-1)
df_filter['zero_zero_next'] = df_filter['zero_zero'].shift(-1)
df_filter['one_one_back'] = df_filter['one_one'].shift(1)
df_filter['zero_zero_back'] = df_filter['zero_zero'].shift(1)
df_filter.fillna(0)

# create columns with intersections of current, back and next samples for 1/1 and 0/0 to be able to filter df
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


# create new 'start' and 'end' columns for filtered df
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

# remove records containing 0 in 'start2' and 'end2' columns (NOT start and end of homozygosity regions)
b = ['0']
df_2_filter = (df_filter_drop.loc[~(df_filter_drop['start2'].isin(b) & df_filter_drop['end2'].isin(b))]).reset_index(drop=True)

# insert 'start' into the first record and 'end' for the last one
df_2_filter.loc[0, 'start2'] = 'start'
df_2_filter.iloc[-1, 7] = 'end'
# print(df_2_filter.to_string(max_rows=300))

# create columns with parameters of ends of hmozygosity regions (from next rows)
df_2_filter['end2_ok'] = df_2_filter['end2'].shift(-1)
df_2_filter['pos_end2ok'] = df_2_filter['POS'].shift(-1)
df_2_filter['one_one_end2ok'] = df_2_filter['one_one'].shift(-1)
df_2_filter['zero_zero_end2ok'] = df_2_filter['zero_zero'].shift(-1)

df_2_filter_drop = df_2_filter.drop(['start', 'end'], axis=1)
# print(df_2_filter_drop.to_string(max_rows=200))

# filter rows containing start in 'start2' and end in 'end2_ok'
df_final = (df_2_filter_drop.loc[(df_2_filter_drop['start2'].isin(a) & df_2_filter_drop['end2_ok'].isin(a))]).\
    reset_index(drop=True).drop(['end2'], axis=1)
# print(df_final.to_string(max_rows=100))

# create df with records as regions of homozygosity (chrom, start, end and samples for 1/1 and 0/0 start and end
# positions
df_final_2 = df_final.drop(['start2', 'end2_ok'], axis=1).rename(columns={'POS':'start', 'pos_end2ok':'end',
                                                                  'one_one': '1/1_start', 'zero_zero': '0/0_start',
                                                                  'one_one_end2ok': '1/1_end', 'zero_zero_end2ok': '0/0_end'}).\
    reindex(columns=['CHROM', 'start', 'end', '1/1_start', '0/0_start', '1/1_end', '0/0_end'])

# print(df_final_2.to_string(max_rows=100))

# create columns with intersection of samples for 1/1 (start and end) and 0/0 (start and end)
all_bins = df_final_2.all(axis=1)

df_final_2.loc[all_bins, 'inter_1/1'] = [list(sorted(set(a).intersection(b))) for a, b in
                                     zip(df_final_2.loc[all_bins, '1/1_start'], df_final_2.loc[all_bins, '1/1_end'])]

df_final_2.loc[all_bins, 'inter_0/0'] = [list(sorted(set(a).intersection(b))) for a, b in
                                     zip(df_final_2.loc[all_bins, '0/0_start'], df_final_2.loc[all_bins, '0/0_end'])]
# print(df_final_2.to_string(max_rows=100))

# create final df with records containing at least 3 samples for 1/1 and 0/0, the same for start and end
df_inter_11_00 = df_final_2[(df_final_2['inter_1/1'].map(len) >= 3) & (df_final_2['inter_0/0'].map(len) >= 3)].\
    reset_index(drop=True).drop(['1/1_start', '0/0_start', '1/1_end', '0/0_end'], axis=1).\
    rename(columns={'inter_1/1':'1/1','inter_0/0':'0/0'})

print(df_inter_11_00.to_string(max_rows=50))

# filter records NOT containing at least 3 samples for 1/1 and 0/0, the same for start and end
df_rest_inter_11_00 = df_final_2[~((df_final_2['inter_1/1'].map(len) >= 3) & (df_final_2['inter_0/0'].map(len) >= 3))]

# create column with intersection of samples for 1/1_start and 0/0_end (crosswise)
df_rest_inter_11_00.loc[all_bins, 'inter_1/1_0/0'] = [list(sorted(set(a).intersection(b))) for a, b in
                                     zip(df_rest_inter_11_00.loc[all_bins, '1/1_start'],
                                         df_rest_inter_11_00.loc[all_bins, '0/0_end'])]

# create columns with intersection of samples for 0/0_start and 1/1_end (crosswise)
df_rest_inter_11_00.loc[all_bins, 'inter_0/0_1/1'] = [list(sorted(set(a).intersection(b))) for a, b in
                                     zip(df_rest_inter_11_00.loc[all_bins, '0/0_start'],
                                         df_rest_inter_11_00.loc[all_bins, '1/1_end'])]

# print(df_rest_inter_11_00.to_string(max_rows=50))

# create final df containing records with at least 3 samples for 1/1 and 0/0 crosswise for start and end
df_inter_cross_11_00 = df_rest_inter_11_00[(df_rest_inter_11_00['inter_1/1_0/0'].map(len) >= 3) &
                                           (df_rest_inter_11_00['inter_0/0_1/1'].map(len) >= 3)].\
    reset_index(drop=True).drop(['1/1_start', '0/0_start', '1/1_end', '0/0_end', 'inter_1/1', 'inter_0/0'], axis=1).\
    rename(columns={'inter_1/1_0/0':'1/1','inter_0/0_1/1':'0/0'})

print(df_inter_cross_11_00.to_string(max_rows=50))

# filter records with 1/1 or 0/0 containing only 2 samples
df_other_11_00 = df_rest_inter_11_00[~((df_rest_inter_11_00['inter_1/1_0/0'].map(len) >= 3) &
                                       (df_rest_inter_11_00['inter_0/0_1/1'].map(len) >= 3))].reset_index(drop=True)
# print(df_other_11_00.to_string())

df_other_11_00_filter = df_other_11_00[(((df_other_11_00['inter_1/1'].map(len) >= 2) &
                                       (df_other_11_00['inter_0/0'].map(len) >= 2)) |
                                       ((df_other_11_00['inter_1/1_0/0'].map(len) >= 2) &
                                       (df_other_11_00['inter_0/0_1/1'].map(len) >= 2)))]
# print(df_other_11_00_filter.to_string())

# create columns with the rest of samples for records where 1/1 or 0/0 contains only 2 samples
df_other_11_00_filter['rest_inter_1/1'] = list(sorted(set(a) ^ set(b)) if len(c) == 2 else [] for a, b, c in
                                     zip(df_other_11_00_filter['1/1_start'],
                                         df_other_11_00_filter['1/1_end'], df_other_11_00_filter['inter_1/1']))

df_other_11_00_filter['rest_inter_0/0'] = list(sorted(set(a) ^ set(b)) if len(c) == 2 else [] for a, b, c in
                                     zip(df_other_11_00_filter['0/0_start'],
                                         df_other_11_00_filter['0/0_end'], df_other_11_00_filter['inter_0/0']))

df_other_11_00_filter['rest_inter_1/1_0/0'] = list(sorted(set(a) ^ set(b)) if len(c) == 2 else [] for a, b, c in
                                     zip(df_other_11_00_filter['1/1_start'],
                                         df_other_11_00_filter['0/0_end'], df_other_11_00_filter['inter_1/1_0/0']))

df_other_11_00_filter['rest_inter_0/0_1/1'] = list(sorted(set(a) ^ set(b)) if len(c) == 2 else [] for a, b, c in
                                     zip(df_other_11_00_filter['0/0_start'],
                                         df_other_11_00_filter['1/1_end'], df_other_11_00_filter['inter_0/0_1/1']))

# merge columns with rest of samples into one
df_other_11_00_filter['rest'] = (df_other_11_00_filter['rest_inter_1/1'] + df_other_11_00_filter['rest_inter_0/0'] +
                                 df_other_11_00_filter['rest_inter_1/1_0/0'] +
                                 df_other_11_00_filter['rest_inter_0/0_1/1']).apply(set).apply(sorted)

# merge columns for records containing only 2 samples for 1/1 or 0/0 with correct samples with 1/1 and 0/0 into one
df_other_11_00_filter['1/1'] = (df_other_11_00_filter['inter_1/1'] + df_other_11_00_filter['inter_1/1_0/0']).apply(set).apply(sorted)
df_other_11_00_filter['0/0'] = (df_other_11_00_filter['inter_0/0'] + df_other_11_00_filter['inter_0/0_1/1']).apply(set).apply(sorted)

# print(df_other_11_00_filter.to_string())

df_other_11_00_filter_final = df_other_11_00_filter.drop(['1/1_start', '0/0_start', '1/1_end', '0/0_end', 'inter_1/1',
                                                          'inter_0/0', 'inter_1/1_0/0', 'inter_0/0_1/1',
                                                          'rest_inter_1/1', 'rest_inter_0/0', 'rest_inter_1/1_0/0',
                                                          'rest_inter_0/0_1/1'], axis=1).reset_index(drop=True).\
    reindex(columns=['CHROM', 'start', 'end', '1/1', '0/0', 'rest'])
print(df_other_11_00_filter_final.to_string())
# df_2_filter.to_csv('P1_start_end_test_chr1_3.csv', sep='\t', index=False)

print("--- %s seconds ---" % (time.time() - start_time))
