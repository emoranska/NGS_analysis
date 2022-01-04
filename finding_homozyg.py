import pandas as pd
import numpy as np
from less_memory_usage import remove_the_same_gt, to_csv


just_gt = remove_the_same_gt('data1.csv', 'data2.csv')
# data3 = to_csv(just_gt, 'data3_concat.csv')

just_gt.astype("string")


# for row in range(len(just_gt)):
#     for x in just_gt.iloc[row, 7:]:
#         print(x)
#     print(row)

# Pandas W3Schools Tutorial
#
# Loop through all values in the "Duration" column.
# If the value is higher than 120, set it to 120:

# for x in just_gt.index:
#     if just_gt.loc[x, 'P1-12'] == '.':
#         just_gt.drop(x, inplace=True)

just_gt_nan = just_gt.replace('.', np.nan)

# print(just_gt_nan.to_string())

# remove records with >=6 samples calls no data ('.' replaced with NaN)
drop_dots = just_gt_nan.dropna(thresh=11, subset=['P1-25', 'P1-93', 'P1-88', 'P1-6', 'P1-89', 'P1-26', 'P1-12', 'P1-92',
                                                'P1-22', 'P1-90', 'P1-28', 'P1-95'])

print(drop_dots.to_string())

# drop_dots.loc[row] --> Series dtype --> we can use regex!!!
for row in drop_dots.index:
    print(drop_dots.loc[row].str.extract(r"(0/1)"))
    # print(drop_dots.loc[row])



# print(drop_dots.str.extract("0/1", expand=True))
# check if in the particular position there are >=3 samples with 0/0 and >= samples with 1/1
# def find_homozyg():


# Delete rows where "Duration" is higher than 120:
# for x in df.index:
#   if df.loc[x, "Duration"] > 120:
#     df.drop(x, inplace = True)