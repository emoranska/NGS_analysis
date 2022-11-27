import pandas as pd
import numpy as np

data = pd.read_csv('chr4_breaks_testdata.csv')

i = data.apply(lambda r: len(set(r[9:])) != 1, axis=1)
just_gt = data[i]

just_gt.astype("string")

# replace missing values signed as '.' with 'NaN'
just_gt_nan = just_gt.replace('.', np.nan)

# remove records with >=6 samples with missing values ('.' replaced with NaN)
drop_dots = just_gt_nan.dropna(thresh=7, subset=['P1-25', 'P1-93', 'P1-88', 'P1-6', 'P1-89', 'P1-26', 'P1-12', 'P1-92',
                                                 'P1-22', 'P1-90', 'P1-28', 'P1-95'])
drop_dots_str = drop_dots.to_string()

# checking if in row there are >= 3 1/1 AND >=3 0/0
lines = drop_dots_str.splitlines()
result = []
for line in lines:
    cols = line.split()
    samples = cols[8:]
    print(samples)
    homo_1 = 0
    homo_2 = 0
    for pos in samples:
        if pos in "1/1":
            homo_1 += 1
        elif pos in "0/0":
            homo_2 += 1
    if homo_1 >= 3 and homo_2 >= 3:
        result.append(line)
    print(homo_1, homo_2)
print(*result, sep="\n")

# open csv file
resultFile = open("chr_test_output.csv", 'w')

# write data to file
for row in result:
    resultFile.write(row + "\n")
resultFile.close()

# # solution from CAPSLOCK
# cond_11x3_true = (drop_dots.iloc[:, 8:] == "1/1").sum(axis=1) >= 3
# cond_00x3_true = (drop_dots.iloc[:, 8:] == "0/0").sum(axis=1) >= 3
# result = drop_dots[cond_11x3_true & cond_00x3_true] # pick
# print(result)
#
# extra_2_rows = drop_dots[(cond_11x3_true & ~cond_00x3_true) | (~cond_11x3_true & cond_00x3_true)].head(2)
# result.append(extra_2_rows)
# print(result)
