import pandas as pd
import csv
from fuc import pyvcf
import re
import numpy as np
from less_memory_usage import remove_the_same_gt, to_csv
from vcf_open import remove_the_same_gt_vcf

# dostosować regex do outputu z .vcf, popatrzeć też na wynik GT strip i ewentualnie zmienić drop.nan


no_the_same_variant = remove_the_same_gt_vcf('P1_chr4_test.vcf')

no_the_same_variant.astype("string")


# replace missing values signed as '.' with 'NaN'
ntsv_nan = no_the_same_variant.replace('.', np.nan)


# remove records with >=6 samples with missing values ('.' replaced with NaN)
drop_dots = ntsv_nan.dropna(thresh=11, subset=['P1-25', 'P1-93', 'P1-88', 'P1-6', 'P1-89', 'P1-26', 'P1-12', 'P1-92',
                                                'P1-22', 'P1-90', 'P1-28', 'P1-95'])

drop_dots_str = drop_dots.to_string()
# print(type(drop_dots), type(drop_dots_str), drop_dots_str)
print(drop_dots_str)


# complete regex in methods.txt, here just part of it for testing
# regex = r"^\d\s{3}chr\d?\s{2}\d+\s+\D+\d+\.\d+\s+\w+"
regex = r"^\d\s{3}(chr\d)\s{2}\d+\s+\D+\d+\.\d+\s+\w+(?:(?:\s{2,}\d/\d)*?\s{2,}([10])/\2\b){3}.*"
# regex = r"(?=(^\d\s{3}(chr\d)\s{2}\d+\s+\D+\d+\.\d+\s+\w+(?=(?:(?:\s+[01]/[01])*?\s+1/1){3})(?:(?:(?:\s+[01]/[01])*?\s+0/0){3})(?:\s+[01]/[01])*\s(^\d\s{3}((chr\d)\s{2}\d+\s+\D+\d+\.\d+\s+\w+(?:\s+[01]/[01])*\s)){1,2}^\d\s{3}(chr\d)\s{2}\d+\s+\D+\d+\.\d+\s+\w+(?=(?:(?:\s+[01]/[01])*?\s+1/1){3})(?:(?:(?:\s+[01]/[01])*?\s+0/0){3})(?:\s+[01]/[01])*\s))|(^\d\s{3}(chr\d)\s{2}\d+\s+\D+\d+\.\d+\s+\w+(?=(?:(?:\s+[01]/[01])*?\s+1/1){3})(?:(?:(?:\s+[01]/[01])*?\s+0/0){3})(?:\s+[01]/[01])*$)"
matches = re.finditer(regex, drop_dots_str, re.MULTILINE)

regions_homozyg = []
for matchNum, match in enumerate(matches, start=1):
    regions_homozyg.append(match.group())
    # print(match.group())
    # print("Match {matchNum} was found at {start}-{end}: {match}".format(matchNum=matchNum, start=match.start(),
    #                                                                     end=match.end(), match=match.group()))

    # for groupNum in range(0, len(match.groups())):
    #     groupNum = groupNum + 1
    #
    #     print("Group {groupNum} found at {start}-{end}: {group}".format(groupNum=groupNum, start=match.start(groupNum),
    #                                                                     end=match.end(groupNum),
    #                                                                     group=match.group(groupNum)))
# print(type(matches), matches)
print(type(regions_homozyg))


# Open File
resultFile = open("P1_chr4_output.csv", 'w')

# Write data to file
for r in regions_homozyg:
    resultFile.write(r + "\n")
resultFile.close()
