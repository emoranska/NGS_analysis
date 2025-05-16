import csv

fasta = {}
with open("../files/RefBeetEL10_cleaned_db_27.11.21.fas", "r") as file:
    for line in file:
        line = line.strip()
        print(type(line), line)
        print(len(line), line[0])
        if line[0] == ">":
            header = line
            fasta[header] = ""
        else:
            data = line
            fasta[header] += data

mite_lengths = []
for k in fasta.keys():
    te_length = k, len(fasta[k])
    mite_lengths.append(te_length)

    # print(f"{k}, {len(fasta[k])}")
print(mite_lengths)

with open('../files/MITE_consensus_seq_lengths.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerows(mite_lengths)
