from scipy.stats import f_oneway

data_file = open("all_ProtExpENSMUSG.txt")
start_line = data_file.readline()

proteins = []
for line in data_file:
    if line == "\n":
        continue
    line = line.strip()
    words = line.split("\t")
    number = words.pop(0)
    ensembl = words.pop(0)
    zt = []
    for i in range(12):
        zt.append([float(i) for i in words[i::12]])
    f, p = f_oneway(*zt)
    if p < 0.001:
        proteins.append((number, ensembl, words, p))

output_file = open("task_a.txt", "w")

for i in proteins:
    output_file.write(i[1] + "\n")