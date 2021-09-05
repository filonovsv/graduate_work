data_file = open("all_ProtExpENSMUSG.txt")
start_line = data_file.readline()

set_genes = set()
dict_genes = dict()

for line in data_file:
    line = line.strip()
    words = line.split("\t")
    ensembl = words[1]
    set_genes.add(ensembl)
    zt = []
    for i in range(12):
        zt.append([float(i) for i in words[i+2::12]])
    dict_genes[ensembl] = zt

ppis_file = open("task_c.txt")
ppis = []
for line in ppis_file:
    line = line.strip()
    words = line.split("\t")
    ppis.append([words[0], words[1]])

files = []
for i in range(12):
    files.append(open("task_d_"+str(i*2)+".txt", "w"))

for ppi in ppis:
    sum = 0.0
    for i, j in zip(dict_genes[ppi[0]], dict_genes[ppi[1]]):
        sum += i[0]*j[0] + i[1]*j[1]
    average = sum/24.0
    for index, item in enumerate(zip(dict_genes[ppi[0]], dict_genes[ppi[1]])):
        i = item[0]
        j = item[1]
        if (i[0]*j[0] + i[1]*j[1])/2 > average:
            files[index].write(ppi[0] + "\t" + ppi[1] + "\n")
