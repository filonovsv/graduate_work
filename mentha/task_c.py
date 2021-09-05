from scipy.stats import f_oneway
import mygene

ppis_file = open("10090")
ppis_file.readline()

genes = []
for line in ppis_file:
    if line == "\n":
        continue
    line = line.strip()
    words = line.split(";")
    gene_A = words[1]
    gene_B = words[3]
    genes.append(gene_A)
    genes.append(gene_B)


mg = mygene.MyGeneInfo()
genes_mg = mg.querymany(genes, scopes='symbol', fields='ensembl.gene', species='mouse', size=1)

genes_ensembl = []

for i in genes_mg:
    try:
        ensembl = i.pop("ensembl", "null")
        if ensembl != "null":
            ensembl = ensembl.pop("gene", "null")
        genes_ensembl.append(ensembl)
    except:
        print(i)
        print("\n")
        genes_ensembl.append("null")

ppis = []
temp = ""
for index, gene in enumerate(genes_ensembl):
    if index % 2 == 0:
        temp = gene
    else:
        if gene != "null" and temp != "null":
            ppis.append((temp, gene))
        else:
            temp = ""

data_file = open("all_ProtExpENSMUSG.txt")
start_line = data_file.readline()

set_genes = set()
dict_genes = dict()

for line in data_file:
    if line == "\n":
        continue
    line = line.strip()
    words = line.split("\t")
    ensembl = words[1]
    set_genes.add(ensembl)
    zt = []
    for i in range(12):
        zt.append([float(i) for i in words[i+2::12]])
    dict_genes[ensembl] = zt

ppis_measured = []
for i in ppis:
    if i[0] in set_genes and i[1] in set_genes:
        ppis_measured.append(i)

output_file = open("task_c.txt", "w")
for ppi in ppis_measured:
    zt = []
    for i, j in zip(dict_genes[ppi[0]], dict_genes[ppi[1]]):
        zt.append([i[0]*j[0], i[1]*j[1]])
    f, p = f_oneway(*zt)
    if p < 0.001:
        output_file.write(ppi[0] + "\t" + ppi[1] + "\n")

