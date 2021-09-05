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
    ensemble = words.pop(0)
    zt = []
    for i in range(12):
        zt.append([float(i) for i in words[i::12]])
    f, p = f_oneway(*zt)
    if p < 0.001:
        proteins.append((number, ensemble, words, p))

times = []

for protein in proteins:
    time = ""
    sum = 0.0
    for i in range(24):
        sum += float(protein[2][i])
    avg = sum / 24
    for i in range(12):
        if (float(protein[2][i]) + float(protein[2][i+12]))/2 > avg:
            time += "\t" + str(i*2)
    times.append((protein[0], protein[1], time))

output_file = open("task_b.txt", "w")
for protein in times:
    output_file.write(protein[1] + protein[2] + "\n")

