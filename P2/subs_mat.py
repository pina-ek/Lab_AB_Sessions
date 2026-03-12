import numpy as np
from Bio import SeqIO
from collections import Counter
c=Counter()
A=0
B=0 
C=0

seqs=[]
total=0
## open file, read file, extract seq and convert it into str
## count how many As, Bs, Cs there are and the total count
with open("baba.fasta", "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        seq=str(record.seq)
        for i in seq:
            if i=="A":
                A+=1
            elif i=="B":
                B+=1
            elif i=="C":
                C+=1
            total+=1
        seqs.append(seq)
## expected freq individual
p = {
    "A": A / total,
    "B": B / total,
    "C": C / total
}

## Pairs observed
pair_total = 0

for i in range(len(seqs)):
    for j in range(i + 1, len(seqs)):
        s1 = seqs[i]
        s2 = seqs[j]
        for k in range(len(s1)):
            a, b = s1[k], s2[k]
            pair = "".join(sorted([a, b]))  # AB y BA -> AB
            c[pair] += 1
            pair_total += 1

# Frecuencias observadas
f_obs = {k: v / pair_total for k, v in c.items()}

# Frecuencias esperadas
f_exp = {
    "AA": p["A"] ** 2,
    "AB": 2 * p["A"] * p["B"],
    "AC": 2 * p["A"] * p["C"],
    "BB": p["B"] ** 2,
    "BC": 2 * p["B"] * p["C"],
    "CC": p["C"] ** 2,
}

ss = ["A", "B", "C"]
subsmat = np.zeros((3, 3))

for i, s1 in enumerate(ss):
    for j, s2 in enumerate(ss):
        pair = "".join(sorted([s1, s2]))
        if pair in f_obs and f_obs[pair] > 0:
            subsmat[i, j] = np.log10(f_obs[pair] / f_exp[pair])
        else:
            subsmat[i, j] = " "

su=np.around(subsmat*10)
print(su)
