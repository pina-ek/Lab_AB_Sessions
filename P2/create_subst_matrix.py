## Estibaliz Pina - BBI - AB

import numpy as np
from Bio import SeqIO
from collections import Counter
c=Counter
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

## expected freq (counted / total)
s=["A", "B", "C"]
frq={
    "A": A/total, 
    "B": B/total, 
    "C": C/total,
    }

exp={
    "AA": frq["A"]**2, 
    "AB": 2*(frq["A"]*frq["B"]), # AB + BA
    "AC": 2*(frq["A"]*frq["C"]), # AC + CA
    "BB": frq["B"]**2,
    "BC": 2*(frq["A"]*frq["C"]), # BC + CB
    "CC": frq["C"]**2, 
    }

### Count pairs (combination per column of AA, AB, etc)
## Comparar una seq con la siguiente: 2 loops (iteracion 1: coge todas las seqs de la lista; iteracion 2: compara esa seq con las demas - no tiene que repetir las comparaciones!!!!!)
pairs_obs=[]
for i in range(len(seqs)): 
    for j in range(i+1, len(seqs)):
        s1=seqs[i] # seq 1 para comparacion - loop 1
        s2=seqs[j] # seq 2 para comparacion - loop 2
        for x in range(len(s1)): # selecciono la columna que quiero comparar
            p1=s1[x]
            p2=s2[x]
            pair = "".join(sorted([p1, p2]))
            pairs_obs.append(pair)

# ya tengo los pares ""hechos"", ahora falta contarlos
count_pairs=c(sorted(pairs_obs))
sp=sorted(count_pairs)
subsmat=np.zeros((3,3))

# pares observados - formula
obs_f={k: v / len(pairs_obs) for k, v in count_pairs.items()}

# subsmat con la formula - log(obs/exp)
lp=[]
for x in sp:
    freqs=obs_f[x]/exp[x]
    log_f=np.log10(freqs)
    l=np.around(log_f*10)
    lp.append(l)
z=0  
      
# solo la mitad de la matriz para que no se repitan valores
for i in range(len(s)):
    for j in range(len(s)):
        if i >=1 and j==0:
            subsmat[i,j]=0
        elif i==2 and j==1:
            subsmat[i,j]=0
        else:
            subsmat[i,j]=lp[z]
            z+=1

print(subsmat)




