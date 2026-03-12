# E. Pina - BBI - AP - Practica 3

import argparse as a
import numpy as np
from Bio.Align import PairwiseAligner


# parametros para el programa
pars=a.ArgumentParser()
pars.add_argument("seq1", help="give 1st string to compare", type=str)
pars.add_argument("seq2", help="give 2nd string to compare", type=str)
pars.add_argument("-ma", "--match", help="score given to matching strings", type=int)
pars.add_argument("-mm", "--mismatch", help="score given to mismatching strings", type=int)
pars.add_argument("-g", "--gap", help="gap penalty score", type=int)
args=pars.parse_args()
# -------------------------------------
# Limites de la matriz (SIN penalty para todos los limites)
s1=" "+args.seq1
s2=" "+args.seq2
matx=args.match
mismatx=args.mismatch
gap=args.gap
c=g=ga=gap
mat=np.zeros((len(s1), len(s2)))

matching_chars_s1=[]
matching_chars_s2=[]

# Alignment matrix: REINICIALIZAR MATRIX POR CADA VALOR NEGATIVO!!!!!

for i in range(1, len(s1)):
    for j in range(1, len(s2)):
        
        up=mat[i-1, j]+gap
        left=mat[i, j-1]+gap
        if s1[i] == s2[j]:
            diag = mat[i-1, j-1] + matx
        else:
            diag = mat[i-1, j-1] + mismatx

        mat[i,j] = max(0, diag, up, left)
            


## A ponerlo bonito
list_seq1=[]
list_seq2=[]
for i in s1:
    list_seq1.append(i)
for j in s2:
    list_seq2.append(j)
sq1=np.array(list_seq1)
sq2=np.array(list_seq2)

corner = " " 

# row 0
if sq1.shape[0] == mat.shape[1]:
    mat1 = np.vstack([sq1, mat])
    
elif sq2.shape[0] == mat.shape[1]:
    mat1 = np.vstack([sq2, mat])
    
else:
    raise ValueError("NOT A MATCH FOR THIS MATRIX -- ROW")

# Column 0
if sq2.shape[0] == mat.shape[0]:
    col = sq2
elif sq1.shape[0] == mat.shape[0]:
    col = sq1
else:
    raise ValueError("NOT A MATCH FOR THIS MATRIX -- COLUMN")

col0 = np.insert(col, 0, corner)[:, None]
mat2  = np.hstack([col0, mat1])


align1 = []
align2 = []

i, j = np.unravel_index(np.argmax(mat), mat.shape)

while mat[i, j] != 0:

    if s1[i] == s2[j]:
        score_diag = mat[i-1, j-1]+ matx
    else:
        score_diag = mat[i-1, j-1]+mismatx

    if mat[i, j] == score_diag:
        align1.append(s1[i])
        align2.append(s2[j])
        i -= 1
        j -= 1

    elif mat[i, j] == mat[i-1, j]+gap:
        align1.append(s1[i])
        align2.append("-")
        i -= 1

    else:
        align1.append("-")
        align2.append(s2[j])
        j -= 1

# reverse the order 
align1 = "".join(reversed(align1))
align2 = "".join(reversed(align2))

aligner=PairwiseAligner()
aligner.mode="local"
aligner.match_score=args.match
aligner.mismatch_score=args.mismatch
aligner.open_gap_score=args.gap
aligner.extend_gap_score=args.gap 
opt_score=aligner.score(args.seq1, args.seq2)

print(" ")
print(f">> Global alignment matrix:"+"\n")
print(mat2)
print(" ")
print("------------")
print(" ")
print(f">> Sequence alignment:"+"\n")
print(f"Optimal score: {opt_score}")
print(" ")
print(align1)
print(align2)