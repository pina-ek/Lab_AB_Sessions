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
# Limites de la matriz (gap penalty para todos los limites)
s1=" "+args.seq1
s2=" "+args.seq2
matx=args.match
mismatx=args.mismatch
gap=args.gap
c=g=ga=gap
mat=np.zeros((len(s1), len(s2)))
arrw=np.zeros((len(s1), len(s2)))
matching_chars_s1=[]
matching_chars_s2=[]

for i in range(len(s1)):
    for j in range(len(s2)):
        if i==0 and j==0:
            continue
        elif i==0 and j>0:
            mat[i,j]=ga
            arrw[i,j]=-1
            ga-=c
        elif i>0 and j==0:
            mat[i,j]=g
            arrw[i,j]=1
            g-=c

# Alignment matrix: BOTH ARROWS AND NUMERS

for i in range(1, len(s1)):
    for j in range(1, len(s2)):
        diag=mat[i-1, j-1]
        up=mat[i-1, j]+gap
        left=mat[i, j-1]+gap

        if s1[i]==s2[j]:
            mat[i,j]=max(diag+matx, up, left)
            arrw[i,j]=0
            
            
            if s1[i-1]!=s2[j-1]:
                
                point=max(diag+matx, up, left)
                if point==diag:
                    arrw[i-1, j-1]=0
                elif point==left:
                    arrw[i, j-1]=0
                elif point==up:
                    arrw[i-1, j]=0

        elif s1[i]!=s2[j]:
            mat[i,j]= max(diag+mismatx, up, left)
            arrw[i,j]=1
            

    
## Hasta aquí funciona pero solo con 1 gap de diferencia
## No funciona para secuencias como FAST - FART, o como FA--T - FAT

# arreglar la matriz de arrows a posteriori (??): ir directamente a los zeros y
# utilizar la matriz mat para ver para donde tirar

for i in range(len(s1)-1, 1,-1):
    for j in range(len(s2)-1, 1,-1):
        diag=mat[i-1, j-1]
        up=mat[i-1, j]
        left=mat[i, j-1]
        point=max(diag, left, up)
        d= arrw[i-1, j-1]
        l= arrw[i, j-1]
        u=arrw[i-1, j]

        if arrw[i,j]==0:
            
            if d==l==u==1:
                if point==diag and point!=left and point!=up:
                    arrw[i-1, j-1]=0
                elif point==left and point!=diag and point!=up:
                    arrw[i, j-1]=0
                    
                    
                elif point==up and point!=diag and point!=left:
                    arrw[i-1, j]=0
                    
                elif point==left==up:
                    if len(s1)>len(s2):
                        arrw[i-1, j]=0
                        
                    elif len(s2)>len(s1):
                        arrw[i, j-1]=0
                        

            else:
                continue
        else:
            continue

## ARREGLAO' !!!!!!!

# Alineamiento de secuencias con el traceback: HAY QUE REVERTIR EL ORDEN!!!


align1 = []
align2 = []

i = len(s1) - 1
j = len(s2) - 1

while i > 0 or j > 0:

    if i > 0 and j > 0:
        score_current = mat[i, j]
        score_diag = mat[i-1, j-1] + (matx if s1[i] == s2[j] else mismatx)
        score_up = mat[i-1, j] + gap
        score_left = mat[i, j-1] + gap

        if score_current == score_diag:
            align1.append(s1[i])
            align2.append(s2[j])
            i -= 1
            j -= 1
            continue

        if score_current == score_up:
            align1.append(s1[i])
            align2.append("-")
            i -= 1
            continue

        if score_current == score_left:
            align1.append("-")
            align2.append(s2[j])
            j -= 1
            continue

    elif i > 0:
        align1.append(s1[i])
        align2.append("-")
        i -= 1
    else:
        align1.append("-")
        align2.append(s2[j])
        j -= 1

align1 = "".join(reversed(align1))
align2 = "".join(reversed(align2))

## Ahora hay que ponerlo bonito
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
    arrw1 = np.vstack([sq1, arrw])
    
elif sq2.shape[0] == mat.shape[1]:
    mat1 = np.vstack([sq2, mat])
    arrw1 = np.vstack([sq2, arrw])
    
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
arrw2 = np.hstack([col0, arrw1])



aligner=PairwiseAligner()
aligner.mode="global"
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
print(f">> Traceback matrix:"+"\n")
print(arrw2)
print(" ")
print("------------")
print(" ")
print(f">> Sequence alignment:"+"\n")
print(f"Optimal score: {opt_score}")
print(" ")
print(align1)
print(align2)

