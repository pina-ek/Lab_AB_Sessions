import numpy as np


seq1="FAST"
seq2="CAT"
mat=np.zeros((len(seq2)+1, len(seq1)+1))
arrw=np.zeros((len(seq2)+1, len(seq1)+1))
gap=-2
g=-2
for i in range(len(seq2)+1):
    for j in range(len(seq1)+1):
        if i==0 and j==0:
            continue
        elif i==0 and j>0:
            mat[i,j]=gap
            arrw[i,j]=-1
            gap-=2
        elif i>0 and j==0:
            mat[i,j]=g
            arrw[i,j]=1
            g-=2
        
print(mat)
print(arrw)
# hasta aquí está bien!!!!

import argparse as a
import numpy as np

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



----------------
arrw[-1, -1]=0
for i in range(len(s1)-1, 1,-1):
    for j in range(len(s2)-1, 1,-1):
        diag=arrw[i-1, j-1]
        up=arrw[i-1, j]
        left=arrw[i, j-1]
        point=max(diag, up, left)

        if point==diag:
            diag=0
            up=1
            left=1
        elif point==left:
            diag=1
            left=0
            up=1
        elif point==up:
            diag=1
            left=1
            up=0            


if j<len(s2)-1 and i<len(s1)-1:
            rd=mat[i+1, j+1]
            z=mat[i, j]
            point=max(z, rd, up)
            if s1[i]!=s2[j]:
                if s1[i+1]==s2[j+1] and s1[i-1]==s2[j-1]:
                    if point==z:
                        arrw[i,j]=0
                    elif point==rd:
                        arrw[i+1, j+1]=0
                        arrw[i,j]=1
                    elif point==up:
                        arrw[i-1, j]=1
                        arrw[i,j]=1




