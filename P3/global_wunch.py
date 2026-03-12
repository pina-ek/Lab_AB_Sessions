# E. Pina - BBI - AP - Practica 3


import numpy as np

# parametros para el programa
seq1="FAST"
seq2="CAT"
m=1 
mm=-1
g=ga=-2
# -------------------------------------
# Limites de la matriz (gap penalty para todos los limites)
s1=" "+seq1
s2=" "+seq2
mat=np.zeros((len(s1), len(s2)))
arrw=np.zeros((len(s1), len(s2)))


for i in range(len(s1)):
    for j in range(len(s2)):
        if i==0 and j==0:
            continue
        elif i==0 and j>0:
            mat[i,j]=ga
            arrw[i,j]=-1
            ga-=2
        elif i>0 and j==0:
            mat[i,j]=g
            arrw[i,j]=1
            g-=2

# ------------------------------------
# Wunch matrix
matx=m
mismatx=mm
gap=g

for i in range(1, len(s1)):
    for j in range(1, len(s2)):
        diag=mat[i-1, j-1]+matx
        up=mat[i-1, j]+gap
        left=mat[i, j-1]+gap

        if s1[i]==s2[j]:
            mat[i,j]=max(diag+2, up, left)
            arrw[i,j]=0
            
            if s1[i-1]!=s2[j-1]:
                
                point=max(diag, up, left)
                if point==diag:
                    arrw[i-1, j-1]=0
                elif point==left:
                    arrw[i, j-1]=0
                elif point==up:
                    arrw[i-1, j]=0

        elif s1[i]!=s2[j]:
            mat[i,j]= max(diag-1, up, left)
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
    fila_usada = sq1
elif sq2.shape[0] == mat.shape[1]:
    mat1 = np.vstack([sq2, mat])
    arrw1 = np.vstack([sq2, arrw])
    fila_usada = sq2
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


print(mat2)
print("-----")
print(arrw2)
