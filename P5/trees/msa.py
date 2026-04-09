from tree_nodes import Node, newick2nodes
from Bio import SeqIO
from Bio.Align import substitution_matrices


## ALIGN PROFILES SCRIPT

subst_mat = substitution_matrices.load("BLOSUM62")

def align_profiles(profile_i, profile_j, gap):
    
    n_i = len(next(iter(profile_i.values())))
    n_j = len(next(iter(profile_j.values())))

    # Initialize scores and traceback matrices
    scores = [[0] * (n_j + 1) for _ in range(n_i + 1)]
    traceback = [[0] * (n_j + 1) for _ in range(n_i + 1)]

    # Initialize edges with gaps
    for i in range(1, n_i + 1):
        scores[i][0] = i * gap
        traceback[i][0] = 1  # up (gap in profile_j)

    for j in range(1, n_j + 1):
        scores[0][j] = j * gap
        traceback[0][j] = -1  # left (gap in profile_i)

    # Fill the matrices
    for i in range(1, n_i + 1):
        for j in range(1, n_j + 1):
            # Average substitution score across all sequence pairs, including gap penalties
            score = 0
            n_pairs = 0  # count of valid pairs (non-gap)
            for seq_i in profile_i.values():
                for seq_j in profile_j.values():
                    aa_i = seq_i[i - 1]
                    aa_j = seq_j[j - 1]
                    # score is an accumulator, add here the scores
                    
                    # FILL IN THE SCORE CALCULATION
                    if aa_i=="-" and aa_j=="-":
                        continue
                    elif (aa_i=="-" and aa_j!="-") or (aa_j=="-" and aa_i!="-"):
                        score+=gap
                        n_pairs+=1
                    else: 
                        n_pairs+=1
                        score+=subst_mat[aa_i.upper()][aa_j.upper()]

            score = score / n_pairs if n_pairs > 0 else 0  # avoid division by zero

            diag = scores[i - 1][j - 1] + score  # match/mismatch
            left = scores[i][j - 1] + gap  # insertion (gap in profile_i)
            up = scores[i - 1][j] + gap  # deletion (gap in profile_j)

            scores[i][j] = max(diag, left, up)

            if diag >= left and diag >= up:
                traceback[i][j] = 0  # diagonal
            elif left >= up:
                traceback[i][j] = -1  # left
            else:
                traceback[i][j] = 1  # up

    # Traceback
    i, j = n_i, n_j
    aln = {name: [] for name in {**profile_i, **profile_j}}

    while i > 0 or j > 0:
        if traceback[i][j] == 0:  # diagonal
            i -= 1
            j -= 1
            for name, seq in profile_i.items():
                aln[name].append(seq[i])
            for name, seq in profile_j.items():
                aln[name].append(seq[j])
        elif traceback[i][j] == -1:  # left (gap in profile_i)
            j -= 1
            for name in profile_i:
                aln[name].append("-")
            for name, seq in profile_j.items():
                aln[name].append(seq[j])
        else:  # up (gap in profile_j)
            i -= 1
            for name, seq in profile_i.items():
                aln[name].append(seq[i])
            for name in profile_j:
                aln[name].append("-")

    return round(scores[-1][-1], 2), {
        name: "".join(reversed(chars)) for name, chars in aln.items()
    }

## FASTA2DICT SCRIPT

def fasta2dict(filename):
   
    records = SeqIO.parse(filename, "fasta")
    d={}
    for r in records:
        d.update({r.id:r.seq})
    return d

## NODES2SPLITS_ALIGN SCRIPT
def node2splits_align(i_node, nodes, id2seq, gap):
    
    i2s=fasta2dict(id2seq)
    
    n=newick2nodes(nodes)

    node = n[i_node]
    left=[]
    right=[]
    # base case
    if node.name != "":
        return [node.name]

    # left
    if node.left != 0:
        left += node2splits_align(node.left, nodes, id2seq, gap)

    # right
    if node.right != 0:
        right += node2splits_align(node.right, nodes, id2seq, gap)

    align_profiles(i2s, left, right, gap)
    return left+right 