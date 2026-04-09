from tree_nodes import Node, newick2nodes
from align_profiles_names import align_profiles

from fasta2dic import fasta2dict

def node2splits_align(i_node, nodes, id2seq, gap):
    """
    >>> node2splits_align(0, newick2nodes("hmgb.dnd"), fasta2dict("hmgb.fasta"), -4)
    ['hmgb_chite', 'hmgl_wheat', 'hmgl_trybr', 'hmgt_mouse']
    >>> node2splits_align(1, newick2nodes("hmgb.dnd"), fasta2dict("hmgb.fasta"), -4)
    ['hmgb_chite', 'hmgl_wheat', 'hmgl_trybr']
    >>> node2splits_align(2, newick2nodes("hmgb.dnd"), fasta2dict("hmgb.fasta"), -4)
    ['hmgb_chite', 'hmgl_wheat']
    
    """
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

a="(((hmgb_chite:0.10,hmgl_wheat:0.25):0.20,hmgl_trybr:0.60):0.25,hmgt_mouse:0.35);"
f="hmgb.fasta"
print(node2splits_align(i_node=0, nodes=a, id2seq=f, gap=-4))
   
