from Bio import Phylo
from io import StringIO

newick_str = "((A,B),((C,(D,E)),F));"
handle = StringIO(newick_str)
tree = Phylo.read(handle, "newick")
Phylo.draw_ascii(tree)

from tree_nodes import Node, newick2nodes

nodes = newick2nodes("((A,B),((C,(D,E)),F));")
nodes
