from tree_nodes import Node, newick2nodes


def node2splits(i_node, nodes):
    """This function should call itself in a recursive way. The idea is to start
    exploring left then rigth. Takes as input the node id to start
    the recursion and a dictionary of node objects.
    >>> node2splits(0, newick2nodes("tree_abcdef.dnd"))
    ['A'] ['B']
    ['D'] ['E']
    ['C'] ['D', 'E']
    ['C', 'D', 'E'] ['F']
    ['A', 'B'] ['C', 'D', 'E', 'F']
    ['A', 'B', 'C', 'D', 'E', 'F']
    >>> node2splits(0, newick2nodes("hmgb.dnd"))
    ['hmgb_chite'] ['hmgl_wheat']
    ['hmgb_chite', 'hmgl_wheat'] ['hmgl_trybr']
    ['hmgb_chite', 'hmgl_wheat', 'hmgl_trybr'] ['hmgt_mouse']
    ['hmgb_chite', 'hmgl_wheat', 'hmgl_trybr', 'hmgt_mouse']
    """
    n=newick2nodes(nodes)
    node = n[i_node]
    left=[]
    right=[]
    # base case
    if node.name != "":
        return [node.name]

    # left
    if node.left != 0:
        left += node2splits(node.left, nodes)

    # right
    if node.right != 0:
        right += node2splits(node.right, nodes)

    print(left, right)
    return left+right 

    
   
print(node2splits(0, "((A,B),((C,(D,E)),F));"))
print(node2splits(0, "(((hmgb_chite:0.10,hmgl_wheat:0.25):0.20,hmgl_trybr:0.60):0.25,hmgt_mouse:0.35);"))