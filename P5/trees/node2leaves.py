from tree_nodes import Node, newick2nodes


def node2leaves(i_node: int, nodes) -> list[str]:
    """This function should call itself in a recursive way. The idea is to start
    exploring left then right. Takes as input the node id to start the recursion
    and a dictionary of nodes. Returns a list of leaf names.

    >>> node2leaves(0, newick2nodes("((A,B),((C,(D,E)),F));"))
    ['F', 'E', 'D', 'C', 'B', 'A']
    >>> node2leaves(1, newick2nodes("((A,B),((C,(D,E)),F));"))
    ['B', 'A']
    >>> node2leaves(2, newick2nodes("((A,B),((C,(D,E)),F));"))
    ['A']
    >>> node2leaves(4, newick2nodes("((A,B),((C,(D,E)),F));"))
    ['F', 'E', 'D', 'C']
    >>> node2leaves(7, newick2nodes("((A,B),((C,(D,E)),F));"))
    ['E', 'D']
    """
    node = nodes[i_node]
    left=[]
    right=[]
    # base case
    if node.name != "":
        return [node.name]

    # left
    if node.left != 0:
        left += node2leaves(node.left, nodes)

    # right
    if node.right != 0:
        right += node2leaves(node.right, nodes)

    return left+right   
    

        

    


print(node2leaves(0, newick2nodes("((A,B),((C,(D,E)),F));")))
print(node2leaves(1, newick2nodes("((A,B),((C,(D,E)),F));")))
print(node2leaves(2, newick2nodes("((A,B),((C,(D,E)),F));")))
print(node2leaves(4, newick2nodes("((A,B),((C,(D,E)),F));")))
print(node2leaves(7, newick2nodes("((A,B),((C,(D,E)),F));")))
