from tree_nodes import Node, newick2nodes


def sisters(nodes):
    
    n=newick2nodes(nodes)
       
    result=[]
    for i in range(len(n)):
        if n[i].name!="":
            if n[i].parent==i-1:
                left_child=n[i-1].left
                right_child=n[i-1].right
                if n[left_child].name!="" and n[right_child].name!="":
                    result.append((left_child, right_child, n[i].parent))
           
    return result

print(sisters("((A,B),((C,(D,E)),F));"))
print(sisters("(((hmgb_chite:0.10,hmgl_wheat:0.25):0.20,hmgl_trybr:0.60):0.25,hmgt_mouse:0.35);"))
print(sisters("((A,B),(((C,D),E),F));"))



# python -m doctest align_profiles.py -v