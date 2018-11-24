from ete3 import Tree, TreeStyle, TextFace, add_face_to_node, SeqMotifFace

tree = Tree("/Users/gabefoley/Dropbox/PhD/Projects/Phylo Island/Species trees/species_tree.nwk")

type1 = open("/Users/gabefoley/Dropbox/PhD/Projects/Phylo Island/Species trees/Type1.txt")
type2a = open("/Users/gabefoley/Dropbox/PhD/Projects/Phylo Island/Species trees/Type2A.txt")
type2b = open("/Users/gabefoley/Dropbox/PhD/Projects/Phylo Island/Species trees/Type2B.txt")
multiple = open("/Users/gabefoley/Dropbox/PhD/Projects/Phylo Island/Species trees/Multiple.txt")



type1_list = []
type2a_list = []
type2b_list = []
multiple_list = []

colours = {"type1": "green", "type2a": "blue", "type2b": "yellow", "multiple": "brown"}


for name in type1:
    type1_list.append(name)

for name in type2a:
    type2a_list.append(name)

for name in type2b:
    type2b_list.append(name)

for name in multiple:
    multiple_list.append(name)

for leaf in tree:
    print (leaf.name)

ts = TreeStyle()
ts.show_leaf_name = False
ts.show_leaf_name = False
ts.branch_vertical_margin = 10  # 10 pixels between adjacent branches

def my_layout(node):
    if node.name == "N0": # Format the root node with the ancestor sequence
        F = TextFace(tight_text=True, fgcolor="black")
        # F.rotation = 330
        add_face_to_node(F, node, column=0, position="branch-top")

    if node.name in tree:
        print (node.name)

        F = TextFace("Gabe", tight_text=True, fgcolor="red")
        F.rotation = 330
        add_face_to_node(F, node, column=0, position="branch-bottom")


ts.layout_fn = my_layout
tree.render(tree_style=ts, h=320, units="mm")
