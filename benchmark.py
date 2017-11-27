from ete3 import Tree

target = 'MYD88_HUMAN'

t1 = Tree(open('results/'+target+'_uniprot.dnd').read())
t2 = Tree(open('results/'+target+'_tree.dnd').read())
t3 = Tree(open('results/'+target+'_mafft.dnd').read())

rf, max_rf, common_leaves, parts_t1, parts_t2, i, j = t1.robinson_foulds(t2, unrooted_trees=True)

print("Perso: RF distance is %s over a total of %s" %(rf, max_rf))
print("Number of branches in tree2 and not in tree1: ", len(parts_t1 - parts_t2))
print("Number of branches in tree1 and not in tree2: ", len(parts_t2 - parts_t1))

print(i)
print(j)

f, max_rf, common_leaves, parts_t1, parts_t2, i, j = t1.robinson_foulds(t3, unrooted_trees=True)

print("MAFFT: RF distance is %s over a total of %s" %(rf, max_rf))
print("Number of branches in tree2 and not in tree1: ", len(parts_t1 - parts_t2))

