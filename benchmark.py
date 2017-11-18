from ete3 import Tree

target = 'LRRD1_HUMAN'

t1 = Tree(open('results/'+target+'_uniprot.dnd').read())
t2 = Tree(open('results/'+target+'_tree.dnd').read())
t3 = Tree(open('results/'+target+'_mafft.dnd').read())

rf, max_rf, common_leaves, parts_t1, parts_t2, i, j = t1.robinson_foulds(t2, unrooted_trees=True)

print("Perso: RF distance is %s over a total of %s" %(rf, max_rf))

rf, max_rf, common_leaves, parts_t1, parts_t2, i, j = t1.robinson_foulds(t3, unrooted_trees=True)

print("MAFFT: RF distance is %s over a total of %s" %(rf, max_rf))

# We can also compare trees sharing only part of their labels

# t1 = Tree('(((a,b),c), ((e, f), g));')
# t2 = Tree('(((a,c),b), (g, H));')
# rf, max_rf, common_leaves, parts_t1, parts_t2 = t1.robinson_foulds(t2)

# print(t1, t2)
# print("Same distance holds even for partially overlapping trees")
# print("RF distance is %s over a total of %s" %(rf, max_rf))
# print("Partitions in tree2 that were not found in tree1:", parts_t1 - parts_t2)
# print("Partitions in tree1 that were not found in tree2:", parts_t2 - parts_t1)