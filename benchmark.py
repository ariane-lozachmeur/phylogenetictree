from ete3 import Tree
import pandas as pd
from Bio import Phylo

def show_tree(newick, name=None):
    if name:
        filename = name+'_tree.dnd'
    else:
        filename = 'tree.dnd'
    handle = open('results/'+filename,'w')
    handle.write(newick)
    handle.close()

    tree = Phylo.read('results/'+filename, "newick")
    # Phylo.draw_ascii(tree)
    Phylo.draw(tree)

prots = ['LRRD1','TRAF6','KCNB2']
methods = ['single_tree','average_tree','centroid_tree','muscle','mafft']
results = pd.DataFrame( )
for p in prots:
    t1 = Tree(open('results/'+p+'_uniprot.dnd').read())
    print(p,'uniprot')
    Phylo.draw(Phylo.read('results/'+p+'_uniprot.dnd', "newick"))
    scores = []
    for m in methods:
        print(p,m)
        t2 = Tree(open('results/'+p+'_'+m+'.dnd').read())
        Phylo.draw(Phylo.read('results/'+p+'_'+m+'.dnd', "newick"))
        rf, max_rf, common_leaves, parts_t1, parts_t2, i, j = t1.robinson_foulds(t2, unrooted_trees=True)
        scores.append(float(rf)/float(max_rf))

    results = results.append(pd.Series(scores,name=p))

results.columns = methods
print(results)
results.to_csv('results/benchmark_results.csv')

