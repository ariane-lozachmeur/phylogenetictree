from ete3 import Tree
import pandas as pd
from Bio import Phylo

# File to benchmark our methods with MUSCLE and MAFFT using the Robinson Foulds metric.
# It create a benchmark_results.csv file in the "results" repository and prints all of the tree that can then be saved for further analysis.

prots = ['LRRD1','TRAF6','KCNB2']
methods = ['single_tree','average_tree','centroid_tree','muscle','mafft','centroidmsa_tree']


def main(prots, methods):
    results = pd.DataFrame( )
    for p in prots:
        try:
            t1 = Tree(open('results/'+p+'_uniprot.dnd').read())
        except:
            print('WARNING: There is no uniprot file for the protein %s' %p)
            continue
        print(p,'uniprot')
        # Phylo.draw(Phylo.read('results/'+p+'_uniprot.dnd', "newick"))
        scores = []
        for m in methods:
            print(p,m)
            try:
                t2 = Tree(open('results/'+p+'_'+m+'.dnd').read())
            except:
                print('WARNING: There is no %s file for the protein %s' %(m,p))
                continue
            # Phylo.draw(Phylo.read('results/'+p+'_'+m+'.dnd', "newick"))
            rf, max_rf, common_leaves, parts_t1, parts_t2, i, j = t1.robinson_foulds(t2, unrooted_trees=True)
            scores.append(float(rf)/float(max_rf))

        results = results.append(pd.Series(scores,name=p))

    results.columns = methods
    results.to_csv('results/benchmark_results.csv')

if __name__ == '__main__':
    main(prots, methods)