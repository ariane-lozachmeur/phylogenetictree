import numpy as np
import math
import pandas as pd
import Bio
from Bio import pairwise2
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.Cluster import treecluster
import matplotlib

from Bio.Cluster import treecluster
import ete3
from Bio import Phylo
import tree as perso
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment


# Input: unaligned sequences or MSA
# - Compute initial matrix of distances (A)
# - Place each sequence in its own cluster
# - While # clusters>1
# o Join the 2 closest clusters (using the matrix of distances)
# o Compute the new centroid representing this new cluster (B)
# o Add this step to the tree (creation of the new link between 2 clusters)
# o Update the matrix of distances
# Output: phylogenetic tree

def init_clusters(sequences=None, MSA=None):
    if sequences:
        print('Starting init clusters')
        return {i : perso.Tree(sequences[i].id, MultipleSeqAlignment([sequences[i]])) for i in range(len(sequences))}
    else:
        return None


def init_distances(sequences=None, MSA=None, distance='blast'):
    matrix = np.matrix([[0 for col in range(len(sequences))] for row in range(len(sequences))])
    if sequences:
        for i in range(len(sequences)):
            s1 = sequences[i]
            for j in range(i+1,len(sequences)):
                s2 = sequences[j]
                alpha = Alphabet.Gapped(IUPAC.protein)
                alignment = pairwise2.align.globalds(s1.seq, s2.seq, blosum62, -10, -0.5)
                
                if distance=='blast':
                    dist = -alignment[0][2]
                
                elif distance=='blosum':
                    seq1 = alignment[0][0]
                    seq2 = alignment[0][1]
                    dist = -blosum_score(seq1,seq2)

                matrix.itemset((i, j), dist)
                matrix.itemset((j,i),dist)

    if MSA:
        pass

    return matrix


def blosum_score(seq1, seq2, gap_creation = -4, gap_extension = -2):
    blosum = pd.read_csv('blosum62.qij',sep='\t', header=0, index_col=0)
    score = 0
    gapped = False
    for i in range(len(seq1)):
        if seq1[i] == '-':
            if gapped == 1:
                score += gap_extension
            else:
                score += gap_creation
                gapped = 1
        elif seq2[i] == '-':
            if gapped == 2:
                score += gap_extension
            else:
                score += gap_creation
                gapped = 2

        else:
            gapped = False
            score += blosum.loc[seq1[i],seq2[i]]

    return score


def compute_distances(clusters, update=None, dist=None):
    if dist is not None and update is not None:
        i,j,new_id = update
        print(i,j)
        matrix = dist.drop([i,j],axis=0)
        matrix = matrix.drop([i,j],axis=1)

        matrix[new_id] = np.zeros(len(matrix))
        matrix = matrix.append(pd.DataFrame([np.zeros(len(matrix)+1)],columns=matrix.columns.values,index=[new_id]))

        for k in clusters.keys():
            if not new_id == k:
                d = clusters[new_id].distance(clusters[k])
                matrix.loc[new_id,k] = d
                matrix.loc[k,new_id] = d

    else:
        matrix = pd.DataFrame([[0 for col in range(len(clusters))] for row in range(len(clusters))])
        for i in clusters.keys():
            for j in clusters.keys():
                if j > i:
                    d = clusters[i].distance(clusters[j])
                    matrix.loc[i,j] =  d
                    matrix.loc[j,i] = d  
    return matrix

  
def show(tree, level):
    ret = "\t"*level + str(tree.id)+'\n'
    if tree.left:
        ret += show(tree.left,level+1)
    if tree.right:
        ret += show(tree.right,level+1)

    return ret


def toNewick(clusters):
    tree = clusters[clusters.keys()[0]]
    if tree.depth()==0:
        return str(tree.id)
    else:
        return '('+toNewick(tree.left)+','+toNewick(tree.right)+')'


def toTreePerso(tree):
    leaves = []
    for s in globins:
        leaves.append(perso.Tree(s.id))

    nodes = {}
    n = -1
    for node in tree:
        if node.left>=0:
            left = leaves[node.left]

        else:
            left = nodes[node.left]

        if node.right>=0:
            right = leaves[node.right]

        else:
            right = nodes[node.right]

        new_node = perso.Tree(n, msa=None, left=left, right=right)
        nodes[n]=new_node
        n+=-1

    return nodes[min(nodes.keys())]


def show_tree(newick):
    handle = open('tree.dnd','w')
    handle.write(newick)
    handle.close()

    tree = Phylo.read("tree.dnd", "newick")
    Phylo.draw_ascii(tree)
    Phylo.draw(tree)

def main(seqs, type, dist):
    
    if type == 'single':
        dist = init_distances(seqs,distance = dist)
        tree = treecluster(distancematrix=dist, method='s')
        tree = toTreePerso(tree)

    elif type == 'average':
        dist = init_distances(seqs,distance = dist)
        tree = treecluster(distancematrix=dist, method='a')
        tree = toTreePerso(tree)

    elif type == 'centroid':
        clusters = init_clusters(seqs)
        dist = compute_distances(clusters)
        cluster_id = -1
        while not len(clusters) == 1:
            # select the minimum of distances
            print(dist)
            maxs = dist.idxmax()
            i1,i2 = -1,-1
            m = 0
            for i in maxs.index:
                if dist.loc[i,maxs[i]]>m:
                    i1 = i
                    i2 = maxs[i]
                    m = dist.loc[i,maxs[i]]

            # define the centroid of the new cluster from the previous clusters
            new_tree = clusters[i1].centroid(clusters[i2],i=cluster_id)

            # remove the clustered we joined of the cluster list 
            del clusters[i1]
            del clusters[i2]

            clusters[cluster_id] = new_tree

            # Update the distances
            dist = compute_distances(clusters, update=(i1,i2,cluster_id), dist=dist)
            cluster_id += -1        

    newick = toNewick(clusters)
    show_tree(newick)


if __name__ == "__main__":

    globins = list(SeqIO.parse("globins.fa", "fasta",alphabet = Alphabet.Gapped(IUPAC.protein)))
    main(globins, type="centroid", dist='blosum')
    
    
