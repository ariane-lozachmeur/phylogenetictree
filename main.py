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

from Bio import AlignIO
from Bio.Cluster import treecluster
import ete3
from Bio import Phylo
import tree as perso
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment

import time


def init_clusters(sequences):
    return {i : perso.Tree(sequences[i].id, MultipleSeqAlignment([sequences[i]])) for i in range(len(sequences))}



def blosum_score(seq1, seq2, gap_creation = -4, gap_extension = -2):
    blosum = pd.read_csv('blosum62.txt',sep='\t', header=0, index_col=0)
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


def compute_distances(clusters, update=None, dist=None, mode='seqs'):
    if dist is not None and update is not None:
        i,j,new_id = update
        matrix = dist.drop([i,j],axis=0)
        matrix = matrix.drop([i,j],axis=1)

        matrix[new_id] = np.zeros(len(matrix))
        matrix = matrix.append(pd.DataFrame([np.zeros(len(matrix)+1)],columns=matrix.columns.values,index=[new_id]))

        for k in clusters.keys():
            if not new_id == k:
                d = clusters[new_id].distance(clusters[k], mode=mode)
                matrix.loc[new_id,k] = d
                matrix.loc[k,new_id] = d
            else:
                matrix.loc[k,k] = -math.inf

    else:
        matrix = pd.DataFrame([[-math.inf for col in range(len(clusters))] for row in range(len(clusters))])
        for i in clusters.keys():
            for j in clusters.keys():
                if j > i:
                    d = clusters[i].distance(clusters[j], mode=mode)
                    matrix.loc[i,j] =  d
                    matrix.loc[j,i] = d  
                elif j == i:
                    matrix.loc[i,i] = -math.inf
    return matrix

  
def show(tree, level):
    ret = "\t"*level + str(tree.id)+'\n'
    if tree.left:
        ret += show(tree.left,level+1)
    if tree.right:
        ret += show(tree.right,level+1)

    return ret


def toNewick(tree):
    if tree.depth()==0:
        return str(tree.id)
    else:
        return '(\n'+toNewick(tree.left)+',\n'+toNewick(tree.right)+'\n)\n'


def toTreePerso(tree, seqs):
    leaves = []
    for s in seqs:
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

def main(file, mode, type, dist):
    name = file.split('/')[-1].split('.')[0]+'_'+type
    if mode == 'seqs':
        seqs = list(SeqIO.parse(file, "fasta",alphabet = Alphabet.Gapped(IUPAC.protein)))

    elif mode == 'msa':
        align = AlignIO.read(file,'fasta')
        seqs = align[:,:]
    else:
        raise ValueError('The mode must be either seqs if the input is a list of unaligned sequences or msa if the input is an MSA')

    if type == 'single':
        clusters = init_clusters(seqs)
        dist = compute_distances(clusters, dist = dist, mode=mode)
        tree = treecluster(distancematrix=dist, method='s')
        tree = toTreePerso(tree, seqs)


    elif type == 'average':
        clusters = init_clusters(seqs)
        dist = compute_distances(clusters,dist = dist, mode=mode)
        tree = treecluster(distancematrix=dist, method='a')
        tree = toTreePerso(tree, seqs)

    elif type == 'centroid':
        # print("Initialisation...")
        clusters = init_clusters(seqs)
        # print("Clusters created")
        dist = compute_distances(clusters)
        # print("Distances initialized")
        cluster_id = -1
        # print("Initialisation over")
        while not len(clusters) == 1:
            # print("There are currently %s clusters" %len(clusters))
            # select the minimum of distances
            mins = dist.idxmin()
            i1,i2 = -1,-1
            m = +math.inf
            for i in mins.index:
                if dist.loc[i,mins[i]]<m:
                    i1 = i
                    i2 = mins[i]
                    m = dist.loc[i,mins[i]]

            # define the centroid of the new cluster from the previous clusters
            new_tree = clusters[i1].centroid(clusters[i2],i=cluster_id, mode=mode)

            # remove the clustered we joined of the cluster list 
            del clusters[i1]
            del clusters[i2]

            clusters[cluster_id] = new_tree

            # Update the distances
            dist = compute_distances(clusters, update=(i1,i2,cluster_id), dist=dist, mode=mode)
            cluster_id += -1 
            tree = clusters[list(clusters.keys())[0]]


    else:
        raise ValueError('The type of clustering must be either single, averzage or centroid')

    newick = toNewick(tree) + ';'
    show_tree(newick, name=name)


if __name__ == "__main__":

    for t in ['single','average','centroid']:
        for file in ['KCNB2','LRRD1','TRAF6']:
            print(t,file)
            start = time.time()
            main("test/"+file+".fa", mode='seqs', type=t, dist='blosum')
            end = time.time()
            print(end-start)
