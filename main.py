import math
import time

import numpy as np
import pandas as pd

import Bio
from Bio import SeqIO, AlignIO, Phylo, Alphabet
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from Bio.Cluster import treecluster

import tree as perso


# The program will loop through all the methods, proteins and type of inputs.
# The input file must be formatted as such: protein.input (for example, KCNB2.msa or TRAF6.fa)

# Choose among centroid, single and average
methods = ['centroid']
# Name of the dataset
proteins = ['TLR']
# Choose among fa (for a list of sequences) and msa
inputs = ['fa']


def init_clusters(sequences):
    return {i : perso.Tree(sequences[i].id, MultipleSeqAlignment([sequences[i]])) for i in range(len(sequences))}


# Creates or updates the matrix of distances between clusters (detailled in step A of the report)
def compute_distances(clusters, update=None, dist=None, type='seqs'):
    
    # If it is NOT the initialisation
    if dist is not None and update is not None:    
        # We delete from the distance matrix the columns and rows corresponding to the deleted trees
        i,j,new_id = update
        matrix = dist.drop([i,j],axis=0)
        matrix = matrix.drop([i,j],axis=1)

        # We add a new column
        matrix[new_id] = np.zeros(len(matrix))
        matrix = matrix.append(pd.DataFrame([np.zeros(len(matrix)+1)],columns=matrix.columns.values,index=[new_id]))

        # We compute the new distances
        for k in clusters.keys():
            if not new_id == k:
                # The distance function is in the tree.py file
                d = clusters[new_id].distance(clusters[k], type=type)
                matrix.loc[new_id,k] = d
                matrix.loc[k,new_id] = d
            else:
                matrix.loc[k,k] = -math.inf

    # If it is the initialisation we compute ALL the distances
    else:
        matrix = pd.DataFrame([[-math.inf for col in range(len(clusters))] for row in range(len(clusters))])
        for i in clusters.keys():
            for j in clusters.keys():
                if j > i:
                    d = clusters[i].distance(clusters[j], type=type)
                    matrix.loc[i,j] =  d
                    matrix.loc[j,i] = d  
                elif j == i:
                    matrix.loc[i,i] = -math.inf
    return matrix


# Converts recursively the Tree object in a string following Newick's format
def to_newick(tree):
    if tree.depth()==0:
        return str(tree.id)
    else:
        return '(\n'+to_newick(tree.left)+',\n'+to_newick(tree.right)+'\n)\n'


# Converts the Tree object returned by Bio.Cluster.treecluster into a custom Tree object (for the single and average likage methods)
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


# Displays the newick tree 
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


# Main function of the algorithm. Type represents the type of input and type the mode of clustering
def main(file, type, mode, dist):
    name = file.split('/')[-1].split('.')[0]+'_'+mode+type
    if type == 'seqs':
        # Create a Tree for each input sequence 
        seqs = list(SeqIO.parse(file, "fasta",alphabet = Alphabet.Gapped(IUPAC.protein)))

    elif type == 'msa':
        # Create a Tree for each input sequence 
        align = AlignIO.read(file,'fasta')
        seqs = align[:,:]
    else:
        raise ValueError('The type must be either seqs if the input is a list of unaligned sequences or msa if the input is an MSA')

    # The single and average modes are using a built-in function of Biopython and are not detailled in the report
    if mode == 'single':
        clusters = init_clusters(seqs)
        dist = compute_distances(clusters, dist = dist, type=type)
        tree = treecluster(distancematrix=dist, method='s')
        tree = toTreePerso(tree, seqs)


    elif mode == 'average':
        clusters = init_clusters(seqs)
        dist = compute_distances(clusters,dist = dist, type=type)
        tree = treecluster(distancematrix=dist, method='a')
        tree = toTreePerso(tree, seqs)

    # Centroid method detailled in the report
    elif mode == 'centroid':
        
        # Initialisation step
        print("Initialisation...")
        clusters = init_clusters(seqs)
        print("Clusters created")
        dist = compute_distances(clusters)
        print("Distances initialized")
        cluster_id = -1
        print("Initialisation over")
        
        # Processing step
        while not len(clusters) == 1:
            print("There are currently %s clusters" %len(clusters))
            # select the minimum of distances
            mins = dist.idxmin()
            i1,i2 = -1,-1
            m = +math.inf
            for i in mins.index:
                if dist.loc[i,mins[i]]<m:
                    i1 = i
                    i2 = mins[i]
                    m = dist.loc[i,mins[i]]

            # create the new tree t3, centroid of t1 and t2
            new_tree = clusters[i1].centroid(clusters[i2],i=cluster_id, type=type)

            # remove the clustered we joined of the cluster list 
            del clusters[i1]
            del clusters[i2]

            # add t3 to the list of clusters
            clusters[cluster_id] = new_tree

            # Update the distances
            dist = compute_distances(clusters, update=(i1,i2,cluster_id), dist=dist, type=type)
            cluster_id += -1 
            tree = clusters[list(clusters.keys())[0]]


    else:
        raise ValueError('The mode of clustering must be either single, averzage or centroid')

    # Termination step
    # Store the tree in Newick format and display it.
    newick = to_newick(tree) + ';'
    show_tree(newick, name=name)


if __name__ == "__main__":

    for m in methods:
        for file in proteins:
            for i in inputs:
                print("Executing the following tree:")
                print("Methods: ",m)
                print("Dataset: ",file)
                print("Type of input: ",i)
                start = time.time()
                if i == 'fa':
                    t = 'seqs'
                else:
                    t = 'msa'
                main("test/"+file+'.'+i, type=t, mode=m, dist='blosum')
                end = time.time()
                print("Time taken: ",end-start)
