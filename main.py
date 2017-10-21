import numpy as np
import math
import pandas as pd
import Bio
from Bio import pairwise2
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat.MatrixInfo import blosum62

# Input: unaligned sequences or MSA
# - Compute initial matrix of distances (A)
# - Place each sequence in its own cluster
# - While # clusters>1
# o Join the 2 closest clusters (using the matrix of distances)
# o Compute the new centroid representing this new cluster (B)
# o Add this step to the tree (creation of the new link between 2 clusters)
# o Update the matrix of distances
# Output: phylogenetic tree

def init_distances(sequences=None,MSA=None):
    matrix = np.matrix([[math.inf for col in range(len(sequences))] for row in range(len(sequences))])
    if sequences:
        for i in range(len(sequences)):
            s1 = sequences[i]
            for j in range(i+1,len(sequences)):
                s2 = sequences[j]

                alignment = pairwise2.align.globalds(s1.seq, s2.seq, blosum62, -10, -0.5)
                print(pairwise2.format_alignment(*alignment[0]))
                score = alignment[0][2]
                matrix.itemset((i, j), score)
        print(matrix)

    if MSA:
        pass

    return matrix

def update_distance(clusters):
    pass


def centroid():
    pass
    

seqs = [
        SeqRecord(Seq('MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFYR')),
        SeqRecord(Seq('MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFYH')),
        SeqRecord(Seq('MVHLTPEEKSAVTALWGKVNVDEVGGAALGRLLVVYPWTQRFYH')),
        ]

if __name__ == "__main__":
    dist = init_distances(seqs)
    clusters = seqs
    while  len(clusters)>=1:
        argmin = np.where(dist == np.min(dist)) 
        i1 = argmin[0][0]
        i2 = argmin[1][0]

        new_centroid = centroid(clusters[i1],clusters[i2])
        del seqs[i1]
        del seqs[i2]
        seqs.append(new_centroid)

        new_cluster = [clusters[i1],clusters[i2]]
        del clusters[i1]
        del clusters[i2-1]
        clusters.append(new_cluster)

        print(clusters)
        dist = update_distances(seqs)