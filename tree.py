import sys
import collections
import subprocess
from copy import copy

import numpy as np
import pandas as pd

from io import StringIO

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline


# The class Tree allows the progressive construction of the tree. Each tree is a cluster.
class Tree:

    # msa is an MultipleSeqAlignment containing all the sequences of the cluster
    # Left and right are the 2 clusters/trees merged to create the current tree.
    def __init__(self, i, msa=None, left=None, right=None):
        self.id = i
        self.msa = msa
        self.left = left
        self.right = right


    def set_children(self, left, right):
        self.left = left[0]
        self.dleft = left[1]
        self.right = right[0]
        self.dright = right[1]


    # The depth of the tree represents the number of children it has 
    def depth(self):
        if self.left == None and self.right == None:
            return 0

        else:
            return 1 + max(self.left.depth(),self.right.depth())

    # Creates the centroid between two trees.
    def centroid(self,tree,i,type):
        for s in self.msa:
            s.id = str(i) + s.id

        for s in tree.msa:
            s.id = str(i) + s.id
        seqs = self.msa.format("fasta") + tree.msa.format("fasta")
       
        # If the input was unaligned sequences, we need to align them using MUSCLE.
        if type == 'seqs':
            with open("tmp.fa","w") as f:
                f.write(seqs)
            
            # print('Load cmd line')
            # muscle_cline = MuscleCommandline(input="tmp.fa")
            # print('Build child')
            # child = subprocess.Popen(str(muscle_cline),stdout=subprocess.PIPE,stderr=subprocess.PIPE, universal_newlines=True, shell=(sys.platform!="win32"))
            # print('Read align')
            # print(child.stdout)
            # align = AlignIO.read(child.stdout, "fasta")
            # print('End read align')

            muscle_cline = MuscleCommandline(input="tmp.fa")
            stdout, stderr = muscle_cline()
            align = AlignIO.read(StringIO(stdout), "fasta")
            print(align)
        
        # Else we use the alignment of the initial MSA
        elif type == 'msa':
            align = self.msa[:,:]
            align.extend(tree.msa)
        else:
            raise ValueError('The type must be either seqs if the input is a list of unaligned sequences or msa if the input is an MSA')
        
        # Return the centroid as a tree
        return Tree(i, align, left=self, right=tree)


    # Computes the distance between the two trees
    def distance(self, tree, type):
        score = 0

        # The self_is_first variable is used to determined which MSA's sequences will be first in the new sorted MSA
        # To make sure the sequences stay grouped in the new MSA we use the Tree's id
        if str(self.id)<str(tree.id):
            self_is_first = True
        else:
            self_is_first = False

        # We score each sequence of the second MSA against the first MSA
        for s in tree.msa:
            score += self.score_seq(s, self_is_first, type=type)

        # We score each sequence of the first MSA against the second MSA
        for s in self.msa:
            score += tree.score_seq(s, not self_is_first, type=type)

        # The distance is the opposite of the score (to have a small distance of close MSA). 
        # We also divide by the number of sequences in the MSA to prevent the large clusters to have larger distances.
        dist = -float(score)/float(len(self.msa)*len(tree.msa))

        return dist


    def score_seq(self, s, self_is_first, type):
        score = 0
        seqs = self.msa.format("fasta") + s.format("fasta")

        # If the input was unaligned sequences, we need to align them using MUSCLE.
        if type == 'seqs':
            with open("tmp.fa","w") as f:
                f.write(seqs)

            muscle_cline = MuscleCommandline(input="tmp.fa")
            child = subprocess.Popen(str(muscle_cline),stdout=subprocess.PIPE,stderr=subprocess.PIPE, universal_newlines=True, shell=(sys.platform!="win32"))
            align = AlignIO.read(child.stdout, "fasta")
            align.sort()
        
        # Else we use the alignment of the initial MSA.
        elif type == 'msa':
            align = self.msa[:,:]
            align.append(s)
            self_is_first = True

        else:
            raise ValueError('The type must be either seqs if the input is a list of unaligned sequences or msa if the input is an MSA')

        # After aligning sequences we score the alignment using the blosum_score function
        if self_is_first:
            s2 = align[-1]
            score = blosum_score(align[:len(self.msa)],s2.seq)
        else:
            s2 = align[0]
            score = blosum_score(align[1:],s2.seq)

        return score
                
# Computes the score of the alignment between a sequence and an MSA (described in scoring section of the report)
def blosum_score(msa, seq, gap_creation = -4):
    # Load the BL62 matrix
    blosum = pd.read_csv('blosum62.txt',sep='\t', header=0, index_col=0)
    score = 0
    profile = []
    # Creation of the MSA's profile
    for i in range(msa.get_alignment_length()):
        profile.append(dict(collections.Counter(msa[:,i])))

    # Scoring according to the formula in the report
    for j in range(msa.get_alignment_length()):
        for residue in profile[j]:
            score+= profile[j][residue]/len(msa) * blosum.loc[residue,seq[j]]

    return score   


if __name__=='__main__':
    a = SeqRecord(Seq("ARCC-EA"), id="Alpha")
    b = SeqRecord(Seq("ARDNQEA"), id="Beta")
    c = SeqRecord(Seq("ARNTEA"), id="Celta")
    d = SeqRecord(Seq("ATNE-A"), id="Delta")
    align1 = MultipleSeqAlignment([a, b])
    align2 = MultipleSeqAlignment([c, d])

    tree1 = Tree(0, align1)
    tree2 = Tree(0, align2)

    print(tree1.distance(tree2))


