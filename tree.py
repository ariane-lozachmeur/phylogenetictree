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
import collections
from copy import copy
import sys

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

import subprocess
from Bio.Align.Applications import MuscleCommandline

import matplotlib.pyplot as plt 
import os

class Tree:
    def __init__(self, i, msa=None, left=None, right=None):
        self.id = i
        self.msa = msa
        self.left = left
        self.right = right


    def set_children(self, left, right):
        self.left = left
        self.right = right


    def depth(self):
        if self.left == None and self.right == None:
            return 0

        else:
            return 1 + max(self.left.depth(),self.right.depth())


    def get_profile(self, seqs = None):
        profile = []
        if seqs == None:
            for i in range(self.msa.get_alignment_length()):
                profile.append(dict(collections.Counter(self.msa[:,i])))

        else:
            for i in range(self.msa.get_alignment_length()):
                profile.append(dict(collections.Counter(self.msa[seqs[0]:seqs[1],i])))            
        return profile        


    def distance(self, tree):
        l1 = len(self.msa)
        l2 = len(tree.msa)
        score = 0
        seqs = self.msa.format("fasta") + tree.msa.format("fasta")
        with open("tmp.fa","w") as f:
            f.write(seqs)
        
        muscle_cline = MuscleCommandline(input="tmp.fa")
        child = subprocess.Popen(str(muscle_cline),stdout=subprocess.PIPE,stderr=subprocess.PIPE, universal_newlines=True, shell=(sys.platform!="win32"))
        align = AlignIO.read(child.stdout, "fasta")
        align.sort()
        for s1 in align[:l1]:
            for s2 in align[l1:]:
                score += blosum_score(s1.seq,s2.seq)
        dist = float(score)/float(l1*l2)
        os.remove('tmp.fa')
        return dist


    def centroid(self,tree,i):
        for s in self.msa:
            s.id = str(i) + s.id

        for s in tree.msa:
            s.id = str(i) + s.id
        seqs = self.msa.format("fasta") + tree.msa.format("fasta")
       
        with open("tmp.fa","w") as f:
            f.write(seqs)
        
        muscle_cline = MuscleCommandline(input="tmp.fa")
        child = subprocess.Popen(str(muscle_cline),stdout=subprocess.PIPE,stderr=subprocess.PIPE, universal_newlines=True, shell=(sys.platform!="win32"))
        align = AlignIO.read(child.stdout, "fasta")
        
        return Tree(i, align, left=self, right=tree)

                

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


