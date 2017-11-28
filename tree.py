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
        self.left = left[0]
        self.dleft = left[1]
        self.right = right[0]
        self.dright = right[1]


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

      # def distance(self, tree):
        # if str(self.id)<str(tree.id):
        #     l1 = len(self.msa)
        # else:
        #     l1 = len(tree.msa)
        # score = 0
        # seqs = self.msa.format("fasta") + tree.msa.format("fasta")
        
        # with open("tmp.fa","w") as f:
        #     f.write(seqs)
        
        # muscle_cline = MuscleCommandline(input="tmp.fa")
        # child = subprocess.Popen(str(muscle_cline),stdout=subprocess.PIPE,stderr=subprocess.PIPE, universal_newlines=True, shell=(sys.platform!="win32"))
        # align = AlignIO.read(child.stdout, "fasta")
        # align.sort()

        # for s1 in align[:l1]:
        #     for s2 in align[l1:]:
        #         score += blosum_score(s1.seq,s2.seq)
        # dist = float(score)/float(len(self.msa)*len(tree.msa))

        
                
        # os.remove('tmp.fa')
        # return dist

    def distance(self, tree, mode):
        d = 0

        if str(self.id)<str(tree.id):
            self_is_first = True
        else:
            self_is_first = False

        for s in tree.msa:
            d += self.score_seq(s, self_is_first, mode=mode)

        for s in self.msa:
            d += tree.score_seq(s, not self_is_first, mode=mode)

        dist = float(d)/float(len(self.msa)*len(tree.msa))

        return -dist


    def score_seq(self, s, self_is_first, mode):
        d = 0
        seqs = self.msa.format("fasta") + s.format("fasta")

        if mode == 'seqs':
            # If the input was unaligned sequences, we need to align them using MUSCLE.
            with open("tmp.fa","w") as f:
                f.write(seqs)

            muscle_cline = MuscleCommandline(input="tmp.fa")
            child = subprocess.Popen(str(muscle_cline),stdout=subprocess.PIPE,stderr=subprocess.PIPE, universal_newlines=True, shell=(sys.platform!="win32"))
            align = AlignIO.read(child.stdout, "fasta")
            align.sort()

        elif mode == 'msa':
            # Else we use the alignment of the initial MSA.
            align = self.msa[:,:]
            align.append(s)

            self_is_first = True
        else:
            raise ValueError('The mode must be either seqs if the input is a list of unaligned sequences or msa if the input is an MSA')

        if self_is_first:
            s2 = align[-1]
            d = blosum_score(align[:len(self.msa)],s2.seq)
        else:
            s2 = align[0]
            d = blosum_score(align[1:],s2.seq)

        return d


    def centroid(self,tree,i,mode):
        for s in self.msa:
            s.id = str(i) + s.id

        for s in tree.msa:
            s.id = str(i) + s.id
        seqs = self.msa.format("fasta") + tree.msa.format("fasta")
       
        if mode == 'seqs':
        # If the input was unaligned sequences, we need to align them using MUSCLE.
            with open("tmp.fa","w") as f:
                f.write(seqs)
            
            muscle_cline = MuscleCommandline(input="tmp.fa")
            child = subprocess.Popen(str(muscle_cline),stdout=subprocess.PIPE,stderr=subprocess.PIPE, universal_newlines=True, shell=(sys.platform!="win32"))
            align = AlignIO.read(child.stdout, "fasta")
        elif mode == 'msa':
            align = self.msa[:,:]
            align.extend(tree.msa)
        else:
            raise ValueError('The mode must be either seqs if the input is a list of unaligned sequences or msa if the input is an MSA')
        
        return Tree(i, align, left=self, right=tree)

                

def blosum_score(msa, seq, gap_creation = -4, gap_extension = -2):
    blosum = pd.read_csv('blosum62.txt',sep='\t', header=0, index_col=0)
    score = 0
    gapped = False
    profile = []
    for i in range(msa.get_alignment_length()):
        profile.append(dict(collections.Counter(msa[:,i])))

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


