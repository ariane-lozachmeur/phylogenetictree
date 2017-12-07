import time
import math
import os

import numpy as np
import pandas as pd

from Bio import SeqIO, pairwise2 
from Bio.Align import MultipleSeqAlignment as MSA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat.MatrixInfo import blosum62

import collections  

from Bio import AlignIO



def score(c1,c2, matrix):
    p1 = dict(collections.Counter(c1))
    p2 = dict(collections.Counter(c2))
    l1 = len(c1)
    l2 = len(c2)

    score = 0
    for r1 in p1:
        for r2 in p2:
            score += blosum.loc[r1,r2] * p1[r1]/l1 * p2[r2]/l2

    return score


def score_matrix(msa1,msa2, matrix, gap, type):
    l1 = msa1.get_alignment_length()
    l2 = msa2.get_alignment_length()

    F = [[({'match':0, 'Iy':0, 'Ix':0},'m') for i in range(l2+1)] for j in range(l1+1)]
    gap_status = None

    for i in range(1,l1+1):
        F[i][0] = ({'match':-math.inf, 'Iy':gap[0] + gap[1]*(i-1), 'Ix':-math.inf},'Iy')

    for j in range(1,l2+1):
        F[0][j] = ({'match':-math.inf, 'Iy':-math.inf, 'Ix':gap[0] + gap[1]*(j-1)},'Ix')

    for i in range(1,l1+1):
        for j in range(1,l2+1):
            match = max(F[i-1][j-1][0].values()) + score(msa1[:,i-1],msa2[:,j-1],matrix=matrix)
            Iy = max( F[i-1][j][0]['match'] + gap[0], F[i-1][j][0]['Iy'] + gap[1])
            Ix = max( F[i][j-1][0]['match'] + gap[0], F[i][j-1][0]['Ix'] + gap[1])
            
            if match>Iy and match>Ix:
                state = 'm'
            elif Iy>match and Iy>Ix:
                state = 'Iy'
            else:
                state = 'Ix'

            F[i][j] = ({'match':match, 'Iy':Iy, 'Ix':Ix},state)

    for i in range(l1+1):
        for j in range(l2+1):
            F[i][j] = (max(F[i][j][0].values()),F[i][j][1])
    
    return F


def needleman_wunsch(msa1,msa2, gap):
    if type(msa1) is SeqRecord:
        msa1 = MSA([msa1])
        msa2 = MSA([msa2])
    elif type(msa1) is Seq:
        msa1 = MSA([SeqRecord(msa1)])
        msa2 = MSA([SeqRecord(msa2)])

    gap_status = None  
    blosum = pd.read_csv('blosum62.txt',sep='\t', header=0, index_col=0)
    F = score_matrix(msa1, msa2, matrix = blosum, gap = gap, type='needleman')
    # print(F)
    
    l1 = msa1.get_alignment_length()
    l2 = msa2.get_alignment_length()

    A = ["" for a in range(len(msa1))]
    B = ["" for b in range(len(msa2))]

    i = l1
    j = l2

    while i>0 or j>0:
        if F[i][j][1] == 'm':
            A = [msa1[a,i-1] + A[a] for a in range(len(msa1))]
            B = [msa2[b,j-1] + B[b] for b in range(len(msa2))]
            i += -1
            j += -1

        elif F[i][j][1] == 'Iy':
            A = [msa1[a,i-1] + A[a] for a in range(len(msa1))]
            B = ['-' + B[b] for b in range(len(msa2))]
            i += -1

        elif F[i][j][1] == 'Ix':
            A = ['-' + A[a] for a in range(len(msa1))]
            B = [msa2[b,j-1] + B[b] for b in range(len(msa2))]
            j += -1

        else:
            print(F[i][j][1])
            return None
    msa = MSA([SeqRecord(Seq(A[a]),id=msa1[a].id) for a in range(len(A))])
    B = MSA([SeqRecord(Seq(B[b]),id=msa2[b].id) for b in range(len(B))])

    msa.extend(B)
    return (msa,F[l1][l2][0]) 


def show(align1, align2, title=None):
    if title:
        print('Alignment : ' + title)

    l1 = align1.get_alignment_length()
    l2 = align2.get_alignment_length()
    l = max(l1,l2)
    i = 0
    for i in range(int(l/60)):
        for s in align1:
            print(s.id+'_1      ',s.seq[i*60:(i+1)*60])
        print()

        for s in align2:
            print(s.id+'_2      ',s.seq[i*60:(i+1)*60])
        print('\n')

    for s in align1:
        print(s.id+'_1      ',s.seq[(i+1)*60:])
    print('\n')

    for s in align2:
        print(s.id+'_2      ',s.seq[(i+1)*60:])
    print('\n')
    print('************************************************************')


    
def compare_align(align1,align2):
    align1.sort()
    align2.sort()

    identity = 0
    l1 = align1.get_alignment_length()
    l2 = align2.get_alignment_length()

    pos_1 = [0 for a in range(len(align1))]
    pos_2 = [0 for b in range(len(align2))]
    
    columns_1 = []
    columns_2 = []

    for i in range(l1):
        for a in range(len(align1)):
            if not align1[a][i] == '-':
                pos_1[a] += 1
        columns_1.append([pos_1[k] for k in range(len(pos_1))])

    for j in range(l2):
        for b in range(len(align2)):
            if not align2[b][j] == '-':
                pos_2[b] += 1
        columns_2.append([pos_2[k] for k in range(len(pos_2))])
        
    TP = 0
    FP = 0

    for c1 in columns_1:
        if c1 in columns_2:
            TP += 1
            columns_2.remove(c1) 
        else:
            FP += 1

    FN = len(columns_2)
    recall = float(TP) / float(TP + FN)
    precision = float(TP) / float(TP + FP)
    show(align1, align2)
    print('FN=',FN, 'FP=', FP, 'TP=', TP)
    print('Recall = %s ; precision = %s' %(recall,precision))
    return (recall,precision)


def identity(s1,s2):
    ident = 0
    for i in range(min(len(s1.seq),len(s2.seq))):
        if s1[i] == s2[i]:
            ident+=1
    return float(ident)/max(len(s1.seq),len(s2.seq))


if __name__ == '__main__':

    gap = (-10,-0.5)
    blosum = pd.read_csv('blosum62.txt',sep='\t', header=0, index_col=0)
    seqs = []
    # print(needleman_wunsch(Seq('SEND'),Seq('AND'),gap=gap)[0])

    results = []
    for filename in os.listdir('test/prefab4/in'):
        print(filename)
        seqs = []
        for record in SeqIO.parse('test/prefab4/in/'+filename, "fasta"):
            if record.id in filename:
                seqs.append(record)
        
        for i in range(1,len(seqs)):
            # print('Aligning sequence 0 with sequence %s' %i)
            start = time.time()
            s = [seqs[0], seqs[i]]
            alignment = needleman_wunsch(s[0],s[1], gap)
            end = time.time()

            real = AlignIO.read('test/prefab4/ref/'+filename, "fasta")
            # for s in real:
            #     s.seq = Seq(str(s.seq).replace('.','-').upper())
            # start2 = time.time()
            # recall,precision = compare_align(alignment[0],real)    
            # end2 = time.time()
            # time_align = end-start
            
            # results.append([filename,recall,precision,time_align])
            results.append([filename,alignment[1]])
            # print('Time alignment:', end-start)
            # print('Time comparison:', end2-start2)
    df_results = pd.DataFrame(results,columns=['seq','score'])
    df = pd.read_csv('results/alignment_benchmark2.csv')
    df = df.merge(df_results,on='seq',how='inner')
    df.to_csv('results/alignment_benchmark3.csv',index=False)