import numpy as np
import pandas as pd

from Bio.Align import MultipleSeqAlignment as MSA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import collections  

def get_score_matrix(s1,s2, blosum, gap):
    F = [[0 for i in range(len(s2)+1)] for  j in range(len(s1)+1)]
    gap_status = None

    for i in range(1,len(s1)+1):
        F[i][0] = gap[0] + gap[1]*(i-1)

    for j in range(1,len(s2)+1):
        F[0][j] = gap[0] + gap[1]*(i-1)

    for i in range(1,len(s1)+1):
        for j in range(1,len(s2)+1):
            match = F[i-1][j-1] + blosum.loc[s1[i-1],s2[j-1]]
            if gap_status == 'delete':
                delete = F[i-1][j] + gap[1]
            else:
                delete = F[i-1][j] + gap[0]

            if gap_status == 'insert':
                insert = F[i][j-1] + gap[1]
            else:
                insert = F[i][j-1] + gap[0]

            if match>delete and match>insert:
                F[i][j] = match
                gap_status = None
            elif delete>match and delete>insert:
                F[i][j] = delete
                gap_status = 'delete'
            else:
                F[i][j] = insert
                gap_status = 'insert'
    return F


def blosum_score(c1,c2):
    p1 = dict(collections.Counter(c1))
    p2 = dict(collections.Counter(c2))
    l1 = len(c1)
    l2 = len(c2)

    score = 0
    for r1 in p1:
        for r2 in p2:
            score += blosum.loc[r1,r2] * p1[r1]/l1 * p2[r2]/l2

    return score


def get_score_matrix_msa(msa1,msa2, blosum, gap):
    l1 = msa1.get_alignment_length()
    l2 = msa2.get_alignment_length()

    F = [[(0,'m') for i in range(l2+1)] for j in range(l1+1)]
    gap_status = None

    for i in range(1,l1+1):
        F[i][0] = (gap[0] + gap[1]*(i-1),'d')

    for j in range(1,l2+1):
        F[0][j] = (gap[0] + gap[1]*(j-1),'i')

    for i in range(1,l1+1):
        for j in range(1,l2+1):
            match = F[i-1][j-1][0] + blosum_score(msa1[:,i-1],msa2[:,j-1])
            if gap_status == 'delete':
                delete = F[i-1][j][0] + gap[1]*len(msa2)
            else:
                delete = F[i-1][j][0] + gap[0]*len(msa2)

            if gap_status == 'insert':
                insert = F[i][j-1][0] + gap[1]*len(msa1)
            else:
                insert = F[i][j-1][0] + gap[0]*len(msa1)

            if match>delete and match>insert:
                F[i][j] = (match,'m')
                gap_status = None
            elif delete>match and delete>insert:
                F[i][j] = (delete,'d')
                gap_status = 'delete'
            else:
                F[i][j] = (insert,'i')
                gap_status = 'insert'
            print(i,j)
            print(match, insert, delete)
            print(F[i][j])
            print('--------')
    return F


def align(s1,s2):
    gap = (-10,-4)
    gap_status = None  
    blosum = pd.read_csv('blosum62.txt',sep='\t', header=0, index_col=0)
    F = get_score_matrix(s1, s2, blosum, gap)
    print(F)
    A = ""
    B = ""
    i = len(s1)
    j = len(s2)
    while i>0 or j>0:
        if i>0 and j>0 and F[i][j] == F[i-1][j-1] + blosum.loc[s1[i-1],s2[j-1]]:
            A = s1[i-1] + A
            B = s2[j-1] + B
            i += -1
            j += -1

        elif i>0 and (F[i][j] == F[i-1][j] + gap[0] or F[i][j] == F[i-1][j] + gap[1]):
            A = s1[i-1] + A
            B = '-' + B
            i += -1

        elif j>0 and (F[i][j] == F[i][j-1] + gap[0] or F[i][j] == F[i][j-1] + gap[1]):
            A = '-' + A
            B = s2[j-1] + B
            j += -1

        else:
            print('Fij',F[i][j])
            print('Fi-1j',F[i-1][j])
            print('Fij-1',F[i][j-1  ])
            print('Fi-1j-1',F[i-1][j-1])
            print('match',blosum.loc[s1[i-1],s2[j-1]])
            return None

    return ([A,B],F[len(s1)][len(s2)])


def align_msa(msa1,msa2):
    gap = (-4,-4)
    gap_status = None  
    blosum = pd.read_csv('blosum62.txt',sep='\t', header=0, index_col=0)
    F = get_score_matrix_msa(msa1, msa2, blosum, gap)
    print(F)
    
    l1 = msa1.get_alignment_length()
    l2 = msa2.get_alignment_length()

    A = ["" for a in range(len(msa1))]
    B = ["" for b in range(len(msa2))]

    i = l1
    j = l2

    while i>0 or j>0:
        if i>0 and j>0 and F[i][j] == F[i-1][j-1] + blosum_score(msa1[:,i-1],msa2[:,j-1]):
            A = [msa1[a,i-1] + A[a] for a in range(len(msa1))]
            B = [msa2[b,j-1] + B[b] for b in range(len(msa2))]
            i += -1
            j += -1

        elif i>0 and (F[i][j] == F[i-1][j] + gap[0] or F[i][j] == F[i-1][j] + gap[1]*l2):
            A = [msa1[a,i-1] + A[a] for a in range(len(msa1))]
            B = ['-' + B[b] for b in range(len(msa2))]
            i += -1

        elif j>0 and (F[i][j] == F[i][j-1] + gap[0] or F[i][j] == F[i][j-1] + gap[1]*l1):
            A = ['-' + A[a] for a in range(len(msa1))]
            B = [msa2[b,j-1] + B[b] for b in range(len(msa2))]
            j += -1

        else:
            print(i,j)
            print('Fij',F[i][j])
            print('Fi-1j',F[i-1][j])
            print('Fij-1',F[i][j-1  ])
            print('Fi-1j-1',F[i-1][j-1])
            print(A)
            print(B)
            return None

    return ([A,B],F[l1][l2])


def align_v2(msa1,msa2):
    gap = (-4,-4)
    gap_status = None  
    blosum = pd.read_csv('blosum62.txt',sep='\t', header=0, index_col=0)
    F = get_score_matrix_msa(msa1, msa2, blosum, gap)
    print(F)
    
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

        elif F[i][j][1] == 'd':
            A = [msa1[a,i-1] + A[a] for a in range(len(msa1))]
            B = ['-' + B[b] for b in range(len(msa2))]
            i += -1

        elif F[i][j][1] == 'i':
            A = ['-' + A[a] for a in range(len(msa1))]
            B = [msa2[b,j-1] + B[b] for b in range(len(msa2))]
            j += -1

        else:
            print(i,j)
            print('Fij',F[i][j])
            print('Fi-1j',F[i-1][j])
            print('Fij-1',F[i][j-1  ])
            print('Fi-1j-1',F[i-1][j-1])
            print(A)
            print(B)
            return None

    return ([A,B],F[l1][l2]) 


def show(alignment):
    s1 = alignment[0][0]
    s2 = alignment[0][1]
    score = alignment[1]

    for i in range(int(len(s1)/60)):
        print('sp|P17970|KCNAB_DROME      ',s1[i*60:(i+1)*60])
        print('sp|Q8BZN2|KCNV1_MOUSE      ',s2[i*60:(i+1)*60],'\n\n\n')

    print('sp|P17970|KCNAB_DROME      ',s1[int(len(s1)/60)*60:])
    print('sp|Q8BZN2|KCNV1_MOUSE      ',s2[int(len(s1)/60)*60:],'\n\n\n')



if __name__ == '__main__':

    s1 = SeqRecord(Seq(s1))
    s2 = SeqRecord(Seq(s2))

    s3 = SeqRecord(Seq('SEND'))
    s4 = SeqRecord(Seq('STND'))
    s5 = SeqRecord(Seq('ASENDZ'))
    s6 = SeqRecord(Seq('ASE-DG'))

    gap = (-4,-4)
    blosum = pd.read_csv('blosum62.txt',sep='\t', header=0, index_col=0)

    alignment = align_v2(MSA([s3,s4]),MSA([s5,s6]))
    print(alignment)

# TODO: visulisation des MSA
#       lecture des sequences par fichier
#        