import sys
import collections
import subprocess
from copy import copy

import numpy as np
import pandas as pd

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline

print('Load cmd line')
muscle_cline = MuscleCommandline(input="tmp_error.fa")
print('Build child')
child = subprocess.Popen(str(muscle_cline),stdout=subprocess.PIPE,stderr=subprocess.PIPE, universal_newlines=True, shell=(sys.platform!="win32"))
print('Read align')
print(child.stdout)
data = child.stdout.read()
print(data)
align = AlignIO.read(child.stdout, "fasta")
print('End read align')