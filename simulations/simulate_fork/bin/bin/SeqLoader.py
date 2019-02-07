
#!/user/bin/env python
import os
import sys
import re
import time
import datetime
import random
from fastaIO import FastaReader, FastaWriter

def rc(seq):
    return ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in seq][::-1])

def load_chasis(inputFile):
    fr=FastaReader(inputFile)
    counter=0
    h=None
    s=None
    
    for header,seq in fr:
        h=header
        s=seq
        counter+=1
    if counter>1:
        raise ValueError("Chasis file must only contain a single reference genome")
    return h,s




        

def getPopGenomeStats(inputFile):
    totleng=0
    totcount=0
    for head,seq in FastaReader(inputFile):
        totcount+=1
        totleng+=len(seq)
    averagegenomesize=float(totleng)/float(totcount)
    popgenomesamples=totcount
    
    return averagegenomesize,popgenomesamples
