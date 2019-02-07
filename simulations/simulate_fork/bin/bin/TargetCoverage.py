#!/usr/bin/env python
import random



def get_uniform_readnumbergenerator(targetcoveragestring,insertsize,readleng,pgll):
    if targetcoveragestring.endswith("rpg"):
        reads=int(targetcoveragestring.rstrip("rpg"))
        return ConstantReadsPerGenome(reads)
    elif targetcoveragestring.endswith("cpg"):
        targetcov=float(targetcoveragestring.rstrip("cpg"))
        return ConstantCoveragePerGenome(pgll,targetcov,readleng)
    elif targetcoveragestring.endswith("pcpg"):
        targetcov=float(targetcoveragestring.rstrip("pcpg"))
        return ConstantPhysicalCoveragePerGenome(pgll,targetcov,insertsize)
    else:
        raise Exception("Unknown definition "+targetcoveragestring)
        
    # either 2.0pc average physical coverage per genome, 2.0c average coverage per genome, 1000r reads per genome
    return 0;


class ConstantPhysicalCoveragePerGenome:
    def __init__(self,pgll,targetcov,insertsize):
        readtable=[]
        for len in pgll:
            reads = int(float(targetcov*len)/float(insertsize))
            readtable.append(reads)
        self.__reads=readtable
        
    def get_reads(self,index):
        return self.__reads[index]


class ConstantCoveragePerGenome:
    def __init__(self,pgll,targetcov,readleng):
        readtable=[]
        for len in pgll:
            reads = int(float(targetcov*len)/float(2*readleng))
            readtable.append(reads)
        self.__reads=readtable
    
    def get_reads(self,index):
        return self.__reads[index]
            


class ConstantReadsPerGenome:
    def __init__(self,readnum):
        self.readnum=readnum
        
    def get_reads(self,index):
        return self.readnum
    
############################################################################################################################################
############################################################ RANDOM ########################################################################
############################################################################################################################################

def get_random_readnumbergenerator(targetcoveragestring,insertsize ,readleng, pgll):
    if targetcoveragestring.endswith("r"):
        reads=int(targetcoveragestring.rstrip("r"))
        return RandomReads(reads,pgll)
    elif targetcoveragestring.endswith("c"):
        targetcov=float(targetcoveragestring.rstrip("c"))
        return RandomCoverage(pgll,targetcov,readleng)
    elif targetcoveragestring.endswith("pcpg"):
        targetcov=float(targetcoveragestring.rstrip("pcpg"))
        return RandomPhysicalCoveragePerGenome(pgll,targetcov,insertsize)
    else:
        raise Exception("Unknown definition "+targetcoveragestring)
        
    # either 2.0pc average physical coverage per genome, 2.0c average coverage per genome, 1000r reads per genome
    return 0;


class RandomReads:
    def __init__(self,readnum,pgll):
        self.__reads=RandomReads.getreadtable(readnum,pgll)
        
    
    @classmethod
    def getreadtable(cls,reads,pgll):
        table=[0,]*len(pgll)
        for i in range(0,reads):
            r =random.randint(0,len(pgll)-1)
            table[r]+=1
        return table
        
    def get_reads(self,index):
        return self.__reads[index]
        

class RandomCoverage:
    def __init__(self,pgll,targetcov,readleng):
        avgen=float(sum(pgll))/float(len(pgll))
        reads=int(float(avgen*targetcov)/float(2*readlen))
        self.__reads=RandomReads.getreadtable(reads,pgll)
    def get_reads(self,index):
        return self.__reads[index]
        


class RandomPhysicalCoveragePerGenome:
    def __init__(self,pgll,targetcov,insertsize):
        reads=0
        for le in pgll:
            reads+=int(float(le*targetcov)/float(insertsize))
        self.__reads=RandomReads.getreadtable(reads,pgll)

    def get_reads(self,index):
        return self.__reads[index]
        
    
    
    