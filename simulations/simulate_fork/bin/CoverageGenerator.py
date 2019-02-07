#!/usr/bin/env python
import random











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
        """
        zero based index of the genome
        """
        return self.__reads[index]
            
    
    