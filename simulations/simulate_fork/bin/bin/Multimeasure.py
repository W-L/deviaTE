#!/usr/bin/env python

class Multimeasure:
    
    def __init__(self,chr,pos,issnp,isoutlier,strandbias,measure):
        self.chr=chr
        self.pos=pos
        self.issnp=issnp
        self.isoutlier=isoutlier
        self.strandbias=strandbias
        self.measure=measure


class MultimeasureReader:
	"""
	A light-weight pileup reader
	
	"""
	def __init__(self,file):
		self.__filename=file
		self.__filehandle=open(file,"r")

	def __iter__(self):
		return self
	
	def close(self):
		self.__filehandle.close()


	def next(self):
		line=""
		while(1):
			line=self.__filehandle.readline()
			if line=="":
				raise StopIteration
			line=line.rstrip('\n')
			if line.startswith("#"):
				continue
			if line != "":
				break
		
		e=self.__parseLine(line)
		return e

	def __parseLine(self,line):
		"""
		(self,chr,pos,issnp,isoutlier,strandbias,measure):
		"""
		a=line.split("\t")
		if(len(a)!=6):
			raise Exception("invalid multimeasure line "+line)
		chr=a[0]
		pos=int(a[1])
		issnp=False
		if(a[2]=="True"):
		    issnp=True
		isoutlier=False
		if(a[3]=="True"):
		    isoutlier=True
		strandbias=float(a[4])
		measure=float(a[5])
		return Multimeasure(chr,pos,issnp,isoutlier,strandbias,measure)

class MultimeasureWriter:
	def __init__(self,file):
		self.__ofh=open(file,"w")
	
	def write(self,mm):
		# (self,chr,pos,issnp,isoutlier,strandbias,measure)
		topr=[mm.chr, mm.pos, mm.issnp, mm.isoutlier,mm.strandbias, mm.measure]
		topr=[str(i) for i in topr]
		self.__ofh.write("\t".join(topr)+"\n")
	
	def close(self):
		self.__ofh.close()