#!/usr/bin/env python
import sys
import random
import collections



class FastqSEBatchWriter:
	"""
	Batch writer for individual PE sequencing; 
	"""
	def __init__(self,prefix,haploid):
		self.__prefix=prefix
		self.__sampleid=0
		self.__batchcounter=0
		self.__haploid=haploid
		self.__writer=None
		

	def write(self,header, seq1,sampleid):
		if(self.__sampleid!=sampleid):
				self.__sampleid=sampleid
				if(sampleid%2==1 or self.__haploid):
						self.__batchcounter+=1
						if(self.__writer is not None):
								self.__writer.close()
						fn1=self.__prefix+str(self.__batchcounter)+".fastq"
						print("Writing to fastq-file {0}".format(fn1))
						self.__writer=FastqWriter(fn1)
		self.__writer.write(header,seq1)

	def close(self):
		self.__writer.close()




class FastqPEBatchWriter:
	"""
	Batch writer for individual PE sequencing; 
	"""
	def __init__(self,prefix,haploid):
		self.__prefix=prefix
		self.__sampleid=0
		self.__batchcounter=0
		self.__haploid=haploid
		self.__writer=None
		

	def write(self,header, seq1,seq2,sampleid):
		if(self.__sampleid!=sampleid):
				self.__sampleid=sampleid
				if(sampleid%2==1 or self.__haploid):
						self.__batchcounter+=1
						if(self.__writer is not None):
								self.__writer.close()
						fn1=self.__prefix+str(self.__batchcounter)+"_1.fastq"
						fn2=self.__prefix+str(self.__batchcounter)+"_2.fastq"
						print("Writing to fastq-files {0} and {1}".format(fn1,fn2))
						self.__writer=FastqPairWriter(fn1,fn2)
		self.__writer.write(header,seq1,seq2)

	def close(self):
		self.__writer.close()


class FastqPairWriter:
	def __init__(self,file1,file2):
		self.__fqw1=FastqWriter(file1)
		self.__fqw2=FastqWriter(file2)
	
	def write(self,header,seq1,seq2):
		self.__fqw1.write(header,seq1)
		self.__fqw2.write(header,seq2)
	
	def close(self):
		self.__fqw1.close()
		self.__fqw2.close()
		



	
class FastqWriter:
	"""
	Write the content to a fastq file
	"""
	def __init__(self,file):
		self.__filename=file
		self.__filehandle=open(file,"w")
		
	def write(self,header,seq):
		fh=self.__filehandle
		fh.write("@"+header+"\n")
		fh.write(seq+"\n")
		fh.write("+"+header+"\n")
		fh.write("I"*len(seq)+"\n")

	def close(self):
		self.__filehandle.close()
	

	
