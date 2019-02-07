#!/usr/bin/env python


class FastqPoolWriter:
	"""
	Write the content to a fasta file
	"""
	def __init__(self,prefix):
		self.__prefix=prefix
		self.__fh1=open(prefix+"_1.fastq","w")
		self.__fh2=open(prefix+"_2.fastq","w")
		self.__counter=1

	def write(self,s1,s2,sampleid):
		fh1=self.__fh1
		fh2=self.__fh2
		fh1.write("@"+str(self.__counter)+"\n")
		fh1.write(s1+"\n")
		fh1.write("+"+str(self.__counter)+"\n")
		fh1.write("I"*len(s1)+"\n")
		
		fh2.write("@"+str(self.__counter)+"\n")
		fh2.write(s2+"\n")
		fh2.write("+"+str(self.__counter)+"\n")
		fh2.write("I"*len(s2)+"\n")

		self.__counter+=1

	def close(self):
		self.__fh1.close()
		self.__fh2.close()
		
		


class FastqIndividualWriter:
	"""
	Write the content to a fasta file
	"""
	def __init__(self,prefix):
		self.__prefix=prefix
		self.__sampleid=0
		self.__batchcounter=0
		self.__fh1=None
		self.__fh2=None
	def write(self,s1,s2,sampleid):
		if(self.__sampleid!=sampleid):
			self.__sampleid=sampleid
			if(sampleid%2==1):
				self.__batchcounter+=1
				if(self.__fh1 is not None):
					self.__fh1.close()
				if(self.__fh2 is not None):
					self.__fh2.close()
				self.__fh1=open(self.__prefix+str(self.__batchcounter)+"_1.fastq")
				self.__fh2=open(self.__prefix+str(self.__batchcounter)+"_2.fastq")
		fh1=self.__fh1
		fh2=self.__fh2
		fh1.write("@"+str(self.__counter)+"\n")
		fh1.write(s1+"\n")
		fh1.write("+"+str(self.__counter)+"\n")
		fh1.write("I"*len(s1)+"\n")
		
		fh2.write("@"+str(self.__counter)+"\n")
		fh2.write(s2+"\n")
		fh2.write("+"+str(self.__counter)+"\n")
		fh2.write("I"*len(s2)+"\n")

		self.__counter+=1

	def close(self):
		self.__fh1.close()
		self.__fh2.close()