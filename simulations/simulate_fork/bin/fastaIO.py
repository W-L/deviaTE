#!/usr/bin/env python
import sys
import random
import collections

class SequenceUtility:
	
	@classmethod
	def rc(cls,sequence):
		# complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N' }
		complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'Y':'R', 'R':'Y', 'W':'W', 'S':'S',
					  'K':'M', 'M':'K', 'D':'H', 'H':'D', 'V':'B', 'B':'V'}
		bases = reversed(list(sequence)) 
		bases = [complement[base] for base in bases] 
		return ''.join(bases)
	
	@classmethod
	def load_chasis(cls,fastafile):
		seqtups=FastaReader.readAllTuples(fastafile)
		if len(seqtups)>1:
			raise Exception("The chasis is a single fasta entry; to avoid confusion the chasis-file must contain only a single fasta entry!")
		return seqtups[0] # first entry and than the second position within tuple (= sequence, the first position is the header)
	
	@classmethod
	def get_length_list(cls,inputFile):
		ll=[]
		fr=FastaReader(inputFile)
		for head,seq in fr:
			ll.append(len(seq))
		fr.close()
		return ll


class FastaBatchWriter:
	"""
	Batch writer for individual PE sequencing; 
	"""
	def __init__(self,prefix,haploid,seqleng):
		self.__prefix=prefix
		self.__sampleid=0
		self.__batchcounter=0
		self.__haploid=haploid
		self.__seqleng=seqleng
		self.__writer=None
		

	def write(self,header, seq1,sampleid):
		if(self.__sampleid!=sampleid):
				self.__sampleid=sampleid
				if(sampleid%2==1 or self.__haploid):
						self.__batchcounter+=1
						if(self.__writer is not None):
								self.__writer.close()
						fn1=self.__prefix+str(self.__batchcounter)+".fasta"
						print("Writing to fasta-file {0}".format(fn1))
						self.__writer=FastaWriter(fn1,self.__seqleng)
		self.__writer.write(header,seq1)

	def close(self):
		self.__writer.close()


class FastaReader:
	"""
	A light-weight fasta reader;
	returns a tuple (header, sequence)
	
	"""
	def __init__(self,file):
		self.__filename=file
		self.__filehandle=open(file,"r")
		self.__prevheader=None

	def __iter__(self):
		return self
	
	def close(self):
		self.__filehandle.close()
	
	def next(self):
		line=""
		header=self.__prevheader
		seq=""
		while(1):
			line=self.__filehandle.readline()
			if line=="":					# file is empty
				if(header is not None):
					self.__prevheader=None		# return last sequence
					return (header,seq)
				else:
					raise StopIteration		# file empty and no last sequence -> STOP
			line=line.rstrip("\n")				# somethin is in the file
			if(line.startswith(">")):			# we have a header
				line=line.lstrip(">")
				if(header is None):			# if it is the first header just set the name of the header
					header=line
				else:
					self.__prevheader=line	# if it is any other header, set the previous to the current and return the sequence
					return(header,seq)
			else:
				seq+=line				# normal line, add to sequence
	
	@classmethod
	def readAllHash(cls,file):
		fh={}
		for n,s in FastaReader(file):
			if n in fh:
				raise ValueError("Invalid sequence IDs. Sequence with ID "+n+" is occuring multiple times")
			fh[n]=s
		return fh
	
	@classmethod
	def readAllTuples(cls,file):
		ft=[]
		for n,s in FastaReader(file):
			ft.append((n,s))
		return ft

class FastaWriter:
	"""
	Write the content to a fasta file
	"""
	def __init__(self,file,seqleng):
		self.__filename=file
		self.__filehandle=open(file,"w")
		self.__seqleng=seqleng
		
	def write(self,n,s):
		sl=self.__seqleng
		fh=self.__filehandle
		fh.write(">"+n+"\n")
		c=0
		while(c<len(s)):
			fh.write(s[c:c+sl]+"\n")
			c+=sl

	def close(self):
		self.__filehandle.close()
	
	@classmethod
	def write_all(cls,file,fastaentries):
		fw=FastaWriter(outputFile)
		for n,s in fastaentries:
			fw.write(n,s)
		fw.close()
		return 1
	
