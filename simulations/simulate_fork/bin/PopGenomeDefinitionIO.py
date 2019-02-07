#!/usr/bin/env python
import sys
import random
import collections
import re


class PopGenDefinitionReader:
	"""
	
	"""
	def __init__(self,file,sc):
		self.__filename=file
		self.__sc=sc
		self.__chasis=""
		self.insertions=None
		self.popsize=None
		self.tuples=self.read_tuples()
	
	def __parse_chasis(self,line):
		a=line.split("=")
		r=a[1].strip(" \"")
		self.__chasis=r
		
	def get_chasis(self):
		return self.__chasis
	
	def read_tuples(self):
			toret=[]
			teleng=None
			sc=self.__sc
			for l in open(self.__filename):
				l= l.rstrip("\n")
				l=re.sub(r"#(.*)","",l) # remove comments
				l=l.strip()
				if l=="": # remove empty lines
					continue
				if(l.startswith("chassis")):
					self.__parse_chasis(l)
				elif "=" in l:
					sc.addDefinition(l)
				else:
					sp=l.split(" ")
					if teleng is None:
						teleng=len(sp)
					if len(sp)!=teleng:
						raise Exception("Invalid popgendef (population genome definition) file; unequal numbers of haploid genomes among insertion sites: "+len(sp)+ " "+teleng)
					toret.append(sp)
			self.insertions=len(toret)
			self.popsize=teleng-1 # minus one since the first column is the position within the chasis
			return toret
		
	def read_transposed(self):
			tuples=self.tuples
			popsize=0
			sc=self.__sc
			if len(tuples)<1:
				raise Exception("Invalid popgendef (population genome definition) file; at least one TE insertion site is required")
			popsize=len(tuples[0])-1
			toret=[]
			for i in range(1,popsize+1):
				tmp=[]
				for k in range(0,len(tuples)):
					toa=(int(tuples[k][0]),tuples[k][i]) # position, tedefinition
					if toa[1]=="*":
						continue
					tmp.append(toa)
				toret.append(tmp)
			return toret
				
				
	


class PopGenDefinitionWriter:
	"""
	Write the content to a fasta file
	"""
	def __init__(self,file,popsize):
		self.__filename=file
		self.__filehandle=open(file,"w")
		self.__popsize=popsize
		
		
	def write_chasis_info(self,header,seq):
		fh=self.__filehandle
		fh.write("# Chasis {0}; Length {1} nt\n".format(header,len(seq)))
		
		
	def write_header(self,tetuples):
		fh=self.__filehandle
		for i,tt in enumerate(tetuples):
			name,seq=tt
			k=i+1
			fh.write("{0}=${0}     # {1}\n".format(k,name))
		return True
	
	def write_popfreq(self,pos,entry,freq):
		"""
		write an entry with the given position and population frequency
		"""
		n=self.__popsize
		fh=self.__filehandle
		count_te=int(freq*n)
		count_empty=n-count_te
		toshuf=[entry for i in range(0,count_te)]
		toshuf.extend("*"*count_empty)
		random.shuffle(toshuf)
		toshuf.insert(0,str(pos))
		tows=" ".join(toshuf)
		fh.write(tows+"\n")
		
	def write_empty(self,pos):
		n=self.__popsize
		fh=self.__filehandle
		tow=[str(pos)]
		tow.extend("*"*n)
		tows=" ".join(tow)
		fh.write(tows+"\n")


	def close(self):
		self.__filehandle.close()
	

	
