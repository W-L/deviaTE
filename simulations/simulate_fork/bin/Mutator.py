#!/usr/bin/env python
import random
import math
	
class PacBioMutator:

	def __init__(self,errorrate,delfrac):
		self.__er=errorrate
		self.__tr={'A':['T','C','G'], 'T':['A','C','G'],'C':['A','T','G'],'G':['A','T','C']}
		self.__ins=['A','T','C','G']
		self.__delfrac=delfrac;
	
	def mutateseq(self,seq):
		if(self.__er<0.0000000001):
			return seq
			
		muts=[]
		for i in range(0,len(seq)):
			if random.random()<self.__er:
				muts.append(i)
		
		ins=self.__ins
		lseq=list(seq)
		
		for i in reversed(muts):
			if random.random()<self.__delfrac:
				#deletion
				del(lseq[i])
		
			else:
				#insertion
				random.shuffle(ins)
				lseq.insert(i,ins[0])
		
		return "".join(lseq)






class PoisonSeqMutator:

	def __init__(self,errorrate):
		self.__er=errorrate
		self.__tr={'A':['T','C','G'], 'T':['A','C','G'],'C':['A','T','G'],'G':['A','T','C']}
	
	def __getPoisson(self,lam):
		L= math.exp(-lam)
		p=1.0
		k=0
		while(True):
			k+=1
			p*=random.random()
			if(p<L):
				break
		return k-1
	
	def mutateseq(self,seq):
		if(self.__er<0.0000000001):
			return seq

		lseq=list(seq)
		for i in range(0,len(lseq)):
			# test every base if it should have an error
			if random.random()<self.__er:
				tr=self.__tr[lseq[i]]
				random.shuffle(tr)
				lseq[i]=tr[0]
		return "".join(lseq)
	
	
	
	
	

	def mutateseq(self,seq):
		if(self.__er<0.0000000001):
			return seq
		aver=float(len(seq))*self.__er
		errors=self.__getPoisson(aver)
		lseq=list(seq) # sequence in a list
		pastpos=set([])  
		for i in range(0,errors):
			pos=random.randint(0,len(lseq)-1)
			if pos in pastpos:
				pos=random.randint(0,len(lseq)-1)
			pastpos.add(pos)
			tr=self.__tr[lseq[pos]]
			random.shuffle(tr)
			lseq[pos]=tr[0]
		return "".join(lseq)


class ExhaustiveSeqMutator:
	
	def __init__(self,errorrate):
		self.__er=errorrate
		self.__tr={'A':['T','C','G'], 'T':['A','C','G'],'C':['A','T','G'],'G':['A','T','C']}
	

	def mutateseq(self,seq):
		if(self.__er<0.0000000001):
			return seq

		lseq=list(seq)
		for i in range(0,len(lseq)):

			# # -- do not mutate base 499 - 501
			# if i in range(499, 501):
			# 	continue
			# # -----
			
			# -- do not mutate base 1999 - 2001
			if i in range(1999, 2001):
				continue
			# -----
			
			# test every base if it should have an error
			if random.random()<self.__er:
				tr=self.__tr[lseq[i]]
				random.shuffle(tr)
				lseq[i]=tr[0]
		return "".join(lseq)
		
