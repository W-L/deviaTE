#!/usr/bin/env python
import random
import math
	

class TEHierarchy:
	def __init__(self,id2fam,id2ord):
		self.__id2fam=id2fam
		self.__id2ord=id2ord
		
	def getFam(self,id):
		return self.__id2fam[id]
	
	def getOrd(self,id):
		return self.__id2ord[id]
		

class TEHierarchyDefault:
	
	def __init__(self):
		pass
	
	def getFam(self,id):
		return id
	
	def getOrd(self,id):
		return id
		

def loadtehier(hierfile):
	id2fam={}
	id2ord={}
	
	firstlinetoread=True
	idcol=0
	famcol=0
	ordcol=0
	for line in open(hierfile):
		a=line.rstrip("\n").split("\t")
		if firstlinetoread:
			firstlinetoread=False
			for i,t in enumerate(a):
				if(t=="id"):
					idcol=i
				elif t=="family":
					famcol=i
				elif t=="order":
					ordcol=i
		else:
			id=a[idcol]
			fam=a[famcol]
			ord=a[ordcol]
			id2fam[id]=fam
			id2ord[id]=ord
	return TEHierarchy(id2fam,id2ord)
			
	
