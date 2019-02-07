#!/usr/bin/env python

class TEInsertionSite:
     def __init__(self,teid,plusstrand):
          self.teid=teid
          self.plusstrand=plusstrand
     def getKey(self):
          strand=""
          if(self.teid>0):
            if(self.plusstrand):
              strand="+"
            else:
              strand="-"
          return str(self.teid)+strand
     
     @classmethod
     def parse_key(cls,key):
        teid=None
        plusStrand=True
        if key.endswith(("+","-")):
            teid=int(key[:-1])
            strandstring=key[-1:]
            if(strandstring=="-"):
                plusStrand=False
        else:
            teid=int(key)
        return TEInsertionSite(teid,plusStrand)	
        

class TEInsertionDefinition:
	def __init__(self,position,comment,tesites):
		self.position=position
		self.comment=comment
		top=[]
		tmph={}
		if "," in comment:
			top=comment.split(",")
		else:
			top.append(comment)
		
		for t in top:
			if "=" not in t:
				continue
			a=t.split("=")
			tmph[a[0]]=a[1]
		
		self.__commenthash=tmph
		self.tesites=tesites
		
	def getSitecount(self):
		return len(self.tesites)
	
	def getTSD(self):
		if "tsd" in self.__commenthash:
			return int(self.__commenthash["tsd"])
		else:
			return 0

class TEDefinitionReader:
	"""
	A light-weight reader of TE definitions;
	returns a TEInsertionDefinition
	
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
		while(1):
			line=self.__filehandle.readline()
			if line=="":
				raise StopIteration
			line=line.rstrip('\n')
			if line.startswith("#"):
				continue
			if line != "":
				break
		a=line.split(" ")
		pos=int(a.pop(0))
		comment=a.pop(0)
		tedefs=[]
		for tmp in a:
			site=TEInsertionSite.parse_key(tmp)
			tedefs.append(site)	
		# position, comment
		e = TEInsertionDefinition(pos,comment,tedefs)
		return e
	
	@classmethod
	def readall(cls,file):
		entries=[]
		for e in TEDefinitionReader(file):
			entries.append(e)
		return entries

class TEDefinitionWriter:
	"""
	Write the content of a TE definition to a file
	"""
	def __init__(self,file):
		self.__filename=file
		self.__filehandle=open(file,"w")
		
	def write(self,tedef):
		position,comment,tedefinitions=tedef.position,tedef.comment,tedef.tesites
		topr=[str(position),comment]
		for ted in tedefinitions:
			topr.append(ted.getKey())
		form=" ".join(topr)
		self.__filehandle.write(form+"\n")
		return 1
	
	def close(self):
		self.__filehandle.close()
	
	@classmethod
	def writeall(cls,file,tedefentries):
		gw=TEDefinitionWriter(file)
		for e in tedefentries:
			gw.write(e)
		gw.close()
		return 1