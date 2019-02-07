#!/usr/bin/env python
import math



class PileupSite:
	def __init__(self,chr,pos,refc,samples):
		self.chr=chr
		self.pos=pos
		self.refc=refc
		self.samples=samples
		
		

class PileupCountSample:
	def __init__(self,counthash):
		self.__ch=counthash
		self.__cov=sum(counthash.values())
		
	def get_A(self):
		return self.__ch['A']
	def get_a(self):
		return self.__ch['a']
	def get_T(self):
		return self.__ch['T']
	def get_t(self):
		return self.__ch['t']
	def get_C(self):
		return self.__ch['C']
	def get_c(self):
		return self.__ch['c']
	def get_G(self):
		return self.__ch['G']
	def get_g(self):
		return self.__ch['g']	
	
	def get_XXX_lower(self,char):
		lc=char.lower()
		if lc not in self.__ch:
			raise Exception("Invalid character "+char)
		return self.__ch[lc]
	
	def get_XXX_upper(self,char):
		uc=char.upper()
		if uc not in self.__ch:
			raise Exception("Invalid character "+char)
		return self.__ch[uc]
		
	
	def get_forward(self):
		return self.get_A()+self.get_T()+self.get_C()+self.get_G()
	
	def get_reverse(self):
		return self.get_a()+self.get_t()+self.get_c()+self.get_g()
		
	
	def get_forward_freq(self):
		return float(self.get_forward())/float(self.get_cov())
		
	def get_strandbias(self):
		ff=self.get_forward_freq()
		bias=math.fabs(ff-0.5)
		return bias

	def get_XXX(self,char):
		return self.get_XXX_lower(char)+self.get_XXX_upper(char)
		
	def get_sumA(self):
		return self.get_A()+self.get_a()
	def get_sumT(self):
		return self.get_T()+self.get_t()
	def get_sumC(self):
		return self.get_C()+self.get_c()
	def get_sumG(self):
		return self.get_G()+self.get_g()
	def get_cov(self):
		return self.__cov
		
	
	def isPolymorphic(self,mincount):
		sitecount=0
		if self.get_sumA()>=mincount:
			sitecount+=1
		if self.get_sumT()>=mincount:
			sitecount+=1
		if self.get_sumC()>=mincount:
			sitecount+=1
		if self.get_sumG()>=mincount:
			sitecount+=1
		if sitecount>=2:
			return True
		else:
			return False

class StrandSyncWriter:
	def __init__(self,file):
		self.__ofh=open(file,"w")
	
	def write(self,ppcs):
		topr=[]
		topr.append(ppcs.chr)
		topr.append(str(ppcs.pos))
		topr.append(str(ppcs.refc))
		for s in ppcs.samples:
			form="{0}:{1}:{2}:{3}:{4}:{5}:{6}:{7}".format(s.get_A(),s.get_a(),s.get_T(),s.get_t(),s.get_C(),s.get_c(),s.get_G(),s.get_g())
			topr.append(form)
		self.__ofh.write("\t".join(topr)+"\n")
	
	def close(self):
		self.__ofh.close()


class PileupCountSite:
	
	def __init__(self,chr,pos,refc,samples):
		self.chr=chr
		self.pos=pos
		self.refc=refc
		self.samples=samples
	
	def count_samples(self):
		return len(self.samples)
	
	def get_siteMinCoverage(self):
		covs=[]
		[covs.append(i.get_cov()) for i in self.samples]
		return min(covs)

	def get_siteMaxCoverage(self):
		covs=[]
		[covs.append(i.get_cov()) for i in self.samples]
		return max(covs)
		
	
	def is_coveragesufficient(self,mincov):
		for s in self.samples:
			if s.get_cov()<mincov:
				return False
		return True
	
	def get_maxstrandbias(self):
		sb=0.0
		for s in self.samples:
			bias=s.get_strandbias()
			if bias> sb:
				sb=bias
		return sb
	
	def is_strandbiased(self,maxbias):
		for s in self.samples:
			bias=s.get_strandbias()
			if bias>maxbias:
				return True
		return False
	

	def get_ATCGhash(self):
		h={'A':0,'T':0,'C':0,'G':0}
		
		for s in self.samples:
			h['A']+=s.get_sumA()
			h['T']+=s.get_sumT()
			h['C']+=s.get_sumC()
			h['G']+=s.get_sumG()
		return h
	
	def get_maj_min(self):
		h=self.get_ATCGhash()
		i=h.items()
		isor=sorted(i,key=lambda i:-i[1])
		maj=isor[0][0]
		min=isor[1][0]
		return maj, min
	
	def get_XXX(self,char):
		toret=[]
		for s in self.samples:
			toret.append(s.get_XXX(char))
		return toret
		
	
	def isPolymorphic(self,mincount):
		"""
		Is the site polymorphic across all samples;
		"""
		h=self.get_ATCGhash()
		sitecount=0
		if h['A'] >= mincount:
			sitecount+=1
		if h['T'] >= mincount:
			sitecount+=1
		if h['C'] >=mincount:
			sitecount+=1
		if h['G'] >=mincount:
			sitecount+=1
		if sitecount>=2:
			return True
		else:
			return False






class PileupCountReader:
	"""
	A pileup count reader
	Returns PileupCountSite 
	
	"""
	def __init__(self,file,minqual):
		self.__plr=PileupLightwightReader(file)
		self.__offset=33
		self.__mq=minqual


	def __iter__(self):
		return self
	
	def close(self):
		self.__plr.close()


	def next(self):
		tp=self.__plr.next()
		offset=self.__offset
		mq=self.__mq
		
		sms=tp.samples
		newsamples=[]
		for i in sms:
			counth={'A':0,'a':0,'T':0,'t':0,'C':0,'c':0,'G':0,'g':0}
			for k in i:
				char=k[0]
				qualchar=k[1]
				qual=ord(qualchar)-offset
				if(qual<mq):
					continue
				
				if char in counth:
					counth[char]+=1
				else:
					pass
			newsamples.append(PileupCountSample(counth))
		return PileupCountSite(tp.chr,tp.pos,tp.refc,newsamples)
		

class StrandSyncReader:
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
		2R	157	A	A:a:T:t:C:c:G:g		A:a:T:t:C:c:G:g
		2R	158	A	A:a:T:t:C:c:G:g		A:a:T:t:C:c:G:g
		"""
		a=line.split("\t")
		if(len(a)<4):
			raise Exception("invalid sync line "+line)
		chr=a.pop(0)
		pos=int(a.pop(0))
		refc=a.pop(0)
		samples=[]
		
		for b in a:
			A,a,T,t,C,c,G,g=[int(i) for i in b.split(":")]
			st={'A':A,'a':a,'T':T,'t':t,'C':C,'c':c,'G':G,'g':g}
			stsa=PileupCountSample(st)
			samples.append(stsa)
			
		return PileupCountSite(chr,pos,refc,samples)

	





class PileupLightwightReader:
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
		2R	157	A	9	.......^F.^F.	IIIIIIIII
		2R	158	A	9	.........	IIIIIIIII
		"""
		a=line.split("\t")
		if(len(a)%3!=0):
			raise Exception("invalid mpileup line "+line)
		chr=a.pop(0)
		pos=int(a.pop(0))
		refc=a.pop(0)
		samples=[]
		while(len(a)>0):
			cov=a.pop(0)
			toparse=a.pop(0)
			qual=list(a.pop(0))
			
			tuplelist=[]
			if(toparse=="*"):
				pass # add empty list to samples in case of the *
			else:
				parsedchars=self.__parsecharstring(toparse,refc)
				if len(parsedchars)!=len(qual):
					raise Exception("character have different length than the quality string")
				tuplelist=zip(parsedchars,qual) # you'v got to love python, zip is great :)
			samples.append(tuplelist)
	
		return PileupSite(chr,pos,refc,samples)

	def __parsecharstring(self,charstring,refc):
		clist=list(charstring)
		lowrc=refc.lower()
		uprc=refc.upper()
		
		toret=[]
		while(len(clist)>0):
			c=clist.pop(0)
			if c==".":
				toret.append(uprc)
			elif c==",":
				toret.append(lowrc)
			elif c=="^":   #          IGNORE   start of read
				clist.pop(0) 
			elif c=="$":   #          IGNORE   end of read
				pass # 
			elif c=="+" or c=="-":  # IGNORE indels
				# pain in the arse; first get the counts of the indel than remove the sequence of the indel
				# -1G	+1C
				digit=clist.pop(0)
				while(len(clist)>0 and clist[0].isdigit()):
					digit+=clist.pop(0)
				count=int(digit)
				for i in range(0,count):
					clist.pop(0)
			else: # all remeining chars *, >, <, ATCGN, atcgn just append
				toret.append(c)
		return toret
		
	
		
		