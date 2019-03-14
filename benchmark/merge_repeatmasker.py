from gtfIO import GTFReader, GTFEntry
from optparse import OptionParser
import collections

parser = OptionParser()
parser.add_option("--input", dest="gtf", help="A gtf file containing the RepeatMasked gtf annotation")
parser.add_option("--output", dest="output", help="An output file for the improved annotation ")
parser.add_option("--match-score", dest="matchscore", help="Match score", default=1.0)
parser.add_option("--mismatch-penalty", dest="mmpenalty", help="The mismatch penalty", default=0.5)
(options, args) = parser.parse_args()


class GTFTEReader:
    def __init__(self, file):
        self.__gtfreader = GTFReader(file)

    def __iter__(self):
        return self

    def next(self):
        for e in self.__gtfreader:
            # 	2L	RepeatMasker	similarity	1724	1849	 6.4	-	.	Target "Motif:DMTRDNA" 1311 1435
            a = e.comment.split(" ")
            target = a[1]
            score = float(e.end - e.start)
            e.sore = score
            target = target.replace('"Motif:', "")
            target = target.rstrip('"')
            e.target = target
            return e
        raise StopIteration


def load_gtfhash(file):
    chrh = collections.defaultdict(lambda: [])
    for e in GTFTEReader(file):
        chrh[e.chr].append(e)
    return chrh


class GTFTEClusterReader:
    def __init__(self, chrh, clusterdistance):
        self.chrh = chrh
        self.clusterdistance = clusterdistance
        self.__listofcluster = self.__createListOfCluster(chrh)

    def __iter__(self):
        return self

    def __createListOfCluster(self, chrh):
        lifc = []
        clud = self.clusterdistance
        for chr, gtfs in self.chrh.items():
            sgtfs = sorted(gtfs, key=lambda g: g.start)

            while(len(sgtfs) > 0):
                active = sgtfs.pop(0)
                highend = active.end
                newclu = [active]
                while(len(sgtfs) > 0 and sgtfs[0].start - clud < highend):
                    act = sgtfs.pop(0)
                    if(act.end > highend):
                        highend = act.end
                    newclu.append(act)
                lifc.append(newclu)
        return lifc

    def next(self):
        while(len(self.__listofcluster) > 0):
            return self.__listofcluster.pop(0)
        raise StopIteration


def load_tefamhash_enablestrand(tecluster):
    tefh = collections.defaultdict(lambda: [])
    for te in tecluster:
        key = te.target + ":" + te.strand
        tefh[key].append(te)
    return tefh


def load_tefamhash_disablestrand(tecluster):
    tefh = collections.defaultdict(lambda: [])
    for te in tecluster:
        key = te.target
        tefh[key].append(te)
    return tefh


def mergeTEentries(entries, matchscore):
    """
    merge overlapping TE entries of the same family.
    New  score is the total length
    """
    merged = []
    tes = sorted(entries, key=lambda e: e.start)
    while(len(tes) > 0):
        a = tes.pop(0)
        start = a.start
        highestend = a.end
        while(len(tes) > 0 and tes[0].start <= highestend):
            b = tes.pop(0)
            if(b.end > highestend):
                highestend = b.end
        ne = GTFEntry(a.chr, a.source, a.feature, start, highestend, matchscore * float(highestend - start), a.strand, a.frame, a.comment)
        ne.target = a.target
        merged.append(ne)
    return merged


def clusterTEentries(entries, mmpenalty):
    """
    cluster non-overlapping TE insertions of the same family.
    a scoring system is used to decide whether non-overlapping TE insertions will be clustered
    """
    clustered = []
    tes = sorted(entries, key=lambda e: e.start)
    while(len(tes) > 0):
        a = tes.pop(0)
        score = a.score
        end = a.end
        while(len(tes) > 0):
            totest = tes[0]
            gap = totest.start - end
            gappen = gap * mmpenalty
            scorewithgap = score - gappen
            if scorewithgap < 0:  # if the score with the gap reaches zero; break
                break
            # novel high score:
            totestscore = scorewithgap + totest.score
            if(totestscore >= score):
                tes.pop(0)
                score = totestscore
                end = totest.end
            else:
                # no new high score
                break

        ne = GTFEntry(a.chr, a.source, a.feature, a.start, end, score, a.strand, a.frame, a.comment)
        ne.target = a.target
        clustered.append(ne)
    return clustered


def iterativeClusterTEentries(entries, mmpenalty):
    clulen = len(entries)
    clustered = clusterTEentries(entries, mmpenalty)
    while(1):
        if(len(clustered) == clulen):
            break
        clulen = len(clustered)
        clustered = clusterTEentries(clustered, mmpenalty)
        # CLUSTERING is finished ##
    return clustered


def resolve_overlapping_te(entries):
    toresolve = []
    for i in entries:
        i.score = float(i.end - i.start)
        toresolve.append(i)
    toresolve = sorted(toresolve, key=lambda ts: -ts.score)
    # start with the highest score=longest
    resolved = []
    for x in toresolve:
        discard = False
        for f in resolved:
            if(x.start > f.end):
                pass
            elif(x.end < f.start):
                pass
            elif(x.start >= f.start and x.end <= f.end):
                discard = True
            elif(x.start < f.start and x.end >= f.start):
                newend = f.start - 1
                x.end = newend
                x.score = newend - x.start
            elif(x.end > f.end and x.start <= f.end):
                newstart = f.end + 1
                x.start = newstart
                x.score = x.end - newstart
            else:
                raise ValueError("This should not happen")
        if(not discard):
            resolved.append(x)
    return resolved


def format_entry(te):
    topr = []
    # chr, source, feature, start, end, score, strand, frame, comment (comment is unparsed)
    topr.append(te.chr)
    topr.append(te.source)
    topr.append(te.feature)
    topr.append(str(te.start))
    topr.append(str(te.end))
    topr.append(str(te.score))
    topr.append(te.strand)
    topr.append(te.frame)
    topr.append(te.comment)
    toret = "\t".join(topr)
    return toret

    # 2L	RepeatMasker	similarity	1724	1849	 6.4	-	.	Target "Motif:DMTRDNA" 1311 1435

##########################################


tefhashloadmethod = load_tefamhash_disablestrand

matchscore, mmpenalty = float(options.matchscore), float(options.mmpenalty)
finalentries = []
for tecluster in GTFTEClusterReader(load_gtfhash(options.gtf), 30000):

    # group TEs by family name and strand
    tefh = tefhashloadmethod(tecluster)

    temp = []
    for tefam, entries in tefh.items():
        merged = mergeTEentries(entries, matchscore)  # merge overlapping TEs of the same family
        clustered = iterativeClusterTEentries(merged, mmpenalty)  # cluster neighbouring TEs of the same family
        temp.extend(clustered)

    nuevo = resolve_overlapping_te(temp)
    finalentries.extend(nuevo)

ofh = open(options.output, "w")
for fe in finalentries:
    f = format_entry(fe)
    ofh.write(f + "\n")



