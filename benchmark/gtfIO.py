# module for repeatmasker cleaning


class GTFLightReader:
    """
    read a GTF/GFF file and obtain an unparsed list (no conversion to int or float); the comment region is also unparsed

    example of GTF
    2L	RepeatMasker	similarity	1	1433	 0.8	-	.	Target "Motif:DM_ROO" 1 1446
    2L	RepeatMasker	similarity	1724	1849	 6.4	-	.	Target "Motif:DMTRDNA" 1311 1435
    """

    def __init__(self, file):
        self.__filename = file
        self.__filehandle = open(file, "r")

    def __iter__(self):
        return self

    def next(self):
        line = ""
        while(1):
            line = self.__filehandle.readline()
            if line == "":
                raise StopIteration
            line = line.rstrip('\n')
            if line != "":
                break

        a = line.split("\t")
        return a


class GTFReader:
    """
    read a GFF/GTF file and obtain GTFEntries (class); the comment is unparsed

    example of GTF:
    2L	RepeatMasker	similarity	1	1433	 0.8	-	.	Target "Motif:DM_ROO" 1 1446
    2L	RepeatMasker	similarity	1724	1849	 6.4	-	.	Target "Motif:DMTRDNA" 1311 1435
    """

    def __init__(self, file):
        self.__filename = file
        self.__filehandle = open(file, "r")

    def __iter__(self):
        return self

    def next(self):
        line = ""
        while(1):
            line = self.__filehandle.readline()
            if line == "":
                raise StopIteration
            line = line.rstrip('\n')
            if line.startswith("#"):
                continue
            if line != "":
                break
        a = line.split("\t")
        # 2L	RepeatMasker	similarity	1724	1849	 6.4	-	.	Target "Motif:DMTRDNA" 1311 1435
        # 0          1		    2		  3       4       5     6       7          8
        e = GTFEntry(a[0], a[1], a[2], int(a[3]), int(a[4]), a[5], a[6], a[7], a[8])
        return e

    @classmethod
    def readall(cls, file):
        entries = []
        for e in GTFReader(file):
            entries.append(e)
        return entries


class GTFWriter:
    """
    Write the content of a gtf file to a file
    """

    def __init__(self, file):
        self.__filename = file
        self.__filehandle = open(file, "w")

    def write(self, entry):
        topr = [entry.chr, entry.source, entry.feature, entry.start, entry.end, entry.score, entry.strand, entry.frame, entry.comment]
        topr = [str(i) for i in topr]
        form = "\t".join(topr)
        self.__filehandle.write(form + "\n")
        return 1

    def close(self):
        self.__filehandle.close()

    @classmethod
    def write_all(cls, file, gtfentries):
        gw = GTFWriter(file)
        for e in gtfentries:
            gw.write(e)
        gw.close()
        return 1


class GTFEntry:
    """
    A minimal representation of a GTF/GFF entry; having the features
    chr, source, feature, start, end, score, strand, frame, comment (comment is unparsed)
    """

    def __init__(self, chr, source, feature, start, end, score, strand, frame, comment):
        self.chr = chr
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        if isinstance(score, str):
            score = score.lstrip(" ").rstrip(" ")
        else:
            pass
        self.score = float(score)
        self.strand = strand
        self.frame = frame
        self.comment = comment

    def __str__(self):
        topr = "{0}:{1}..{2}-{3}".format(self.chr, self.start, self.end, self.feature)
        return topr

