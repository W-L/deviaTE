#!/usr/bin/env python3

import sys
import pysam

limit = sys.argv[1]
out = sys.argv[2]

# remove reads under alignment length limit
inpfile = pysam.AlignmentFile('tmp', 'r')
outfile = pysam.AlignmentFile(out, 'w', template=inpfile)

for read in inpfile.fetch():
    al_len = read.query_alignment_length
    if al_len >= int(limit):
        outfile.write(read)

inpfile.close()
outfile.close()
