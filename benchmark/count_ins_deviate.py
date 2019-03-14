#!/usr/bin/env python3

import argparse
from glob import glob

parser = argparse.ArgumentParser()
parser.add_argument('--inp')
args = parser.parse_args()


class TE:

    def __init__(self, deviaTE_file):
        self.file = deviaTE_file
        self.cov = float()
        self.hqcov = float()

    def readTE(self):
        with open(self.file, 'r') as infile:
            while True:
                nline = infile.readline()
                cl = nline.split(' ')
                if cl[0][0] == ('#'):
                    if cl[1] == 'insertions/haploid:':
                        self.insertions = float(cl[2])
                        self.hq_insertions = float(cl[4])
                        infile.readline()  # skip last headerline
                        break

            for line in infile:
                l = line.split(' ')
                self.name = l[0]
                self.cov += float(l[8])
                self.hqcov += float(l[10])


def write_table(TElist, outfile):
    with open(outfile, 'w') as out:

        for i in TElist:
            line = (i.name
                    + ' ' + str(round(i.insertions, 3))
                    + ' ' + str(round(i.hq_insertions, 3))
                    + ' ' + str(round(i.cov, 3))
                    + ' ' + str(round(i.hqcov, 3))
                    + '\n')
            out.write(line)


TEs = []

files = glob(args.inp + '.*')

for f in files:
    if f[-4:] != '.out':
        te = TE(f)
        te.readTE()
        TEs.append(te)

write_table(TElist=TEs, outfile=args.inp + '.out')

