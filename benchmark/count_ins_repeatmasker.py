#!/usr/bin/env python3

import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('--inp')
parser.add_argument('--lib')
parser.add_argument('--assembly')
args = parser.parse_args()


class GTFline:

    def __init__(self, l):
        # 2R	RepeatMasker	similarity	2156697	2158764	2067	+	.	Target "Motif:DMIFACA" 1673 3796
        lsplit = l.split('\t')
        self.sites = int(lsplit[4]) - int(lsplit[3])
        name = lsplit[8].split(' ')[1]
        name = name.replace('"Motif:', '').rstrip('"')
        self.name = name


class TE:

    def __init__(self, name):
        self.name = name
        self.insertions = 0
        self.sites = 0

    def proportion(self, genomesize):
        self.proportion = self.sites / genomesize


def get_genome_size(genome):
    size = 0
    with open(genome, 'r') as gfile:
        for line in gfile:
            if line.startswith('>') is False:
                line = line.rstrip('\n')
                nuc = len(line)
                size += nuc
    return(size)


def write_dict(dic, outname):
    with open(outname, 'w') as outfile:
        for i in dic.items():
            te = i[1]
            line = te.name + ' ' + str(te.insertions) + ' ' + str(te.sites) + ' ' + str(te.proportion) + '\n'
            outfile.write(line)

#################################################################


TEs = defaultdict(TE)

with open(args.lib, 'r') as library:
    for line in library:
        if line.startswith('>'):
            cl = line.replace('>', '').rstrip('\n').rstrip(' ')
            TEs[cl] = TE(name=cl)


with open(args.inp, 'r') as infile:
    for line in infile:
        gffline = GTFline(l=line)

        # add sites and increase insertion counter
        currTE = TEs[gffline.name]
        currTE.insertions += 1
        currTE.sites += gffline.sites

gs = get_genome_size(genome=args.assembly)

for TE in TEs.items():
    TE[1].proportion(genomesize=gs)

write_dict(dic=TEs, outname=args.inp + '.out')
