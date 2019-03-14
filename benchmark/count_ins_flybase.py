#! /usr/bin/env python3

import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('--inp')
parser.add_argument('--table')
args = parser.parse_args()


def write_dict(dic, outname):
    with open(outname, 'w') as outfile:
        for i in dic.items():
            line = i[0] + ' ' + str(i[1]) + '\n'
            outfile.write(line)


def get_names(table):
    TEdict = defaultdict(str)

    with open(table, 'r') as tetable:
        tetable.readline()  # skip header
        for line in tetable:
            cl = line.split('\t')
            new_name = cl[0]
            old_name = cl[2]
            TEdict[old_name] = new_name

    return(TEdict)


# get dictionary for names first
TEdict = get_names(table=args.table)

# collect all insertions
TEs = defaultdict(int)

with open(args.inp, 'r') as infile:
    for line in infile:
        l = line.split(' ')
        n = l[3]
        if 'name=' in n:
            if '{}' in n:
                n_sliced = n.split('=')[1].split('{}')[0]
                # search for n_sliced in dictionary
                if n_sliced in TEdict.keys():
                    newname = TEdict[n_sliced]
                    TEs[newname] += 1
                else:
                    # first try: element was not in the verbatim - try with appending -element
                    ntry = n_sliced + '-element'
                    if ntry in TEdict.keys():
                        newname = TEdict[ntry]
                        TEs[newname] += 1
                        continue

                    # all the exceptions:
                    if ntry == 'H-element':
                        newname = 'DMHFL1'
                        TEs[newname] += 1
                        continue

                    if ntry == 'ninja-Dsim-like-element':
                        newname = 'DSRN'
                        TEs[newname] += 1
                        continue

                    if ntry == 'DM88-element':
                        newname = 'DM88'
                        TEs[newname] += 1
                        continue

                    if ntry in ['HeT-Tag-element', 'Xanthias-element', 'Y-element']:
                        newname = 'unknown'
                        TEs[newname] += 1
                        continue

                    print("Old name not in dict")
            else:
                print('different format of annotation line')
        else:
            print('wrong field of the annotation header')

write_dict(dic=TEs, outname=args.inp + '.out')





