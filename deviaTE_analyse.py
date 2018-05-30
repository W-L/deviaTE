#! /usr/local/bin/python3

import argparse
import sys
import pysam
import warnings
from bin import deviaTE_pileup as pileup

deviaTE = '/'.join(sys.argv[0].split('/')[:-1])


# set up parser and arguments
parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, required=True, dest='bamfile', default=None, help='alignment file to be analyzed')
parser.add_argument('--family', type=str, required=True, dest='family_sel', default=None, help='TE family to visualize; a header in the reference library')
parser.add_argument('--library', type=str, dest='te_library', default=deviaTE + '/library/te_library', help='path to reference sequence library')
parser.add_argument('--output', type=str, default=None, help='name of output table')
parser.add_argument('--sample_id', type=str, default=None, help='sample identifier')
parser.add_argument('--annotation', type=str, default=None, help='annotation in gff-format')
parser.add_argument('--log', type=str, default=None, help='logfile from preparation script')
args = parser.parse_args()

if args.sample_id is None:
    args.sample_id = args.bamfile
    
# set default output if none given
if args.output is None:
    args.output = args.sample_id + '.' + args.family_sel

# look for the refseq of chosen family
refs = open(args.te_library, 'r')
refseq = None

for line in refs:
    if line.startswith('>'):
        family = line.replace('>', '').rstrip('\n')
        if family == args.family_sel:
            refseq = refs.readline().rstrip('\n')
            refs.close()
            break

if refseq is None:
    raise ValueError('Selected family is not in references. Typo? ' + args.family_sel)
else:
    # set family as class variable
    pileup.Site.name = args.family_sel


# get annotation for the selected family
fam_anno = []
if args.annotation is not None:
    anno = open(args.annotation, 'r')

    for line in anno:
        if line.startswith(args.family_sel):
            entry = line.split('\t')
            fam_anno.append(tuple(entry[2:5]))

    if len(fam_anno) is 0:
        warnings.warn('no annotaions found for: ' + args.family_sel)

    anno.close()

# search the logfile, if provided, for normalization factor
norm_fac = 1
if args.log is not None:
    log = open(args.log, 'r')

    for line in log:
        if '#total_read_length:' in line:
            norm_fac = int(line.rstrip('\n').split(' ')[-1])
            norm_fac = norm_fac / (10 ** 6)
            break

    log.close()


# create a list of Site instances
sites = []
c = 0
for base in refseq:
    sites.append(pileup.Site(pos=c, refbase=base, sid=args.sample_id, fam=args.family_sel))
    c += 1


# open the alignment file for data extraction
bamfile_op = pysam.AlignmentFile(args.bamfile, 'rb')

pileup.perform_pileup(sitelist=sites, bam=bamfile_op, min_int_del_len=20,
                      min_trunc_len=10, min_indel_len=2)

bamfile_op.close()

# process the extracted data
for site in sites:
    site.sum_coverage()
    site.is_snp(min_count=5, min_freq=0.05, A=site.A, C=site.C, G=site.G, T=site.T, cov=site.cov)
    site.filter_IND(att='int_del', min_count=3, norm_fac=norm_fac)
    site.filter_IND(att='ins', min_count=5, norm_fac=norm_fac)
    site.filter_IND(att='delet', min_count=5, norm_fac=norm_fac)
    site.filter_trunc(min_trunc_count=3)
    site.check_annotation(anno=fam_anno)
    if norm_fac is not 1:
        site.normalize(norm_factor=norm_fac)

# write the results
pileup.Site.write_frame(sites=sites, out=args.output)
print('analysis completed; table written to: ' + args.output)
