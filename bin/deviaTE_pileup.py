#! /usr/local/bin/python3

import warnings
import pandas
from collections import Counter

# set up constants
MATCH = 0
INSERTION = 1
DELETION = 2
REF_SKIP = 3
SOFT_CLIP = 4
EQUAL = 7

uniq_nuc = {'A', 'C', 'G', 'T'}
ambig_nuc = {'V', 'Y', 'H', 'D', 'N'}

class Site:

    def __init__(self, pos, refbase, sid, fam):

        if refbase not in uniq_nuc:
            warnings.warn('Reference sequence contains ambiguous nucleotide: ' + refbase)

        self.TEfam = fam
        self.sample_id = sid
        self.pos = pos
        self.refbase = refbase
        self.A = 0
        self.C = 0
        self.G = 0
        self.T = 0
        self.cov = 0
        self.snp = False
        self.refsnp = False
        self.int_del = []
        self.ins = []
        self.delet = []
        self.trunc_left = 0
        self.trunc_right = 0
        self.annotation = []

    def sum_coverage(self):
        self.cov = self.A + self.C + self.G + self.T

    def is_snp(self, min_count, min_freq, A, C, G, T, cov):

        nuc = {'A' : A, 'C' : C, 'G' : G, 'T' : T}

        # check if site is a polymorphic SNP
        # only if coverage is higher than most abundant base
        if cov > max(nuc.values()):
            # check counts
            alt_counts = {base : count for base, count in nuc.items() if base is not self.refbase}
            alt_freqs = {base : (count / cov) for base, count in alt_counts.items()}

            if any(x >= min_count for x in alt_counts.values()):
                if any(x >= min_freq for x in alt_freqs.values()):
                    self.snp = True

        # check if site is reference snp
        # needs 0 at ref, and coverage equal to most abundant base, but ne 0
        if self.refbase in uniq_nuc:
            #ref_count = getattr(self, self.refbase)
            ref_count = nuc[self.refbase]
            if ref_count == 0:
                if cov is not 0:
                    if cov == max(nuc.values()):
                        self.refsnp = True
                        
        # return val for unit test
        return (self.snp, self.refsnp)

    def filter_IND(self, att, min_count, norm_fac):
        # grab attribute to filter
        attr = getattr(self, att)

        if len(attr) is not 0:
            cnt = Counter(attr)
            feat = list(cnt.items())
            keep = []

            for i in feat:
                if i[1] >= min_count:
                    keep.append(i)

            if len(keep) is 0:
                setattr(self, att, 'NA')
            else:
                upt = Site.reformat_tuple(keep, norm_factor=norm_fac)
                setattr(self, att, upt)
        else:
            setattr(self, att, 'NA')

    def filter_trunc(self, min_trunc_count):
        if self.trunc_left < min_trunc_count:
            self.trunc_left = 'NA'
        if self.trunc_right < min_trunc_count:
            self.trunc_right = 'NA'

    def check_annotation(self, anno):
        # if no annotation was provided
        # or no elements found for this family
        if len(anno) is 0:
            self.annotation = 'NA'
        else:
            for i in anno:
                if self.pos >= int(i[1]) - 1 and self.pos <= int(i[2]) - 1:
                    self.annotation = i[0]

        # if there were annotations, but site is not in one
        if len(self.annotation) is 0:
            self.annotation = 'intergenic'

    def normalize(self, norm_factor):
        # normalizes all counts for sequencing depth
        self.A = self.A / norm_factor
        self.C = self.C / norm_factor
        self.G = self.G / norm_factor
        self.T = self.T / norm_factor
        self.cov = self.cov / norm_factor

        if self.trunc_left is not 'NA':
            self.trunc_left = self.trunc_left / norm_factor
        if self.trunc_right is not 'NA':
            self.trunc_right = self.trunc_right / norm_factor

    @staticmethod
    def write_frame(sites, out):
        # create a list of all object instances
        # and turn into a pandas frame
        site_list = [x.__dict__ for x in sites]
        fr = pandas.DataFrame(site_list)

        # order the columns, introduce hash and print
        fr = fr[['TEfam', 'sample_id', 'pos', 'refbase', 'A', 'C', 'G', 'T', 'cov', 'snp', 'refsnp',
                 'int_del', 'trunc_left', 'trunc_right', 'ins', 'delet', 'annotation']]

        fr = fr.rename(columns={'TEfam': '#TEfam'})

        fr.to_csv(out, index=False, sep=' ', mode='w')

    @staticmethod
    def reformat_tuple(tup, norm_factor):
        # takes a list of tuples
        # returns reformatted string for printing in tsv
        r = ''
        for i in tup:
            r = r + str(i[0][0]) + ':' + str(i[0][1]) + ':' + str(i[1] / norm_factor) + ','

        return r[:-1]


def perform_pileup(sitelist, bam, min_int_del_len, min_trunc_len, min_indel_len):
    # initiate lookup sets
    readdump_int_del = set()
    readdump_trunc = set()
    readdump_indels = set()

    # for all covered positions in the reference
    for pileupcolumn in bam.pileup(Site.name, truncate=True): 
        # for each read at this pos
        for pileupread in pileupcolumn.pileups:

            # coverage and bases
            # only get nucleotides if no deletion or splice
            if (pileupread.is_del == 0) and (pileupread.is_refskip == 0):

                # count the base at the queryposition of this read
                nt = pileupread.alignment.query_sequence[pileupread.query_position]
                ntu = nt.upper()

                if ntu in ambig_nuc:
                    # warnings.warn('ignoring ambiguous base in read: ' + ntu)
                    pass
                elif ntu in uniq_nuc:
                    site = sitelist[pileupcolumn.pos]
                    setattr(site, ntu, getattr(site, ntu) + 1)

                else:
                    warnings.warn('ignoring unknown base in read: ' + ntu)

            # internal deletions
            # extract cigarstring of read to check if spliced
            read_cigar = pileupread.alignment.cigarstring
            read_name = pileupread.alignment.query_name

            if ('N' in read_cigar or 'D' in read_cigar) and (read_name not in readdump_int_del):
                readdump_int_del.add(read_name)
                read_tuple = pileupread.alignment.cigartuples
                read_start = pileupread.alignment.reference_start
                new_int_del = eval_cigartuple_int_del(read_tuple, read_start, min_int_del_len)
                # for every int_del detected in this read
                for i in new_int_del:
                    site = sitelist[i[0]]
                    curr = site.int_del
                    curr.append(i)
                    site.int_del = curr

            # truncations/soft clips
            # similar to above
            if 'S' in read_cigar and (read_name not in readdump_trunc):
                readdump_trunc.add(read_name)
                read_tuple = pileupread.alignment.cigartuples
                read_start = pileupread.alignment.reference_start
                read_end = pileupread.alignment.reference_end - 1
                new_truncs = eval_cigartuple_trunc(read_tuple, read_start, read_end, min_trunc_len)
                # for every trunc in this read increment directional counter
                for i in new_truncs:
                    if i[1] is 'l':
                        site = sitelist[i[0]]
                        site.trunc_left = site.trunc_left + 1

                    if i[1] is 'r':
                        site = sitelist[i[0]]
                        site.trunc_right = site.trunc_right + 1

            # indels
            if ('I' in read_cigar or 'D' in read_cigar) and (read_name not in readdump_indels):
                readdump_indels.add(read_name)
                read_tuple = pileupread.alignment.cigartuples
                read_start = pileupread.alignment.reference_start - 1
                new_indels = eval_cigartuple_indel(read_tuple, read_start, min_indel_len, min_int_del_len)
                for i in new_indels:
                    if i[-1] is 'D':
                        site = sitelist[i[0] - 1]
                        curr = site.delet
                        curr.append((i[0], i[1]))
                        site.delet = curr

                    if i[-1] is 'I':
                        site = sitelist[i[0] - 1]
                        curr = site.ins
                        curr.append((i[0], i[1]))
                        site.ins = curr


def eval_cigartuple_int_del(cigartuple, read_start, min_int_del_len):
    int_del_shift = 0
    int_dels = []

    for cig_id, length in cigartuple:
        if cig_id == DELETION or cig_id == REF_SKIP:
            int_del_start = read_start + int_del_shift
            int_del_end = int_del_start + int(length)

            if length >= min_int_del_len:
                int_dels.append((int_del_start, int_del_end))

            int_del_shift += length

        elif cig_id == SOFT_CLIP or cig_id == INSERTION:
            continue
        elif cig_id == MATCH or cig_id == EQUAL:
            int_del_shift += length
        else:
            raise ValueError('Cigarstring contains unusual symbol: ' + cig_id)

    return int_dels


def eval_cigartuple_trunc(cigartuple, read_start, read_end, min_trunc_len):
    truncs = []

    left_id = cigartuple[0][0]
    left_len = cigartuple[0][1]
    right_id = cigartuple[-1][0]
    right_len = cigartuple[-1][1]

    if left_id is SOFT_CLIP and left_len >= min_trunc_len:
        truncs.append((read_start, 'l'))

    if right_id is SOFT_CLIP and right_len >= min_trunc_len:
        truncs.append((read_end, 'r'))

    return truncs


def eval_cigartuple_indel(cigartuple, read_start, min_indel_len, min_int_del_len):
    indel_shift = 0
    indels_read = []

    for cig_id, length in cigartuple:
        if cig_id == DELETION:
            del_start = read_start + indel_shift
            del_end = del_start + int(length)

            if (length >= min_indel_len) and (length < min_int_del_len):
                indels_read.append((del_start, del_end, 'D'))
            
            indel_shift += length

        elif cig_id == INSERTION:
            ins_start = read_start + indel_shift
            ins_end = ins_start + int(length)

            if (length >= min_indel_len):
                indels_read.append((ins_start, ins_end, 'I'))


        elif cig_id == SOFT_CLIP or cig_id == REF_SKIP:
            continue
        elif cig_id == MATCH or cig_id == EQUAL:
            indel_shift += length
        else:
            print(cig_id)
            raise ValueError('Cigarstring contains unusual symbol: ' + cig_id)

    return indels_read
