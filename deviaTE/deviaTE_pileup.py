#!/usr/bin/env python3

import warnings
import pandas
import pysam
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


class Sample:

    def __init__(self, name, fam, lib, anno, bam):
        self.name = name
        self.sites = list()
        self.int_dels = list()
        self.fam = fam
        self.lib = lib
        self.anno = anno
        self.bam = bam

    def get_ref(self):
        refs = open(self.lib, 'r')
        self.refseq = None

        for line in refs:
            if line.startswith('>'):
                family = line.replace('>', '').rstrip('\n')
                if family == self.fam:
                    refseq = refs.readline().rstrip('\n')
                    refs.close()
                    self.refseq = refseq
                    break

    def get_anno(self):
        self.fam_anno = []
        anno_file = open(self.anno, 'r')

        for line in anno_file:
            if line.startswith(self.fam):
                entry = line.split('\t')
                self.fam_anno.append(tuple(entry[2:5]))
        anno_file.close()
            
    def get_norm_fac_rpm(self):
        # count number of reads
        bamfile_op = pysam.AlignmentFile(self.bam, 'rb')
        c = bamfile_op.mapped + bamfile_op.unmapped
        bamfile_op.close()
        return(c / (10**6))
        

    def perform_pileup(self, hq_threshold):
        # initiate lookup sets
        readdump_int_del = set()
        readdump_trunc = set()
        readdump_indels = set()

        bamfile_op = pysam.AlignmentFile(self.bam, 'rb')

        # for all covered positions in the reference
        for pileupcolumn in bamfile_op.pileup(contig=self.fam, truncate=True, max_depth=1000000):
            # for each read at this pos
            for pileupread in pileupcolumn.pileups:
                
                pr = Pileupread(isdel=pileupread.is_del, isref=pileupread.is_refskip,
                                qseq=pileupread.alignment.query_sequence, qpos=pileupread.query_position,
                                colpos=pileupcolumn.pos, cig_string=pileupread.alignment.cigarstring,
                                qname=pileupread.alignment.query_name, cig_tuples=pileupread.alignment.cigartuples,
                                ref_start=pileupread.alignment.reference_start, ref_end=pileupread.alignment.reference_end,
                                mapq=pileupread.alignment.mapping_quality)

                # base
                if (pr.is_del == 0) and (pr.is_refskip == 0):
                    pr.count_nucleotide(sample_sites=self.sites)
                    pr.count_hq_coverage(sample_sites=self.sites, hqt=hq_threshold)

                # internal deletions
                if ('N' in pr.cigar_string or 'D' in pr.cigar_string):
                    if pr.query_name not in readdump_int_del:
                        readdump_int_del.add(pr.query_name)
                        pr.eval_int_del(sample_sites=self.sites)

                # truncations
                if 'S' in pr.cigar_string:
                    if pr.query_name not in readdump_trunc:
                        readdump_trunc.add(pr.query_name)
                        pr.eval_trunc(sample_sites=self.sites)

                # indels
                if ('I' in pr.cigar_string or 'D' in pr.cigar_string):
                    if pr.query_name not in readdump_indels:
                        readdump_indels.add(pr.query_name)
                        pr.eval_indel(sample_sites=self.sites)

        bamfile_op.close()

    def mean_read_length(self):
        bamfile_op = pysam.AlignmentFile(self.bam, 'rb')
        read_lengths = 0
        c = 0
        for read in bamfile_op:
            if read.is_unmapped is False:
                read_lengths += read.query_length
                c += 1
        bamfile_op.close()
        return(read_lengths / c)

    def collect_int_dels(self):
        for s in self.sites:
            if s.int_del is not 'NA':
                if ',' in s.int_del:
                    multi_csd = s.int_del.split(',')
                    for csd in multi_csd:
                        j = csd.split(':')
                        self.int_dels.append(Int_del(start=int(j[0]), end=int(j[1]), abundance=float(j[2])))
                else:
                    j = s.int_del.split(':')
                    self.int_dels.append(Int_del(start=int(j[0]), end=int(j[1]), abundance=float(j[2])))

    def calc_phys_cov(self):
        # add abundance of deletion to any spanned site
        for site in self.sites:
            for int_del in self.int_dels:
                if site.pos in int_del.range:
                    site.phys_cov += int_del.abundance
                    
    def get_norm_fac_scg(self, genes):
        bamfile_op = pysam.AlignmentFile(self.bam, 'rb')
        # dict of refs and their len
        ref_dict = dict(zip(bamfile_op.references,bamfile_op.lengths))
        
        genes = genes.split(',')
        av_cov_genes = list()
        for g in genes:
            # calc average cov of single copy gene
            sum_cov = sum([len(pileupcolumn.pileups) for pileupcolumn in bamfile_op.pileup(contig=g, truncate=True)])
            if sum_cov > 0:
                av_cov_genes.append(sum_cov / ref_dict[g])
        
        # average across multiple genes
        norm_fac = sum(av_cov_genes) / len(av_cov_genes)
        bamfile_op.close()
        return(norm_fac)
        
    def estimate_insertions(self, norm_factor):
        # average cov of TE
        av_cov = average_cov(sitelist=self.sites, start=1, end=len(self.sites))
        # normalize with single copy gene to obtain copy number per haploid
        ihat = av_cov / norm_factor
        return(ihat)
        
    def write_frame(self, out, insertions, command, t, norm):
        # create a list of all object instances
        # and turn into a pandas frame
        site_list = [x.__dict__ for x in self.sites]
        fr = pandas.DataFrame(site_list)

        # order the columns, add hash and print
        fr = fr[['TEfam', 'sample_id', 'pos', 'refbase', 'A', 'C', 'G', 'T', 'cov',
                 'phys_cov', 'hq_cov', 'snp', 'refsnp', 'int_del', 'int_del_freq',
                 'trunc_left', 'trunc_right', 'ins', 'delet', 'annotation']]
        fr = fr.rename(columns={'TEfam': '#TEfam'})
        # add a line with the insertions
        with open(out, 'w') as outfile: 
            outfile.write('# ' + t + ', command: ' + command + ', norm: ' + norm + '\n')
            outfile.write('# insertions/haploid: ' + str(insertions) + '\n') 
        fr.to_csv(out, index=False, sep=' ', mode='a')


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
        self.phys_cov = 0
        self.hq_cov = 0
        self.snp = False
        self.refsnp = False
        self.int_del = []
        self.int_del_freq = 'NA'
        self.ins = []
        self.delet = []
        self.trunc_left = 0
        self.trunc_right = 0
        self.annotation = []

    def sum_coverage(self):
        self.cov = self.A + self.C + self.G + self.T

    def is_snp(self, min_count, min_freq, A, C, G, T, cov):
        # sets self.snp of self.refsnp to True
        # awkward args for unit test
        nuc = {'A': A, 'C': C, 'G': G, 'T': T}

        # polymorphic SNP = coverage is higher than most abundant base
        if cov > max(nuc.values()):
            alt_counts = {base: count for base, count in nuc.items() if base is not self.refbase}
            alt_freqs = {base: (count / cov) for base, count in alt_counts.items()}

            if any(x >= min_count for x in alt_counts.values()):
                if any(x >= min_freq for x in alt_freqs.values()):
                    self.snp = True

        # reference snp = 0 at ref nuc, and coverage equal to most abundant base, but not 0
        if self.refbase in uniq_nuc:
            if nuc[self.refbase] == 0:
                if cov is not 0 and cov is not 0.0:
                    if cov == max(nuc.values()):
                        self.refsnp = True

        # return val for unit test
        return (self.snp, self.refsnp)

    def filter_IND(self, att, min_count):
        # grab attribute to filter
        attr = getattr(self, att)

        if len(attr) is not 0:
            cnt = Counter(attr)
            feat = list(cnt.items())  # list of tuples, ((start, end), count)
            keep = []

            for i in feat:
                if i[1] >= min_count:
                    keep.append(i)

            if len(keep) is 0:
                setattr(self, att, 'NA')
            else:
                upt = reformat_tuple(keep)
                setattr(self, att, upt)
        else:
            setattr(self, att, 'NA')

    def filter_trunc(self, min_trunc_count):
        if self.trunc_left < min_trunc_count:
            self.trunc_left = 'NA'
        if self.trunc_right < min_trunc_count:
            self.trunc_right = 'NA'

    def check_annotation(self, anno):
        if len(anno) is 0:
            self.annotation = 'NA'
        else:
            for i in anno:
                if self.pos >= int(i[1]) - 1 and self.pos <= int(i[2]) - 1:
                    self.annotation = i[0]

        # if site is not in annotation
        if len(self.annotation) is 0:
            self.annotation = 'intergenic'

    def normalize(self, norm_factor):
        # normalizes all counts for sequencing depth
        # int_del abundance, ins, del
        self.A = round(self.A / norm_factor, 3)
        self.C = round(self.C / norm_factor, 3)
        self.G = round(self.G / norm_factor, 3)
        self.T = round(self.T / norm_factor, 3)
        self.cov = round(self.cov / norm_factor, 3)
        self.phys_cov = round(self.phys_cov / norm_factor, 3)
        self.hq_cov = round(self.hq_cov / norm_factor, 3)
        if self.trunc_left is not 'NA':
            self.trunc_left = round(self.trunc_left / norm_factor, 3)
        if self.trunc_right is not 'NA':
            self.trunc_right = round(self.trunc_right / norm_factor, 3)
            
        self.norm_feature(attr='int_del', norm_factor=norm_factor)
        self.norm_feature(attr='ins', norm_factor=norm_factor)
        self.norm_feature(attr='delet', norm_factor=norm_factor)
        
    def norm_feature(self, attr, norm_factor):
        feat = getattr(self, attr)
        if feat is not 'NA':
            feat_new = ''
            if ',' in feat:
                multi_csd = feat.split(',')
                for csd in multi_csd:
                    j = csd.split(':')
                    feat_new = str(feat_new) + str(j[0]) + ':' + str(j[1]) + ':' + str(round((float(j[2]) / float(norm_factor)), 3)) + ','
                setattr(self, attr, feat_new[:-1])
                
            else:
                j = feat.split(':')
                setattr(self, attr, str(j[0]) + ':' + str(j[1]) + ':' + str(round((float(j[2]) / float(norm_factor)), 3)))


class Int_del:

    def __init__(self, start, end, abundance):
        self.start = start
        self.end = end
        self.pos = (start, end)
        self.abundance = abundance
        self.range = range(start + 1, end)

    def est_freq(self, sites, corr_factor):
        # estimate frequency as abundance div by average coverage spanned
        mean_cov = average_cov(sitelist=sites, start=self.start + 1, end=self.end - 1)
        self.freq = (self.abundance / mean_cov) ** corr_factor

    def write_freq(self, sites):
        # parse new column at start pos of int del
        id_site = sites[self.start]
        if id_site.int_del_freq == 'NA':
            id_site.int_del_freq = str(self.start) + ':' + str(self.end) + ':' + str(self.freq)
        else:
            id_site.int_del_freq = id_site.int_del_freq + ',' + str(self.start) + ':' + str(self.end) + ':' + str(self.freq)


class Pileupread:

    def __init__(self, isdel, isref, qseq, qpos, colpos, cig_string, qname, cig_tuples, ref_start, ref_end, mapq):
        self.is_del = isdel
        self.is_refskip = isref
        self.query_seq = qseq
        self.query_pos = qpos
        self.column_pos = colpos
        self.cigar_string = cig_string
        self.query_name = qname
        self.cigar_tuple = cig_tuples
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.mapq = mapq

    def count_nucleotide(self, sample_sites):
        # count the base at the queryposition of this read
        # increase counter of nucleotide at column pos
        nt = self.query_seq[self.query_pos].upper()

        if nt in ambig_nuc:
            # warnings.warn('ignoring ambiguous base in read: ' + nt)
            pass
        elif nt in uniq_nuc:
            site = sample_sites[self.column_pos]
            setattr(site, nt, getattr(site, nt) + 1)

        else:
            warnings.warn('ignoring unknown base in read: ' + nt)
            
    def count_hq_coverage(self, sample_sites, hqt):
        # count coverage only above threshold
        site = sample_sites[self.column_pos]
        mapping_qual = site.hq_cov
        
        if self.mapq >= hqt:
            site.hq_cov = site.hq_cov + 1
        else:
            pass

    def eval_int_del(self, sample_sites):
        int_del_shift = 0

        for cig_id, length in self.cigar_tuple:
            if cig_id == DELETION or cig_id == REF_SKIP:
                int_del_start = self.ref_start + int_del_shift
                int_del_end = int_del_start + int(length)

                if length >= 20:
                    site = sample_sites[int_del_start]
                    curr = site.int_del
                    curr.append((int_del_start, int_del_end))  # append tuples of (start, end)
                    site.int_del = curr

                int_del_shift += length

            elif cig_id == SOFT_CLIP or cig_id == INSERTION:
                continue
            elif cig_id == MATCH or cig_id == EQUAL:
                int_del_shift += length
            else:
                raise ValueError('Cigarstring contains unusual symbol: ' + cig_id)

    def eval_trunc(self, sample_sites):
        truncs = []
        left_id = self.cigar_tuple[0][0]
        left_len = self.cigar_tuple[0][1]
        right_id = self.cigar_tuple[-1][0]
        right_len = self.cigar_tuple[-1][1]

        if left_id is SOFT_CLIP and left_len >= 10:
            site = sample_sites[self.ref_start]
            site.trunc_left = site.trunc_left + 1

        if right_id is SOFT_CLIP and right_len >= 10:
            site = sample_sites[self.ref_end - 1]
            site.trunc_right = site.trunc_right + 1

    def eval_indel(self, sample_sites):
        indel_shift = 0
        indels_read = []

        for cig_id, length in self.cigar_tuple:
            if cig_id == DELETION:
                del_start = self.ref_start - 1 + indel_shift
                del_end = del_start + int(length)

                if (length >= 2) and (length < 20):
                    site = sample_sites[del_start - 1]
                    curr = site.delet
                    curr.append((del_start, del_end))
                    site.delet = curr

                indel_shift += length

            elif cig_id == INSERTION:
                ins_start = self.ref_start - 1 + indel_shift
                ins_end = ins_start + int(length)

                if (length >= 2):
                    site = sample_sites[ins_start - 1]
                    curr = site.ins
                    curr.append((ins_start, ins_end))
                    site.ins = curr

            elif cig_id == SOFT_CLIP or cig_id == REF_SKIP:
                continue
            elif cig_id == MATCH or cig_id == EQUAL:
                indel_shift += length
            else:
                raise ValueError('Cigarstring contains unusual symbol: ' + cig_id)


def reformat_tuple(tup):
    # takes a list of tuples
    # returns reformatted string for printing in tsv
    r = ''
    for i in tup:
        r = r + str(i[0][0]) + ':' + str(i[0][1]) + ':' + str(i[1]) + ','
    return r[:-1]


def average_cov(sitelist, start, end):
    # returns mean cov in specified region
    sum_cov = 0
    region = sitelist[start:end + 1]
    for s in region:
        sum_cov = sum_cov + s.cov + s.phys_cov
    return(sum_cov / len(region))


def correction_factor(x):
    # coefficients from model in R
    coefs = [2.198139e-01, 4.182264e-03, -7.480947e-06, 4.379848e-09]
    y = coefs[3] * (x ** 3) + coefs[2] * (x ** 2) + coefs[1] * x + coefs[0]
    return(y)








