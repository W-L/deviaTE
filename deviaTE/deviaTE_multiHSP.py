#!/usr/bin/env python3

import warnings
import itertools

# set up constants
MATCH = 0
INSERTION = 1
DELETION = 2
REF_SKIP = 3
SOFT_CLIP = 4

cig_dict = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'S'}

class Multihit:

    def __init__(self, read_id, hsp_list, fam):
        self.id = read_id
        self.hsps = list(hsp_list)
        self.fam = fam

        # sort the segments based on their first aligned position
        self.hsps.sort(key=lambda s: s.al_start)

    def create_MACs(self):
        self.MACs = []

        for L in range(1, len(self.hsps) + 1):
            for subset in itertools.combinations(self.hsps, L):
                self.MACs.append(MAC(read_id=self.id, hsp_list=subset, fam=self.fam))

    def find_hMAC(self):
        # find highest-scoring MAC for a Multihit
        hs = 0
        self.hMAC = []

        for mac in self.MACs:
            if mac.valid is True:
                if mac.score > hs:
                    hs = mac.score
                    self.hMAC = mac

        self.hMAC_score = hs


class MAC(Multihit):

    def __init__(self, read_id, hsp_list, fam):

        self.id = read_id
        self.hsps = list(hsp_list)
        self.n_hsp = len(hsp_list)
        self.fam = fam
        self.valid = True
        self.read_ranges = []
        self.ref_ranges = []
        self.score = 0

        # sort the segments based on their first aligned position
        self.hsps.sort(key=lambda s: s.al_start)

    def construct(self):
        for hsp in self.hsps:
            self.read_ranges.append(hsp.read_range)
            self.ref_ranges.append(hsp.ref_range)

    def check_overlap(self, limit):
        # checks for overlaps in the read and ref ranges
        if self.n_hsp is 1:
            return()
        
        elif self.n_hsp > 1:
            for subset in itertools.combinations(self.ref_ranges, 2):
                ref_overlap = set.intersection(*subset)
                
                # check for chiasma reads
                if min(subset[0]) >= (max(subset[1]) - 5):
                    self.valid = False

                if len(ref_overlap) > limit:
                    self.valid = False

            for subset in itertools.combinations(self.read_ranges, 2):
                read_overlap = set.intersection(*subset)

                if len(read_overlap) > limit:
                    self.valid = False

        else:
            warnings.warn('mac contains no hsp')

    def check_distance(self, limit):
        # checks for too long distances between hsps in the MAC
        # if one of them is too long -> valid = False
        if self.n_hsp is 1:
            return()

        elif self.n_hsp > 1:
            # flatten the read_ranges to find min/max
            position_list = [pos for hsp in self.read_ranges for pos in hsp]
            full_span = list(range(min(position_list), max(position_list) + 1))

            # find non-overlap of segments with NAND logic
            non_overlap = nand(full_span, position_list)
            non_overlap_ranges = get_ranges(non_overlap)

            for cont in non_overlap_ranges:
                if len(cont) > limit:
                    self.valid = False

        else:
            warnings.warn('mac contains no hsp')
            
    def score_MAC(self):
        # evaluates a score for every mac - depending on non-overlapping al_len
        if self.n_hsp is 1:
            self.score = self.hsps[0].al_len
            
        elif self.n_hsp > 1:
            score_set = set()
            for range_set in self.read_ranges:
                score_set = score_set.union(range_set)
            
            self.score = len(score_set)
        
        else:
            warnings.warn('mac contains no hsp')
            
        
    def build_cigar(self):
        # used on hMAC to generate cigarstring
        cigarbox = ''

        first_clip = ''
        if self.hsps[0].cigartuples[0][0] == SOFT_CLIP:
            first_clip = 'S' * self.hsps[0].cigartuples[0][1]

        last_clip = ''
        if self.hsps[-1].cigartuples[-1][0] == SOFT_CLIP:
            last_clip = 'S' * self.hsps[-1].cigartuples[-1][1]

        cigarbox = cigarbox + first_clip

        # initiate first position in read and ref
        cig_read_pos = self.hsps[0].al_start
        cig_ref_pos = self.hsps[0].ref_start

        for hsp in self.hsps:
            # print(hsp.cigartuples)
            # for each segment get a cigarstack first
            cig_stack = ''

            for cig_id, length in hsp.cigartuples:
                cig_stack = cig_stack + cig_dict[cig_id] * length

            cig_stack = list(enumerate(list(cig_stack)))
            cig_stack = [(pos, s) for (pos, s) in cig_stack if s is not 'S']

            # then check whether to add to cigarbox straight away
            # or pop overlapping elements first
            hsp_read_pos = hsp.al_start
            hsp_ref_pos = hsp.ref_start

            if hsp_read_pos < cig_read_pos:
                # print("overlap in read")
                while hsp_read_pos < cig_read_pos:
                    popped = cig_stack.pop(0)
                    if popped[1] == 'M':
                        hsp_ref_pos += 1
                    hsp_read_pos += 1

            elif hsp_read_pos > cig_read_pos:
                # print("gap in read")
                diff_read = hsp_read_pos - cig_read_pos
                ins = 'I' * diff_read
                cig_read_pos = hsp_read_pos
                cigarbox = cigarbox + ins

            if hsp_ref_pos < cig_ref_pos:
                # print("overlap in ref")
                while hsp_ref_pos < cig_ref_pos:
                    if len(cig_stack) is not 0:
                        cig_stack.pop(0)
                        last_clip = last_clip + 'S'
                        cig_read_pos += 1
                    hsp_ref_pos += 1

            elif hsp_ref_pos > cig_ref_pos:
                # print("gap in ref")
                diff_ref = hsp_ref_pos - cig_ref_pos
                skip = 'N' * diff_ref
                cig_ref_pos = hsp_ref_pos
                cigarbox = cigarbox + skip

            # then append the rest of the cigar_stack to the cigarbox
            while cig_stack:
                curr = cig_stack.pop(0)
                sym = curr[1]

                if sym == 'M':
                    cigarbox = cigarbox + sym
                    cig_read_pos += 1
                    cig_ref_pos += 1

                elif sym == 'I':
                    cigarbox = cigarbox + sym
                    cig_read_pos += 1

                elif sym == 'D':
                    cigarbox = cigarbox + sym
                    cig_ref_pos += 1

                else:
                    print("some other symbol in cigar")

            # print("segment done")

        cigarbox = cigarbox + last_clip
        # print("used all segments")

        # transform cigarbox to cigarstring
        cig_groups = itertools.groupby(cigarbox)
        cig_list = [(sum(1 for _ in count), label) for label, count in cig_groups]
        cig_str = [str(count) + str(symbol) for count, symbol in cig_list]
        cigar_string = ''.join(cig_str)

        self.cig = cigar_string
        return(cigar_string)

    def write_read(self, f):
        # replaces the reference_start and the cigarstring
        read_string = self.hsps[0].orig_container
        read_list = read_string.split('\t')

        ref_starts = [hsp.ref_start for hsp in self.hsps]
        read_list[3] = str(min(ref_starts) + 1)
        read_list[5] = self.cig
        
        if read_list[1] == '256':
            read_list[1] = '0'
        elif read_list[1] == '272':
            read_list[1] = '16'
            
        # check if len CIG == len seq
        cig_split = iter([''.join(x) for _, x in itertools.groupby(self.cig, key=str.isdigit)])
        cig_len = sum([int(x) for x in cig_split if next(cig_split) not in 'ND'])
        seq_len = len(read_list[9])
        
        if cig_len != seq_len:
            print('CIGAR and read of unequal length')
            print(self.cig)
            print(read_list)

        read_out = '\t'.join(read_list)
        # print(read_out)
        f.write(read_out + '\n')


class HSP:
    # translator for the AlignmentSegment object from pysam
    def __init__(self, cigartuples, al_start, al_end, ref_start, ref_end, orig_container):

        self.cigartuples = cigartuples
        self.al_start = al_start
        self.al_end = al_end
        self.al_len = self.al_end - self.al_start
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.read_range = set(range(self.al_start, self.al_end))
        self.ref_range = set(range(self.ref_start, self.ref_end))
        self.orig_container = orig_container


def get_ranges(inp):
    # input is a list of integers
    # returns a list of connected lists
    ranges = []
    for key, group in itertools.groupby(enumerate(sorted(inp)), lambda tup: tup[0] - tup[1]):
        group = list(group)
        glist = [x[1] for x in group]
        ranges.append(glist)
    return(ranges)


def nand(a, b):
    res = [x for x in a if x not in b]
    return(res)
    