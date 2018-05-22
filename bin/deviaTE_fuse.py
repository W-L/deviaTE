#! /usr/local/bin/python3

import pysam
import sys
import argparse
import subprocess
import deviaTE_multiHSP as multiHSP
from collections import Counter, defaultdict


parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, required=True)
parser.add_argument('--output', type=str)
parser.add_argument('--si', action='store_true')
args = parser.parse_args()

if args.output is None:
    args.output = args.input + '.fused'

inpfile = pysam.AlignmentFile(args.input, 'r')
outfile = pysam.AlignmentFile(args.output, 'w', template=inpfile)
outfile.close()
outfile = open(args.output, 'a+')

read_dump = set()

for read in inpfile.fetch():
    rid = read.query_name
    if rid not in read_dump:
        read_dump.add(rid)
        seg_list = []
        inp_search = pysam.AlignmentFile(args.input, 'r')
        for read_search in inp_search.fetch():
            if rid == read_search.query_name:
                seg_list.append(read_search)
        
        if len(seg_list) > 1:
            #print(seg_list[0].query_name)
            # process the seg_list of this rid
            Multihit_list = []
            fams = defaultdict(list)
            
            for seg in seg_list:
                if seg.is_reverse:
                    fam_strand = seg.reference_name + '-'
                else:
                    fam_strand = seg.reference_name + '+'

                # for every family represented in this rid
                # and separated by strand
                # HSP is a translator-class for pysams AlignmentSegment
                fams[fam_strand].append(multiHSP.HSP(cigartuples=seg.cigartuples,
                                                     al_start=seg.query_alignment_start,
                                                     al_end=seg.query_alignment_end,
                                                     ref_start=seg.reference_start,
                                                     ref_end=seg.reference_end))
                
            # Multihit_list contains all Multihits per family/strand for the current rid
            for fam, hsp_list in fams.items():
                Multihit_list.append(multiHSP.Multihit(read_id=rid, hsp_list=hsp_list,
                                                       orig_container=seg_list[0].tostring(inpfile),
                                                       fam=fam))
            
            for mh in Multihit_list:
                mh.create_MACs()
                
                for mac in mh.MACs:
                    mac.construct()
                    mac.check_overlap(limit=5)
                    mac.check_distance(limit=20)
                    mac.score_MAC()
                    
                mh.find_hMAC()
                
            # sort list by hMAC-score - only keep first
            Multihit_list.sort(key=lambda h: h.hMAC_score, reverse=True)
            hMAC_read = Multihit_list[0].hMAC
            
            # construct new cigar
            hMAC_read.build_cigar()
            hMAC_read.write_read(f=outfile)
            
        else:
            outfile.write(read.tostring(inpfile) + '\n')
            
        inp_search.close()
outfile.close()
inpfile.close()


if args.si:
    # sort and index
    deviaTE = '/'.join(sys.argv[0].split('/')[:-2])
    
    sam = subprocess.Popen('source ' + deviaTE + '/CONFIG', shell=True, stdout=subprocess.PIPE, encoding='utf-8')
    stdout = sam.communicate()
    samtools = stdout[0].split('\n')[0]
    
    bam = [samtools, 'view', '-b', args.output, '-o', args.output + '.bam']
    sort = [samtools, 'sort', args.output + '.bam', '-o', args.output + '.sort.bam']
    index = [samtools, 'index', args.output + '.sort.bam', args.output + '.sort.bam.bai']
    
    subprocess.run(' '.join(bam), shell=True)
    subprocess.run(' '.join(sort), shell=True)
    subprocess.run(' '.join(index), shell=True)

    if '*' not in args.output:
        rm = ['rm', args.output, args.output + '.bam']
        subprocess.run(' '.join(rm), shell=True)


