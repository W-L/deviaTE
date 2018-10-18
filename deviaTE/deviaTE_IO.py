#!/usr/bin/env python3

import subprocess
import pysam
import os



_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    return os.path.join(_ROOT, path)


def execute(command):
    running = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               encoding='utf-8', shell=True)
    stdout, stderr = running.communicate()
    print(stdout)
    print(stderr)
    
    
def map_bwa(command, outfile):
    mapping = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               encoding='utf-8', shell=True)
    m = open(outfile, 'w')

    while True:
        chunk = mapping.stdout.read(4096)
        if len(chunk) is 0:
            m.close()
            break
        else:
            m.write(chunk)

    print(str(mapping.stderr))
        

class fq_file:
    def __init__(self, inp):
        self.path = inp
        
    def prep(self, lib, qual_tr, min_rl, min_al, read_ty, thr):
        args = ['deviaTE_prep',
                '--input', self.path,
                '--qual_threshold', qual_tr,
                '--min_read_length', min_rl,
                '--min_alignment_length', min_al,
                '--quality_encoding', read_ty,
                '--threads', thr]

        if lib:
            args.append('--library')
            args.append(lib)

        execute(command=' '.join(args))
    
        
class bam_file:
    def __init__(self, inp, orig_name=None):
        self.path = inp
        self.orig_name = orig_name
                    
    def fuse(self):
        args = ['deviaTE_fuse',
                '--input', self.path]
        
        execute(command=' '.join(args))
        
    def analyze(self, lib, fam, sid, out, anno, corr, hqt, scgs, rpm):
        args = ['deviaTE_analyse',
                '--input', self.path,
                '--family', fam,
                '--sample_id', sid,
                '--output', out,
                '--hq_threshold', hqt]
        
        if lib:
            args.append('--library')
            args.append(lib)
        if anno:
            args.append('--annotation')
            args.append(anno)
        if corr:
            args.append('--no_freq_corr')
        if scgs:
            args.append('--single_copy_genes')
            args.append(scgs)
        if rpm:
            args.append('--rpm')
            
        execute(command=' '.join(args))    
        
    
class analysis_table:
    
    def __init__(self, inp):
        self.path = inp
    
    def plot(self, out, free_y):
        args = ['deviaTE_plot',
                '--input', self.path,
                '--output', out]
        
        if free_y:
            args.append('--free_yaxis')
        
        execute(command=' '.join(args))    


def filter_alignment_length(inp, outp, lim):
    # remove reads under alignment length limit
    inpfile = pysam.AlignmentFile(inp, 'r')
    outfile = pysam.AlignmentFile(outp, 'w', template=inpfile)

    for read in inpfile.fetch():
        al_len = read.query_alignment_length
        if al_len >= int(lim):
            outfile.write(read)

    inpfile.close()
    outfile.close()




