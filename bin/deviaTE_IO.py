#! /usr/local/bin/python3
import subprocess


def execute(command):
    running = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               encoding='utf-8', shell=True)
    stdout, stderr = running.communicate()
    print(stdout)
    print(stderr)
        
        
class fq_file:
    def __init__(self, inp):
        self.path = inp
        
    def prep(self, script_loc, lib, qual_tr, min_rl, min_al, read_ty, thr):
        args = [script_loc + '/deviaTE_prep.sh',
                '-i', self.path,
                '-q', qual_tr,
                '-r', min_rl,
                '-a', min_al,
                '-y', read_ty,
                '-t', thr]

        if lib:
            args.append('-l')
            args.append(lib)

        execute(command=' '.join(args))
    
        
class bam_file:
    def __init__(self, inp, from_fq):
        self.path = inp
        self.from_fq = from_fq
        
        if from_fq:
            self.log = self.path.split('.')[0] + '.fastq.log'
            
    def fuse(self, script_loc):
        args = [script_loc + '/bin/deviaTE_fuse.py',
                '--input', self.path,
                '--si']
        
        execute(command=' '.join(args))
        
    def analyze(self, script_loc, lib, fam, sid, out, anno):
        args = [script_loc + '/deviaTE_analyse.py',
                '--input', self.path,
                '--family', fam,
                '--sample_id', sid,
                '--output', out]
        
        if lib:
            args.append('--library')
            args.append(lib)
        if anno:
            args.append('--annotation')
            args.append(anno)
        if self.from_fq:
            args.append('--log')
            args.append(self.log)
        
        execute(command=' '.join(args))    
        
    
class analysis_table:
    
    def __init__(self, inp):
        self.path = inp
    
    def plot(self, script_loc, out, free_y, col_ref):
        args = [script_loc + '/deviaTE_plot.R',
                '--input', self.path,
                '--output', out]
        
        if free_y:
            args.append('--free_yaxis')
        if col_ref:
            args.append('--color_reference')
        
        execute(command=' '.join(args))    

            