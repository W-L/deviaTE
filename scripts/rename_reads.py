from pathlib import Path
import gzip
import sys

from deviaTE.utils import readfq



def modify_read_names(fastq_path):
    # check whether fastq is gzipped
    if fastq_path.name.endswith(('.gz', '.gzip')):
        fh = gzip.open(fastq_path, 'rt')
    else:
        fh = open(fastq_path, 'rt')

    i = 0
    for desc, name, seq, qual in readfq(fh):
        name = desc.split(' ')[0]
        name_num = f'@{name}_{i}'
        print(name_num)
        print(seq)
        print('+')
        print(qual)
        i += 1

    fh.close()


fastq_file = Path(sys.argv[1])
modify_read_names(fastq_file)

