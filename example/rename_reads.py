from pathlib import Path
import gzip
import sys


def _readfq(fp):
    last = None
    while True:
        if not last:
            for ll in fp:
                if ll[0] in ">@":
                    last = ll[:-1]
                    break
        if not last:
            break
        desc, name, seqs, last = last[1:], last[1:].partition(" ")[0], [], None
        for ll in fp:
            if ll[0] in "@+>":
                last = ll[:-1]
                break
            seqs.append(ll[:-1])
        if not last or last[0] != "+":
            yield desc, name, "".join(seqs), None
            if not last:
                break
        else:
            seq, leng, seqs = "".join(seqs), 0, []
            for ll in fp:
                seqs.append(ll[:-1])
                leng += len(ll) - 1
                if leng >= len(seq):
                    last = None
                    yield desc, name, seq, "".join(seqs)
                    break
            if last:
                yield desc, name, seq, None
                break




def modify_read_names(fastq_path):
    # check whether fastq is gzipped
    if fastq_path.name.endswith(('.gz', '.gzip')):
        fh = gzip.open(fastq_path, 'rt')
    else:
        fh = open(fastq_path, 'rt')

    i = 0
    for desc, name, seq, qual in _readfq(fh):
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
