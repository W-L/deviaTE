import argparse
from typing import TextIO
from itertools import groupby

import numpy as np
from numpy.typing import NDArray


class Seq2Int:

    def __init__(self, seq: str):
        """
        Translator class to turn nucleotide sequences into integer arrays
        :param seq: input sequence
        """
        self.seq = seq
        transDict = {'A': '0', 'C': '1', 'G': '2', 'T': '3', 'N': '4'}
        self.base2int = str.maketrans(transDict)


    def translate(self) -> NDArray:
        """
        :return: translated nucleotide sequence
        """
        read_integer = self.seq.translate(self.base2int)
        int_seq = np.frombuffer(read_integer.encode(), 'u1') - ord('0')
        return int_seq


def init_logger(logfile: str, args: argparse.Namespace) -> None:
    """
    Initialize the logger with the given logfile and log the arguments.

    :param logfile: The path to the logfile.
    :param args: The arguments to log.
    """
    with open(logfile, 'w'):
        pass
    import logging
    logging.basicConfig(format='%(asctime)s %(message)s',
                        level=logging.INFO,
                        handlers=[logging.FileHandler(f"{logfile}"), logging.StreamHandler()])

    logging.info("deviaTE")
    logging.info('\n')
    for a, aval in args.__dict__.items():
        logging.info(f'{a} {aval}')
    logging.info('\n')



def reverse_complement(dna: str) -> str:
    '''
    Return the reverse complement of a dna string.
    :param dna: string of characters of the usual alphabet
    :return: reverse complemented input dna
    '''
    trans = str.maketrans('ATGC', 'TACG')
    rev_comp = dna.translate(trans)[::-1]
    return rev_comp



def read_fa(fh: TextIO):
    """
    Generator for fasta files: yields all headers and sequences in the file.

    :param fh: The file handle.
    :yield: A tuple containing the header and sequence.
    """
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        headerStr = header.__next__().strip().split(' ')[0]  # drop the ">"
        # join all sequence lines to one
        seq = "".join(s.strip() for s in faiter.__next__())
        yield headerStr, seq



def readfq(fp: TextIO):
    """
    Read a fastq file and yield the entries.

    :param fp: File handle for the fastq file.
    :yield: A tuple containing the fastq read header, read ID, and sequence.
    """
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for ll in fp:  # search for the start of the next record
                if ll[0] in ">@":  # fasta/q header line
                    last = ll[: -1]  # save this line
                    break
        if not last:
            break
        desc, name, seqs, last = last[1:], last[1:].partition(" ")[0], [], None
        for ll in fp:  # read the sequence
            if ll[0] in "@+>":
                last = ll[: -1]
                break
            seqs.append(ll[: -1])
        if not last or last[0] != "+":  # this is a fasta record
            yield desc, name, "".join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = "".join(seqs), 0, []
            for ll in fp:  # read the quality
                seqs.append(ll[: -1])
                leng += len(ll) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield desc, name, seq, "".join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield desc, name, seq, None  # yield a fasta record instead
                break




