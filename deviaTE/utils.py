import argparse
from pathlib import Path

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



def is_gzipped(filepath: str | Path) -> bool:
    """
    Check if a file is gzipped.
    :param filepath: string of filepath to be checked.
    :return: boolean indicating if the file is gzipped.
    """
    with open(filepath, 'rb') as fh:
        return fh.read(2) == b'\x1f\x8b'



def rapidgzip_count_lines(filepath: str | Path) -> int:
    '''
    Count the number of lines in a file.
    :param filepath:  string or path to a file.
    :return:
    '''
    import rapidgzip
    result = 0
    with rapidgzip.open(str(filepath)) as file:
        while chunk := file.read(1024 * 1024):
            result += chunk.count(b'\n')
    return result


def rawcount(filename: str | Path) -> int:
    '''
    Get the number of lines in a file.
    :param filename: string or path to a file.
    :return: number of lines in the file.
    '''
    sf = open(str(filename), 'rb')

    lines = 0
    buf_size = 1024 * 1024
    read_f = sf.raw.read

    buf = read_f(buf_size)
    while buf:
        lines += buf.count(b'\n')
        buf = read_f(buf_size)

    sf.close()
    return int(lines)



class QualTrans:

    def __init__(self):
        '''
        translator class to shift the colon quality symbol
        this makes sure that the paf tag parsing still works
        '''
        self.qtrans = str.maketrans({':': '9'})


    def shift_qual(self, qual_string: str) -> str:
        '''
        apply the translation of the quality string
        :param qual_string: input quality string
        :return: translated quality string
        '''
        return qual_string.translate(self.qtrans)





