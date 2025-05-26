import argparse
from pathlib import Path
from importlib.metadata import version

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

    logging.info(f"deviaTE {version('deviaTE')}")
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


def translate_name(name: str) -> str:
    """
    Translate a name by replacing invalid characters with dashes.
    :param name: The input name string.
    :return: The translated name string.
    """
    invalid_chars = '/\`*|;":. ' + "'"
    # make sure names don't contain any invalid characters
    if any(x in name for x in invalid_chars):
        for c in invalid_chars:
            name = name.replace(c, '-')
    return name.strip()



def find_blocks_generic(arr: NDArray, x: int, min_len: int) -> NDArray:
    """
    Find blocks in the array that match x.

    :param arr: The input array.
    :param x: The value to find blocks of.
    :param min_len: The minimum length of blocks to report.
    :return: An array containing the start and end positions of the blocks.
    """
    # find run starts
    x_pos = np.where(arr == x)[0]

    if x_pos.shape[0] == 0:
        return np.array([])

    # diff between neighboring loc
    x_diff = np.diff(x_pos)
    # if diff > 1: new block
    big_dist = np.where(x_diff > 1)[0]
    # the first entry is a block start and then all other where a big change happens
    # also each change is a block end, and the last entry of course as well
    block_ranges = np.concatenate((np.array([x_pos[0]]), x_pos[big_dist + 1],
                                   x_pos[big_dist] + 1, np.array([x_pos[-1] + 1])))
    blocks = block_ranges.reshape(big_dist.shape[0] + 1, 2, order='F')
    # only report blocks longer than min_len
    blocks_filt = blocks[np.where(blocks[:, 1] - blocks[:, 0] > min_len)[0], :]
    return blocks_filt




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





