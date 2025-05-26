from io import StringIO
from collections import defaultdict
from pathlib import Path

import numpy as np

from deviaTE.utils import translate_name



class PafLine:
    '''
    @DynamicAttrs
    parse a single alignment from a PAF into a flexible container
    '''

    def __init__(self, line: str, tags: bool = True):
        """
        Parse a line from a PAF file
        :param line: string representation of a line in a PAF file
        :param tags: boolean indicator whether to parse tags
        """
        self.line = line
        fields = ['qname', 'qlen', 'qstart', 'qend',
                  'strand', 'tname', 'tlen', 'tstart', 'tend',
                  'num_matches', 'alignment_block_length',
                  'mapq']
        core = 12
        record = line.strip().split("\t")
        # convert the fields to their actual types
        f = PafLine.format_records(record[:core])
        for i in range(core):
            setattr(self, fields[i], f[i])
        # make sure query and target name are strings
        self.qname = str(self.qname)
        self.tname = str(self.tname)
        self.tname = translate_name(self.tname)
        # marker for reverse sequences
        self.rev = 0 if self.strand == '+' else 1
        # parse the tags only if needed
        if tags:
            tags_parsed = PafLine.parse_tags(record[core:])
            self.align_score = int(tags_parsed.get("AS", 0))
            self.cigar = tags_parsed.get("cg", None)
            self.s1 = tags_parsed.get("s1", 0)
            prim = tags_parsed.get("tp", None)
            self.primary = 1 if prim == 'P' else 0
            self.seq = tags_parsed.get("sq", None)
            self.qual = tags_parsed.get("ql", None)


    @staticmethod
    def format_records(record: list) -> list:
        """
        Helper function to make fields of a PafLine the right type
        :param record: Split string of PAFline into list
        :return: Same list but with types converted to int
        """
        return [PafLine.conv_type(x, int) for x in record]


    @staticmethod
    def parse_tags(tags: list) -> dict:
        """
        Parse tags of a PAFline into a dictionary
        :param tags: List of SAM style tags
        :return: Dict of SAM style tags
        """
        c = {"i": int, "A": str, "f": float, "Z": str}
        return {
            key: PafLine.conv_type(val, c[tag])
            for key, tag, val in (x.split(":") for x in tags)
        }


    @staticmethod
    def conv_type(s: str, func: callable):
        """
        Generic converter, to change strings to other types
        :param s: Input string to convert to a different type
        :param func: Target type of input string
        :return: Either converted or original type
        """
        try:
            return func(s)
        except ValueError:
            return s





# shorthand typehint used in many places
paf_dict_type = dict[str, list[PafLine]]


class Paf:

    def __init__(self):
        pass


    @staticmethod
    def parse_PAF(paf_file: str | StringIO, min_len: int = 1) -> dict:
        """
        Parse the contents of a PAF file into a dictionary of records
        :param paf_file: Can be either a string or a StringIO object
        :param min_len: minimum alignment length to consider an entry
        :return: Dict of parsed PAF file
        """
        if isinstance(paf_file, str) and Path(paf_file).is_file():
            with open(paf_file, 'r') as paff:
                paf_dict = Paf._parse_content(fh=paff, min_len=min_len)
        elif isinstance(paf_file, StringIO):
            paf_dict = Paf._parse_content(fh=paf_file, min_len=min_len)
        else:
            raise ValueError("need file path or StringIO")
        return paf_dict


    @staticmethod
    def _parse_content(fh, min_len: int) -> dict:
        """
        Parser for PAF files into defaultdicts
        :param fh: Filehandle of PAF
        :param min_len: minimum alignment block length
        :return: parsed dict with PAF entries
        """
        paf_dict = defaultdict(list)
        # iterate records in the paf file
        for record in fh:
            paf = PafLine(record)
            # FILTERING of PAF ENTRIES
            if paf.alignment_block_length < min_len:
                continue
            if not paf.primary:
                continue
            # add this entry to paf dict
            paf_dict[str(paf.qname)].append(paf)
        return paf_dict


    @staticmethod
    def choose_best_mapper(records: list) -> list:
        """
        Structured array to decide between ties, by using the score of the DP algorithm
        :param records: List of multiple mappers to decide from
        :return: Best mapper according to attributes
        """
        mapq = [(record.mapq, record.align_score) for record in records]
        custom_dtypes = [('q', int), ('dp', int)]
        mapping_qualities = np.array(mapq, dtype=custom_dtypes)
        sorted_qual = np.argsort(mapping_qualities, order=["q", "dp"])
        record = [records[sorted_qual[-1]]]
        return record


    @staticmethod
    def single_rec(paf_dict: paf_dict_type) -> paf_dict_type:
        """
        Iterate a paf dict and find the best mapping if there are multiple
        :param paf_dict: Input paf dictionary
        :return: paf dict with one record per input read
        """
        for rid, rec_list in paf_dict.items():
            if len(rec_list) > 1:
                chosen_rec = Paf.choose_best_mapper(rec_list)
                paf_dict[rid] = chosen_rec
        return paf_dict




