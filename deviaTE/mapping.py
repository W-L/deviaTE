from pathlib import Path

import mappy

from deviaTE.utils import readfq


class Mapper:

    def __init__(self, ref: str, preset: str):
        """
        Initialize a Mapper object; wrapper for minimap2's mappy implementation
        For the default case of mapping against some linear references
        :param ref: The path to the reference FASTA file.
        :param preset: Parameter preset for mm2 for read type
        """
        self.preset = preset
        # check that the given reference exists
        if not Path(ref).is_file():
            raise FileNotFoundError("Given reference file does not exist")
        valid_ref_types = ['.mmi', '.fa', '.fasta', '.fa.gz', '.fasta.gz']
        if not Path(ref).suffix in valid_ref_types:
            raise ValueError(f"file type of the reference file {ref} must be one of {valid_ref_types}")
        # initialise the mappy aligner
        self.aligner = mappy.Aligner(fn_idx_in=ref, preset=preset)



    def map_file(self, seq_file: str) -> str:
        """
        Map a full batch of reads and return the PAF formatted mapping hits.
        :param seq_file: path to a file with sequences
        :return: tuple of output file name and number of input reads
        """
        if not Path(seq_file).is_file():
            raise FileNotFoundError("Given sequence file does not exist")
        # collect results in a list
        results = []
        with open(seq_file, 'r') as sf:
            for desc, name, seq, quals in readfq(sf):
                hits = self.aligner.map(seq)
                for hit in hits:
                    results.append(f"{name}\t{len(seq)}\t{hit}")
        # transform to a single string
        alignments = '\n'.join(results)
        # write alignments to file
        outfile = f"{Path(seq_file).name}.paf"
        with open(outfile, 'w') as paf_out:
            paf_out.write(alignments)
        return outfile


