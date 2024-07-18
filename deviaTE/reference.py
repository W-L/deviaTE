import logging
import gzip
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
import pandas as pd

from deviaTE.utils import Seq2Int, readfq, is_gzipped
from deviaTE.analyse import CoverageConverter, incr_type
from deviaTE.plot import visualise
from deviaTE.config import Config
import deviaTE.mapping
import deviaTE.paf



class InputFile:

    def __init__(self, conf: Config, infile: Path):
        """
        Initialise object representing a single input file
        :param conf: Configuration object of the experiment
        :param infile: file that is analysed in this class
        """
        self.conf = conf
        self.input = infile

        logging.info(f"starting analysis of input file: {self.input}")

        # initialise a mapper given the configuration
        self.mapper = deviaTE.mapping.Mapper(ref=conf.args.library, preset=conf.args.preset)


    def analyse_coverage(self) -> None:
        """
        Run analysis for this input file and all selected families
        :return:
        """
        # only map if the alignment file does not already exist
        afile = f'{self.input.name}.paf'
        if not Path(afile).is_file():
            # map all reads against the reference library
            logging.info(f"Starting mapping of {self.input.name} to {self.conf.args.library}")
            alignment_file = self.mapper.map_file(seq_file=str(self.input))
        else:
            logging.info(f'Alignment file {afile} already exists. Skipping read mapping')
            alignment_file = afile
        # filter the alignments for length
        logging.info(f"Parsing alignment file: {alignment_file}")
        paf_dict = deviaTE.paf.Paf.parse_PAF(alignment_file, min_len=self.conf.args.min_align_len)
        # if multiple mappings for a read: choose best one
        self.paf_dict = deviaTE.paf.Paf.single_rec(paf_dict)
        # convert the coverage - need to load the reads into memory for this
        self.incr, self.incr_hq = self._convert_to_increments()
        # if selected, get norm factor of single copy genes
        if self.conf.args.single_copy_genes:
            self.scg_norm_fac = single_gene_norm_fac(scgs=self.conf.args.single_copy_genes,
                                                sequence_lib=self.conf.sequences,
                                                increments=self.incr)
        else:
            self.scg_norm_fac = 0


    def analyse_families(self) -> None:
        """
        Run the analysis of all selected families after analysing coverage
        :return:
        """
        # run further analysis for all selected references
        self.results_files = []
        logging.info(f"Starting analysis of TE sequeces")
        for family in self.conf.args.families:
            rf = self._analyse_family(f=family)
            self.results_files.append(rf)


    def _convert_to_increments(self, hq_limit: int = 15) -> tuple[incr_type, incr_type]:
        """
        Convert the mappings to coverage increments of the reference sequences
        :param hq_limit: Limit at which to consider coverage as high quality
        :return: tuple of dicts of coverage counts per reference
        """

        # use converted object to get the counts
        c = CoverageConverter()
        c_hq = CoverageConverter(qt=hq_limit)
        incr = c.convert_records(paf_dict=self.paf_dict, seqs=self.seqs, quals=self.quals)
        incr_hq = c_hq.convert_records(paf_dict=self.paf_dict, seqs=self.seqs, quals=self.quals)
        return incr, incr_hq


    def _analyse_family(self, f: str) -> str:
        """
        Analyse single family
        :param f: Identifier of the reference sequence
        :return:
        """
        # create a reference object
        fam = Reference(
            sid=f,
            seq=self.conf.sequences[f],
            incr=self.incr[f],
            incr_hq=self.incr_hq[f]
        )
        # normalize by rpm
        if self.conf.args.rpm:
            logging.info('Normalization: reads per million')
            fam.normalize_rpm(nreads=len(self.seqs))
        # or normalize with set of single copy genes
        elif self.conf.args.single_copy_genes:
            logging.info('Normalization: single copy genes')
            fam.normalize_scg(scg_factor=self.scg_norm_fac)
        # find polymorphism and non-reference sites
        fam.mark_snps()
        # put results in a table
        results_file = fam.write_results(infile=self.input)
        return results_file


    def visualise(self) -> None:
        """
        Produce visualised results
        :return:
        """
        logging.info(f"Starting visualisation of TE sequences")
        for rf in self.results_files:
            visualise(rf, self.conf.args.annotation)



class Reference:

    def __init__(self, sid: str, seq: str, incr: list = None, incr_hq: list = None):
        """
        This is a class to represent reference sequences, ie a family to be analysed
        :param sid: Identifier of the reference
        :param seq: Sequence of the reference
        :param incr: coverage increments after mapping and converting
        :param incr_hq: high-quality coverage increments after mapping and converting
        """
        self.sid = sid
        self.seq = seq
        self.length = len(seq)
        self.cov = np.zeros(shape=(self.length, 5), dtype=int)
        self.cov_hq = np.zeros(shape=(self.length, 5), dtype=int)
        self.ihat = 'NA'
        self.ihat_hq = 'NA'
        self.snp = np.zeros(shape=self.length, dtype=bool)
        self.annotations = []
        if incr:
            self._update_cov(c="cov", incr=incr)
        if incr_hq:
            self._update_cov(c="cov_hq", incr=incr_hq)


    def _update_cov(self, c: str, incr: list) -> None:
        """
        Apply the counts of coverage to the attributes of this class
        :param c: which coverage to increment
        :param incr: List of converted increments
        :return:
        """
        if len(incr) < 1:
            return
        # temporary container for coverage
        tmp_cov = np.zeros(shape=self.cov.shape, dtype="uint16")
        for (start, end, query_arr, addition) in incr:
            # range to index into the bit we want to increment
            indices = np.arange(query_arr.shape[0])
            # add the "addition" to the corresponding bases
            np.add.at(tmp_cov[start: end], (indices, query_arr), addition)
        # add the coverage to the previous one
        setattr(self, c, getattr(self, c) + tmp_cov)


    def _cov_sum(self) -> NDArray:
        """
        Sum of coverage counts for all nucleotides and gaps
        :return: NDArray with coverage counts
        """
        return np.sum(self.cov, axis=1)


    def _cov_sum_hq(self) -> NDArray:
        """
        Sum of high-quality coverage counts
        :return: NDArray with high-quality coverage counts
        """
        return np.sum(self.cov_hq, axis=1)


    def mean_cov(self) -> float:
        """
        Mean coverage of this reference sequence
        :return: Mean coverage after summing across all nucleotides
        """
        return np.mean(self._cov_sum())


    def _normalize(self, norm_factor: float) -> None:
        """
        Normalise the coverage counts by given factor
        :param norm_factor: Can be either from rpm or scg
        :return:
        """
        # once we normalise coverage can't be ints
        self.cov = self.cov.astype(float)
        self.cov_hq = self.cov.astype(float)
        self.cov /= norm_factor
        self.cov_hq /= norm_factor


    def normalize_rpm(self, nreads: int) -> None:
        """
        Normalise by reads per million
        :param nreads: Number of total reads of the concerned input file
        :return:
        """
        rpm_fac = nreads / (10 ** 6)
        if not rpm_fac:
            logging.info("Normalisation factor is 0 or NA. Skipping normalisation and estimation of insertions")
        else:
            self._normalize(norm_factor=rpm_fac)


    def normalize_scg(self, scg_factor: float) -> None:
        """
        Normalise using the single copy genes, and estimate the number of insertions
        :param scg_factor: Mean coverage of the single copy genes
        :return:
        """
        if not scg_factor:
            logging.info("Normalisation factor is 0 or NA. Skipping normalisation and estimation of insertions")
        else:
            self._estimate_insertions(norm_factor=scg_factor)
            self._normalize(norm_factor=scg_factor)


    def _estimate_insertions(self, norm_factor: float) -> None:
        """
        Estimate the number of insertions of this reference
        :param norm_factor: Mean coverage of the single copy genes
        :return:
        """
        self.ihat = self.mean_cov() / norm_factor
        self.ihat_hq = np.mean(self._cov_sum_hq()) / norm_factor


    def mark_snps(self) -> None:
        """
        Determine which sites show variation
        :return:
        """
        covsum = self._cov_sum()
        meanc = self.mean_cov()
        cov_freq = np.divide(self.cov.astype(float), covsum[:, np.newaxis].astype(float),
                             out=np.zeros_like(self.cov.astype(float)),
                             where=covsum[:, np.newaxis] != 0)
        # polymorphic SNP -> 3 conditions
        # coverage is higher than most abundant base
        # minimum count of meancov * 0.1
        # minimum frequency of 0.1
        polys = np.where((covsum > np.max(self.cov, axis=1)) & (covsum > meanc * 0.1) & (np.max(cov_freq, axis=1) < 0.9))[0]
        self.snp[polys] = True
        # alternative option: the majority of counts are for nucleotide that is not reference
        intseq = Seq2Int(seq=self.seq).translate()
        notref = np.where(np.argmax(self.cov, axis=1) != intseq)[0]
        self.snp[notref] = True


    def write_results(self, infile: Path) -> str:
        """
        Write the results to a file
        :param infile: Input file that was analysed
        :return: name of the output file to pass onto visualisation
        """
        # put results together into a dataframe
        frame = {
            '#TEfam': [self.sid] * self.length,
            'sample_id': [infile.name] * self.length,
            'pos': list(range(self.length)),
            'refbase': list(self.seq),
            'A': self.cov[:, 0],
            'C': self.cov[:, 1],
            'G': self.cov[:, 2],
            'T': self.cov[:, 3],
            'cov': np.sum(self.cov, axis=1),
            'hq_cov': np.sum(self.cov_hq, axis=1),
            'snp': self.snp,
            'delet': self.cov[:, 4]
        }
        df = pd.DataFrame(frame)
        # add headers
        ihat = f'{self.ihat} or {self.ihat_hq} (hq coverage only)'
        outf = f'{infile.name}.{self.sid}.deviate'
        with open(outf, 'w') as outfile:
            outfile.write('# insertions/haploid: ' + ihat + '\n')
        df.to_csv(outf, index=False, sep=' ', mode='a')
        return outf



def single_gene_norm_fac(scgs: list, sequence_lib: dict[str, str], increments: incr_type) -> float:
    """
    Calculate the normalisation factor (mean coverage) of the single copy genes
    :param scgs: List of names of the genes that should be used
    :param sequence_lib: Sequences in the reference library file as dict
    :param increments: Converted coverage counts
    :return: Normalisation factor
    """
    scg_means = []
    for s in scgs:
        scg = Reference(sid=s, seq=sequence_lib[s], incr=increments[s])
        scg_means.append(scg.mean_cov())
    scg_norm_fac = np.mean(scg_means)
    return scg_norm_fac


