# deviaTE

deviaTE is a python tool for the analysis and visualization of mobile genetic element sequences.

## Update v2 (2024-04)

Since the previous python base reached end of life, deviaTE needed an update. This update became quite substantial - hence bumped to version 2:

- changed dependencies: now purely python
- installation via pip only: no conda environment needed
- biggest dependency change: mapping via minimap2 instead of bwa (much faster and allows for modern long reads)
- arguments and usage simplified (see new instructions below)
- full coverage with pytest, docstrings, typehints, ...


Feature deprecation:

- Detection of structural variants is currently not included. Please use the previous version (at this [github link](https://github.com/W-L/deviaTE/tree/10d2b7063b2fef7fcaa24b0a45fa655a0c4d7565)) for this and let me know via GitHub issue that this is of interest, then I'll include it in the next version


## Installation

deviaTE needs python >=3.10 and pip:

`pip install deviaTE`



## Usage 


```
usage: deviaTE [-h] [--input INPUT] [--preset {sr,map-ont,map-pb,map-hifi}] [--library LIBRARY] [--annotation ANNOTATION] [--min_align_len MIN_ALIGN_LEN] [--families [FAMILIES ...]]
[--rpm | --single_copy_genes [SINGLE_COPY_GENES ...]]

options:
-h, --help                  show this help message and exit
--input INPUT               Input file(s) to be analysed. Can be *.fastq, *.fa, or directory of files. Optionally gzipped.
--preset {sr,map-ont,map-pb,map-hifi}   Minimap2 mapping preset. (sr, map-ont, map-pb, map-hifi) [sr]
--library LIBRARY           Path to reference library. Defaults to drosophila transposons from https://github.com/bergmanlab/drosophila-transposons
--annotation ANNOTATION     Path to annotation (gff) of sequences in library. Defaults to drosophila TE annotation from https://github.com/bergmanlab/drosophila-transposons
--min_align_len MIN_ALIGN_LEN           Minimum length of valid alignments
--families [FAMILIES ...]   Which transposon families to analyse. Default: all sequences in library.
--rpm                       normalize all abundances by reads per million
--single_copy_genes [SINGLE_COPY_GENES ...]     space-separated names of single-copy genes in reference to use for normalisation
```


DeviaTE is a command-line program that analyzes and visualizes the diversity of mobile genetic elements from sequencing data without the need for an assembled genome of the host species. 
The only required argument is `--input`. For this, it takes sequencing data (`--input` single file or directory of files). It can be used with short and long reads (`--preset`, minimap2 parameter preset for short reads [sr], nanopore reads [map-ont] or pacbio [map-pb, map-hifi]).
It also requires mobile genetic element consensus sequences (`--library`, fasta file). If no library is given it will use the Drosphila transposon sequences from https://github.com/bergmanlab/drosophila-transposons.
TEs to be analysed are selected with `--families`. These can be multiple (space-separated) or if not specified, all reference sequences in the library are used.


Available arguments are listed with `-h` or `--help`.


## Example

An example is available for testing. The sequences are from the Drosophila 12 Genomes Consortium et al. 2007. Evolution of genes and genomes on the Drosophila phylogeny. *Nature*. 450(7167):203-218.

We can analyse the TE jockey (DMLINEJA) and get a visualization using:

`deviaTE --input ../data/jockey_dmel.fastq --families FBte0000088`

this produces an alignment file called `jockey_dmel.fastq.paf`, creates the output table `jockey_dmel.fastq.FBte0000088.deviate` with information about coverage and estimated insertions (if selected), and the visualisation `jockey_dmel.fastq.FBte0000088.deviate.pdf`. 



Manual and Walkthrough of previous versions can be found (at this [github link](https://github.com/W-L/deviaTE/tree/10d2b7063b2fef7fcaa24b0a45fa655a0c4d7565))


## Description of results table

The table starts with some header-lines denoted by #. This header contains the estimated number of TE insertions (if selected) and column names. Each row corresponds to one position of the TE sequence. Since version 2, `hq_cov` reports coverage of high quality bases instead of high-quality mappings, since that's more interesing e.g. for nanopore data. 


| Column      | Description                                |
|-------------|--------------------------------------------|
| `TEfam`     | Name of the analysed TE family             |
| `sample_id` | input file name                            |
| `pos`       | position in the reference sequence         |
| `refbase`   | Nucleotide in the reference sequence at this position |
| `A C G T`   | counts of each nucleotide at this position |
| `cov`       | total coverage at this position            |
| `hq_cov`    | coverage of high-quality bases only (>Q15) |
| `snp`       | indicator for variant position             |
| `delet`     | count of gap observations                  |





## Normalization methods

By default no normalization is performed and reported counts are raw abundances, which are not suitable for comparing TEs between samples. Therefore two different strategies are implemented, normalization per million mapped reads and normalization by single-copy genes.

- Per million mapped reads: Normalize all counts per million mapped reads to account for different sequencing depths, selected by `--rpm`.
- Single-copy gene normalization: relate all counts to the number of reads mapping to one or more single-copy genes. Accounts for differential sequencing depth and additionally estimates the insertion/copy number of the analysed mobile genetic element. For this normalization method add the sequence of multiple single-copy genes to the file containing the TE consensus sequences used as `--library`. Then add `--single_copy_genes GENE1 GENE2 GENE3 ...`, where GENE1 etc. are the headers in the library file. The estimated copy number per haploid genome is written to the header-section of the resulting output table.



## Special use-case: *Drosophila*

If you are analyzing TEs in *Drosophila* specifying a `--library` or `--annotation` of reference sequences is optional. By default deviaTE automatically downloads and uses the TE library from https://github.com/bergmanlab/drosophila-transposons if no library and annotation are given.

For single-copy gene normalization in Drosophila five genes are automatically added to the library (Dmel_rpl32, Dmel_piwi, Dmel_Act5C, Dmel_RpII140 and Dmel_p53), which can be used for normalisation: 

`--single_copy_genes Dmel_rpl32 Dmel_piwi ...`


## Special use-case: Paired-end reads 

You can use DeviaTE for paired-end reads by mapping them in single read mode. This can be done, for example, by using a single concatenated fastq file that contains both read pairs (read1 and read2). In order to prevent some issues the reads then need to be given unique names (e.g. readID_1, readID_2, ..., readID_n), which can be achieved using a script found at: `scripts/rename_reads.py`  (Thanks Anna for pointing out this issue)

```
python rename_reads.py sample.fq >sample_uniq.fq
```


## Citation

A paper describing deviaTE is available here: https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13030

```
@article{weilguny2019,
  title = {{{DeviaTE}}: {{Assembly-free}} Analysis and Visualization of Mobile Genetic Element Composition},
  author = {Weilguny, Lukas and Kofler, Robert},
  year = {2019},
  journal = {Molecular Ecology Resources},
  volume = {19},
  number = {5},
  pages = {1346--1354},
  doi = {10.1111/1755-0998.13030}
}
```

## Bugs, issues, questions, suggestions ...
If you find any problems, have questions or ideas for further improvement please use the issue tracker on this repository, thanks!


## License
deviaTE is licensed under the GPLv3 License


## Testing

The code is covered by pytests. To run these install: `pip install pytest pytest-cov`. Then run tests: `cd tests; pytest --cov --cov-report html`.
To test local builds: `pip install dist/deviate-2.0.0-py3-none-any.whl --force-reinstall --no-deps`

