# Walkthrough for deviaTE

## Getting started

Right after installing and skipping through the manual we recommend to test that the tool works properly. 
We have prepared a FASTQ file containing only reads of *Drosophila melanogaster* that map to the transposable element *jockey* 
(can be downloaded [here](https://github.com/W-L/deviaTE/blob/master/example/jockey_dmel.fastq) [~3Mb, in the example folder of this repository]).
These sequences are sanger reads from Drosophila 12 Genomes Consortium et al. 2007. Evolution of genes and genomes on the Drosophila phylogeny. *Nature*. 450(7167):203-218.

You can now use the convenient wrapper script `deviaTE` to analyse the TE jockey (DMLINEJA) and get a visualization:

```bash
deviaTE --input_fq jockey_dmel.fastq --read_type sanger --families DMLINEJA
```

Alternatively you can run all three steps of deviaTE manually with: 

```bash
deviaTE_prep --input jockey_dmel.fastq --quality_encoding sanger
deviaTE_analyse --input jockey_dmel.fastq.fused.sort.bam --family DMLINEJA
deviaTE_plot --input jockey_dmel.fastq.fused.sort.bam.DMLINEJA
```

this will first trim, map and filter the reads, which yields an alignment file called `jockey_dmel.fastq.fused.sort.bam`
then the tool collects the quantitative information and calculates estimators, which produces the output table `jockey_dmel.fastq.fused.sort.bam.DMLINEJA`.
The name is simply a combination of the input to the tool and the analysed TE family.
Finally a visualization is produced, resulting in an illustration of jockey that should look like [this](https://github.com/W-L/deviaTE/blob/master/example/jockey_dmel.fastq.DMLINEJA.pdf)





## Basic usage

<p align="center">
  <img src="https://github.com/W-L/deviaTE/blob/master/doc/workflow.png" alt="Architecture of deviaTE" width="600"/>
</p>



### Using the single-command wrapper script

To analyze and visualize a transposable element families in from sequencing reads, you can use

```deviaTE --input_fq foo.fastq --families TEfamily1,TEfamily2,... --library TE_consensus_sequences.fasta```

where foo.fastq contains sequencing reads and TEfamily1 etc. are headers of TE reference sequences in the file defined by `--library`. Any fasta file used as library must first be indexed with:

```
bwa index TE_consensus_sequences.fasta
```

Multiple families can be selected, they just need to be separated by commas without spaces. DeviaTE can also be applied to multiple fastq files in a directory using `--input_fq_dir` when running the program from within that directory. 

Other available arguments can be seen with `deviaTE -h/--help` and are documented in the [Manual](https://github.com/W-L/deviaTE/blob/master/doc/MANUAL.md) 

This will result in a table of quantitative information and a visualization for each selected TE families in each of the specified samples.



### Using the three sequential steps

Similar to the example presented above, you can manually run the three steps in the workflow of DeviaTE for more flexibility. The first step involves the script called `deviaTE_prep`, which trims, maps and filters the input sequencing reads and performs the detection of internally deleted variants. Basic usage involves specifying the input sequences (single sample) and the consensus sequences of TEs, e. g.:

```
deviaTE_prep --input foo.fastq --library TE_consensus_sequences.fasta
```

Further arguments are available, which specify the quality encoding and the parameters used for trimming and filerting.

From this step you obtain an alignment file of the form `foo.fastq.fused.sort.bam` and an index of this file. For preparation no TE family needs to be selected, thus reads are mapped to all reference sequences present in the library. This means that the next step of analysing TE families can be performed multiple times on this alignment file. The basic command for analysing a TE family is of this form:

```
deviaTE_analyse --input foo.fastq.fused.sort.bam --family TEfamily --library TE_consensus_sequences.fasta
```

Where the input is the alignment file from the previous step, the library contains the reference sequences and the selected family is a header within that library. Additional arguments can be used to provide an annotation of the TE sequences, name the sample and the output, set a threshold for unambiguously mapped reads or to select a normalization method.

This script produces the output table containing the quantitative information, which by default will be named as: `input.TEfamily`. This file is then used to produce a visualization, e. g. in basic scenarios:

```
deviaTE_plot --input foo.fastq.fused.sort.bam.TEfamily
```

Here, optional arguments can be used to change the output name and its format (pdf or eps), change fontsize or switch the y-axis to be free between multiple plots in a grid, e. g. when more than one sample is plotted and the coverage between them is highly different. For a description on how to plot multiple samples and TE families at once, see below.


## Description of the output table

Besides the visualization an output table is produced, which contains the quantitative information about the TE family. Generally every row in the table corresponds to one position of the TE sequence.
The table starts with some header-lines denoted by #. The first line contains a timestamp and the command used to generate the file, while the following line states the estimated number of TE insertions (only when single gene normalization was selected, see below). The final line of the header section contains the column names of the table. 

The following is a list of these columns in the output table of deviaTE:

column name | example value | description
--- | --- | ---
`TEfam` | DMLINEJA | Name of the analysed TE family
`sample_id` | jockey_dmel.fastq | either the input file name or an identifier provided by `--sample_id`
`pos` | 709 | position in the reference sequence
`refbase` | T | base in the reference at this position
`A C G T` | 124 | counts of each nucleotide at this position
`cov` | 140 | total coverage at this position
`phys_cov` | 18 | physical coverage at this position, e.g. spanned by deletion
`hq_cov` | 140 | coverage above the specified threshold for unambigously mapped reads (`--hq_threshold`)
`snp` | True | denotes a single nucleotide polymorphism at this position
`refsnp` | False | denotes a fixed difference to the reference sequence
`int_del` | 709:771:25 | describes an internal deletion at this position in the form of `start:end:abundance`
`int_del_freq` | 709:771:0.4 | frequency of internal deletions in the form of `start:end:frequency` (e.g. 0.4 means 40% of reads have this deletion)
`trunc_left` | 17 | absolute number of reads that show a truncation to the left of this position
`trunc_right` | 17 | absolute number of reads that show a truncation to the right of this position
`ins` | 709:711:2 | describes an insertion of two bp length at this position, occuring twice. `start:end:abundance`
`delet` | 709:710:3 | similarly for small deletions, `start:end:abundance`
`annotation` | CDS | if an annotation file was provided this column contains the feature in which this base falls



## Normalization methods

By default no normalization is performed and reported counts are raw abundances, which are not suitable for comparing TEs between samples. Therefore two different strategies are implemented, normalization per million mapped reads and normalization by single-copy genes.

### Per million mapped reads

This option normalizes all counts (coverage, polymorphisms and structural variants) per million mapped reads and thus accounts for different depth of sequencing when comparing two or more samples. This option can be selected in the wrapper script `deviaTE` or in the second step of the workflow `deviaTE_analyse` by simply adding the argument `--rpm`

### Single-copy gene normalization

This normalization method relates all counts to the number of reads mapping to one or more genes, which are present only once in the genome of the investigated species (i. e. single-copy genes). Therefore it accounts for the sequencing depth between samples and additionally obtains an estimate of the insertion/copy number of the analysed TE. For this normalization method you can add the sequence of multiple single-copy genes to the TE consensus sequences used as `--library` and then activate the normalization by adding `--single_copy_genes GENE1,GENE2,GENE3...` to either `deviaTE` or `deviaTE_analyse`, where GENE1 etc. are the headers with which these genes appear in the library file. The estimated copy number per haploid genome can then be found in the header-section of the resulting output table.



### Special use-case: *Drosophila*

If you are analyzing TEs in *Drosophila* you do not have to specify a `--library` or `--annotation` of reference sequences, since we provide consensus sequences with our tool. When choosing which `--families` to analyze any elements from the first column (ID) of the [available TE consensus sequences in *Drosophila*](https://github.com/W-L/deviaTE/blob/master/deviaTE/lib/te_table) can be selected.

For single-copy gene normalization we have also added five genes (Dmel_rpl32, Dmel_piwi, Dmel_act5C, Dmel_ras and Dmel_p53), which can be used in the following way: `--single_copy_genes Dmel_rpl32,Dmel_piwi...`


### Special use-case: Plot multiple TE families from one or more samples

DeviaTE can handle plotting of an arbitrary number of TE families in one or more samples and automatically aligns plots by TE (column) and samples (row). To produce such a grid of plots, you can simply concatenate multiple output tables from `deviaTE_analyse` using the standard Unix tool `cat`. The sequence of your files for concatenation is irrelevant.

```
cat sample1.TE1 sample2.TE2 sample1.TE2 sample2.TE1 > allSamples_allTEs
deviaTE --input allSamples_allTEs
```

Depending on how many plots your figure consists of you may wish to increase the default text size in plots from 14 using `--fontsize`.


### Special use-case: Already mapped reads

You might have already mapped the sequening reads of your experiments and want to analyze TEs in these samples without having to remap your fastq files. In this case you can either use the wrapper script `deviaTE` and just substitute the input command to deviaTE with `--input_bam` for a single alignment file or with `--input_bam_dir` for a directory containing multiple files.
Alternatively, you can also run `deviaTE_analyse` with your alignment file, basically skipping the first step of the workflow and then plot outputs individually. In any case you still need to provide the reference sequences used for mapping as `--library` to the script.





