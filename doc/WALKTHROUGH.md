# Walkthrough for deviaTE

## Getting started

Right after installing and skipping through the manual we recommend to test that the tool works properly. 
We have prepared a FASTQ file containing only reads of *Drosophila melanogaster* that map to the transposable element *jockey* 
(can be downloaded [here](https://github.com/W-L/deviaTE/blob/master/example/jockey_dmel.fastq) [~3Mb, in the example folder of this repository]).
These sequences are sanger reads from Drosophila 12 Genomes Consortium et al. 2007. Evolution of genes and genomes on the Drosophila phylogeny. *Nature*. 450(7167):203-218.

You can now use the convenient wrapper script `deviaTE` to analyse the TE jockey (DMLINEJA) and get a visualization:

```bash
deviaTE --input_fq jockey_dmel.fastq --read_type phred+33 --families DMLINEJA
```

Alternatively you can run all three steps of deviaTE manually with: 

```bash
deviaTE_prep --input jockey_dmel.fastq --quality_encoding phred+33
deviaTE_analyse --input jockey_dmel.fastq.fused.sort.bam --family DMLINEJA
deviaTE_plot --input jockey_dmel.fastq.fused.sort.bam.DMLINEJA
```

this will first trim, map and filter the reads, which yields an alignment file called `jockey_dmel.fastq.fused.sort.bam`
then the tool collects the quantitative information and calculates estimators, which produces the output table `jockey_dmel.fastq.fused.sort.bam.DMLINEJA`. For a detailed description of this output, see below.
Finally a visualization is produced, resulting in an illustration of jockey that should look like [this](https://github.com/W-L/deviaTE/blob/master/example/jockey_dmel.fastq.DMLINEJA.pdf)


**If deviaTE does not work at this point, please consult the Troubleshooting section at the bottom of this document or do not hesitate to open an issue on this GitHub**


## Quality encoding of the FASTQ file

Depending on the sequencing platform used to produce your sequences, quality scores are encoded with certain score offset.  The following parameters should be used for `--read_type`:

* Sanger and Illumina version 1.8+ reads are encoded in: `phred+33`
* Illumina reads from version 1.3 to pre 1.8: `phred+64`


This is a rough and simple way to find out about the encoding: If the quality scores contain character 0 it is phred+33. If the quality scores do not contain 0, it is phred+64. More information can be found on the [FASTQ-wikipage](https://en.wikipedia.org/wiki/FASTQ_format#Encoding)



## Basic usage

<p align="center">
  <img src="https://github.com/W-L/deviaTE/blob/master/doc/workflow.png" alt="Architecture of deviaTE" width="600"/>
</p>



### Using the single-command wrapper script

To analyze and visualize transposable element families from sequencing reads the most basic usage is:

```
deviaTE --input_fq foo.fastq --families TEfamily1,TEfamily2,... --library TE_consensus_sequences.fasta
```

where foo.fastq contains sequencing reads and TEfamily1 etc. are headers of TE reference sequences in the file defined by `--library`. 

Any fasta file used as library must first be indexed with:

```
bwa index TE_consensus_sequences.fasta
```

Multiple families can be selected, they simply need to be separated by commas without spaces. DeviaTE can also be applied to multiple fastq files in a directory using `--input_fq_dir` and running the program from within that directory. 

Other available arguments can be seen with `deviaTE -h/--help` and are documented in the [Manual](https://github.com/W-L/deviaTE/blob/master/doc/MANUAL.md) 

The result will be a table of quantitative information and a visualization for each selected TE families in each of the specified samples.



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

Where the input is the alignment file from the previous step, the library contains the reference sequences and the selected family is the header of a sequence within that library. Additional arguments can be used to provide an annotation of the TE sequences, name the sample and the output, set a quality threshold for unambiguously mapped reads or to select a normalization method.

This script produces the output table containing the quantitative information, which by default will be named as: `input.TEfamily`. This file is then used to produce a visualization, e. g. in basic scenarios:

```
deviaTE_plot --input foo.fastq.fused.sort.bam.TEfamily
```

Here, optional arguments can be used to change the output name and its format (pdf or eps), change fontsize or allow the y-axis range to be set individually for multiple plots in a grid, e. g. when more than one sample is plotted and the coverage between them is highly different. For a description on how to plot multiple samples and TE families at once, see below.


## Description of the output table

Besides the visualization an output table is produced, which contains the quantitative information about the TE family. 
The table starts with some header-lines denoted by #. The first line contains a timestamp and the command used to generate the file, while the following line states the estimated number of TE insertions (only when single gene normalization was selected, see below). The final line of the header section contains the column names of the table. And each subsequent row in the table corresponds to one position of the TE sequence.

The following is a list describing the columns in the output table of deviaTE:

column name | example value | description
--- | --- | ---
`TEfam` | DMLINEJA | Name of the analysed TE family
`sample_id` | jockey_dmel.fastq | either the input file name or an identifier set by `--sample_id`
`pos` | 709 | position in the reference sequence
`refbase` | T | Nucleotide in the reference sequence at this position
`A C G T` | 124 | counts of each nucleotide at this position
`cov` | 140 | total coverage at this position
`phys_cov` | 18 | physical coverage at this position, e.g. spanned by deletion
`hq_cov` | 140 | coverage above the specified threshold for unambigously mapped reads (`--hq_threshold`)
`snp` | True | denotes a single nucleotide polymorphism at this position
`refsnp` | False | denotes a fixed difference to the reference sequence
`int_del` | 709:771:25 | describes an internal deletion at this position in the form of `start:end:abundance`
`int_del_freq` | 709:771:0.4 | frequency of internal deletions in the format `start:end:frequency` (e.g. 0.4 means 40% of reads have this deletion)
`trunc_left` | 17 | number of reads that show a truncation to the left of this position
`trunc_right` | 17 | number of reads that show a truncation to the right of this position
`ins` | 709:711:2 | describes insertions at this position with the format `start:end:abundance`
`delet` | 709:710:3 | describes short deletions, `start:end:abundance`
`annotation` | CDS | if an annotation file was provided this column contains the feature in which this position lies



## Normalization methods

By default no normalization is performed and reported counts are raw abundances, which are not suitable for comparing TEs between samples. Therefore two different strategies are implemented, normalization per million mapped reads and normalization by single-copy genes.

### Per million mapped reads

This option normalizes all counts (coverage, polymorphisms and structural variants) per million mapped reads and thus accounts for different depth of sequencing when comparing two or more samples. This option can be selected in the wrapper script `deviaTE` or in the second step of the workflow `deviaTE_analyse` by simply adding the switch `--rpm`.

### Single-copy gene normalization

This normalization method relates all counts to the number of reads mapping to one or more genes, which are present only once in the genome of the investigated species (i. e. single-copy genes). Therefore it accounts for the sequencing depth between samples and additionally obtains an estimate of the insertion/copy number of the analysed TE. For this normalization method you can add the sequence of multiple single-copy genes to the file containing the TE consensus sequences used as `--library`. Then you can activate the normalization by adding `--single_copy_genes GENE1,GENE2,GENE3...` to either `deviaTE` or `deviaTE_analyse`, where GENE1 etc. are the headers with which these genes appear in the library file. If more than one gene is specified, the average of their coverage is used for normalization. The estimated copy number per haploid genome can then be found in the header-section of the resulting output table.



### Special use-case: *Drosophila*

If you are analyzing TEs in *Drosophila* specifying a `--library` or `--annotation` of reference sequences is optional, since we provide consensus sequences with our tool (from https://github.com/bergmanlab/transposons). When choosing which `--families` to analyze any elements from the first column (ID) of the [available TE consensus sequences in *Drosophila*](https://github.com/W-L/deviaTE/blob/master/deviaTE/lib/te_table) can be selected.

For single-copy gene normalization we have also added five genes (Dmel_rpl32, Dmel_piwi, Dmel_act5C, Dmel_RpII140 and Dmel_p53), which can be used in the following way: `--single_copy_genes Dmel_rpl32,Dmel_piwi...`


### Special use-case: Plot multiple TE families from one or more samples output tables

DeviaTE can handle plotting of an arbitrary number of TE families in one or more samples and automatically aligns plots by TE (column) and sample (row). To produce such a grid of plots, you can simply concatenate multiple output tables from `deviaTE_analyse` using the standard Unix tool `cat`, and then feed these merged output tables to the plotting script. The order of your concatenated files is irrelevant.

```
cat sample1.TE1 sample2.TE2 sample1.TE2 sample2.TE1 > TwoSamples_TwoTEs
deviaTE_plot --input TwoSamples_TwoTEs
```

Depending on how many plots your figure consists of you may wish to increase the default text size in plots using `--fontsize`.


### Special use-case: Already mapped reads

You might have already mapped sequencing reads of your sample to a library of TEs and want to analyze another TE without remapping your fastq files. In this case you can either use the wrapper script `deviaTE` and simply substitute the input command with `--input_bam` for a single alignment file or with `--input_bam_dir` for a directory containing multiple files.
Alternatively, you can also run `deviaTE_analyse` with your alignment file, basically skipping the first step of the workflow and then plot outputs individually. In any case you need to provide the reference sequences used for mapping as `--library` to the scripts.

It is, however, not possible to analyze TEs from an alignment file, which contains sequencing reads mapped to a reference other than a library of TE consensus sequences. This is meant to decrease the time needed to analyze TEs from sequencing reads that have already been mapped to TE reference sequences.

### Special use-case: Analyse all reference sequences in the library

By specifying `--families ALL` it is possible to run `deviaTE` on all consensus sequences in a fasta library.


### Sepcial use-case: Skip the detection of internal deletions

If you already have an alignment file, in which internal deletions have been detected from a previous run (`FILE.fused.sort.bam`), or you simply are not interested in internally deleted TEs, you can deactivate the detection with the option `--no_delet_detect`.
This can be useful in cases where reads have successfully been mapped and deletions have been detected and a quick reanalysis and visualization is desired.

### Special use-case: Paired-end reads 

You can use DeviaTE for paired-end reads by mapping them in single read mode. This can be done, for example, by using a single concatenated fastq file that contains both read pairs (read1 and read2). In order to prevent some issues the reads then need to be given unique names (e.g. readID_1, readID_2, ..., readID_n), which can be achieved using a little script found at: example/rename_reads.py

```
python rename_reads.py sample.fq >sample_uniq.fq
```

Thanks to Anna for pointing out this issue!


## Troubleshooting section

### bash: deviaTE: command not found... OR EnvironmentError: (Errno 13) Permission denied

```
bash: deviaTE: command not found...
Could not install packages due to an EnvironmentError: [Errno 13] Permission denied:
```

If you encounter one of these errors during or after an installation with pip, the automatic transfer of deviaTE's executable scripts has probably failed. 
This can, for example, happen on a shared computing environment with restricted permissions.

To properly install deviaTE in such cases, you can use the `--user` flag with pip after declaring the installation directory (`~/deviaTE` here):

```
export PYTHONUSERBASE=~/deviaTE
pip3 install --user deviaTE
```


Afterwards, please run and add the following line to your `~/.bashrc` or `~/.bashprofile` to make deviaTE executable from anywhere in your system

```
export PATH=~/deviaTE/bin:$PATH
```

or run these lines to automatically add the location to your PATH:

```
printf "\nexport PATH=~/deviaTE/bin:$PATH\n" >> ~/.bash_profile
source ~/.bash_profile
```



### Error during visualization: unused argument 'cmd='

```
Error in fread(cmd = paste("grep -v ^#", inp_arg), data.table = FALSE,  : 
  unused argument (cmd = paste("grep -v ^#", inp_arg))
Execution halted
```

If you see an error like this, you probably have an outdated version of the R package data.table. Please update to a newer version (>1.11.6).


### libcrypto missing libraries

On MacOS, samtools sometimes encounters an error similar to:

*samtools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file: No such file or directory*

This can usually be solved by making sure that openssl is installed, and then linking the missing libraries manually:

```
brew update
brew install openssl

ln -s /usr/local/opt/openssl/lib/libcrypto.1.0.0.dylib /usr/local/lib/
ln -s /usr/local/opt/openssl/lib/libssl.1.0.0.dylib /usr/local/lib/
```

### Conda not accessing executables

Usually conda adds its environment to the user's PATH during the initial installation of conda/anaconda. If there are issues with that the following lines can be added to the user's .bashrc/.bash_profile.


miniconda:
```
export PATH="/home/YOURUSER/miniconda3/bin:$PATH"
export LD_LIBRARY_PATH=/home/YOURUSER/miniconda3/lib:/homes/YOURUSER/miniconda3/lib64:$LD_LIBRARY_PATH
```

anaconda:
```
export PATH="/home/YOURUSER/anaconda3/bin:$PATH"
export LD_LIBRARY_PATH=/home/YOURUSER/anaconda3/lib:/homes/YOURUSER/anaconda3/lib64:$LD_LIBRARY_PATH
```



