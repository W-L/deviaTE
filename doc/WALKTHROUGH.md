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

Additionally an output table is produced in both cases, which has the TE family name as file ending. 
This is a space-separated file with a header denoted by #. The first line contains a timestamp, the command used to generate the file,
the estimated number of TE insertions on the following line (only when single gene normalization was selected, see below)
and finally the actual header of the table. 

## Description of output

The following table is a list of the columns in the output of deviaTE

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


