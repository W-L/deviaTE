# Usage Manual for deviaTE

DeviaTE is a command-line program that analyzes and visualizes the diversity of mobile genetic elements from sequencing data without the need for an assembled genome of the host species. Besides coverage, it displays single nucleotide polymorphisms, indels, truncations and internal deletions and allows for easy comparative analyses. 

The tool operates in three individual steps (see Fig.). It first prepares sequencing reads (`deviaTE_prep`), analyses a transposon chosen from the consensus sequence library (`deviaTE_analyse`) and finally produces a visualization from the output of the previous step (`deviaTE_plot`). Alternatively `deviaTE` is a convenient, single-command wrapper. Available arguments for all scripts are listed by adding `-h` or `--help`.

<p align="center">
  <img src="https://github.com/W-L/deviaTE/blob/master/doc/workflow.png" alt="Architecture of deviaTE" width="600"/>
</p>

## `deviaTE` 

### Description

This script executes the full pipeline to produce a table describing each position of the TE and a visualization of the TE using that table. The input can either be sequencing reads (FASTQ), already mapped reads (BAM) or directories of either of those.

### Arguments

#### Required

Argument name | Type | default | Description
--- | --- | --- | ---
`--input_fq` | string | none | Sequencing reads to be analysed. Part of mutually exclusive argument group `--input_*`
`--input_bam` | string | none | Alternatively already mapped reads can be analysed. Part of mutually exclusive argument group `--input_*`
`--input_fq_dir` |  string | none | Use all fastq files in a directory instead of an individual file. Part of mutually exclusive argument group `--input_*`
`--input_bam_dir` | string | none | Use all bam files in a directory instead of an individual file. Part of mutually exclusive argument group `--input_*`
`--families` | string | none | Comma separated list of TE families to be analysed, no spaces. E.g. PPI251,DMLINEJA,DMCOPIA

#### Optional

Argument name | Type | default | Description
--- | --- | --- | ---
`--library` | string | included consensus set for *D. melanogaster* | FASTA file of reference consensus sequences of TE families. Items used in the list for `--families` must be headers in this file
`--read_type` | string | sanger | type of sequencing read to determine quality encoding for trimming (either `illumina` or `sanger`)
`--rpm` | bool | false | normalize abundances per million mapped reads. Mutually exclusive with single copy genes
`--single_copy_genes` | string | false | Comma separated list of single copy genes to normalize abundances with. Must be present in file defined by `--library`. Mutually exclusive with rpm
`--annotation`| string | included annotations for *D. mel.* consensus set | GFF3 file with annotations of the TE consensus sequences
`--min_read_len` | int | 1 | minimum length of reads to be retained
`--min_alignment_len` | int | 1 | minimum length of aligned segments to be considered
`--quality_threshold` | int | 15 | minimum base quality for filtering, uses modified Mott algorithm
`--hq_threshold`| int | 20 | threshold of which mapped segments to consider ambiguous or unambiguous
`--no_freq_corr` | bool | false | deactivates the correction of estimated frequencies of internal deletions
`--free_yaxis`| bool | false | frees the y axis during plotting, e.g. to make very different y axis ranges visible
`--threads`| int | 1 | use this many threads for whichever subroutines possible



## `deviaTE_prep`

### Description

First step of the pipeline. Trims, maps and filters reads. Also performs the detection of internal deletions within the mapped reads. Input is a file with sequencing reads in FASTQ format. Outputs an alignment file ready to be processed in the subsequent steps.

### Arguments

#### Required

Argument name | Type | default | Description
--- | --- | --- | ---
`--input` | string | none | Sequencing reads to be analysed. In FASTQ format

### Optional 

Argument name | Type | default | Description
--- | --- | --- | ---
`--library` | string | included consensus set for *D. melanogaster* | FASTA file of reference consensus sequences of TE families
`--quality_encoding` | string | sanger | type of sequencing read to determine quality encoding for trimming (either `illumina` or `sanger`)
`--qual_threshold` | int | 15 | minimum base quality for filtering, uses modified Mott algorithm
`--min_read_length` | int | 1 | minimum length of reads to be retained
`--min_alignment_length` | int | 1 | minimum length of aligned segments to be considered
`--threads`| int | 1 | use this many threads for whichever subroutines possible      
`--nofuse`| bool | false | skip the detection of internal deletions by novel algorithm



## `deviaTE_analyse`

### Description 

Second and main step of the deviaTE workflow. Uses the previously prepared bam file and produces quantitative information and estimates at nucleotide resolution for the chosen TE family.

### Arguments 

#### Required

Argument name | Type | default | Description
--- | --- | --- | ---
`--input` | string | none | Prepared bam file to be analysed
`--family` | string | none | TE family for which information should be collected. Must be one of the headers of reference sequences in `--library`

### Optional 

Argument name | Type | default | Description
--- | --- | --- | ---
`--library` | string | included consensus set for *D. melanogaster* | FASTA file of reference consensus sequences of TE families
`--annotation`| string | included annotations for *D. mel.* consensus set | GFF3 file with annotations of the TE consensus sequences
`--output` | string | input + TE family name | User specified name for the output table
`--sample_id` | string | input | name that is given to the sample in the output table and also displayed in the resulting visualization
`--no_freq_corr` | bool | false | deactivates the correction of estimated frequencies of internal deletions
`--hq_threshold`| int | 20 | threshold of which mapped segments to consider ambiguous or unambiguous
`--rpm` | bool | false | normalize abundances per million mapped reads. Mutually exclusive with single copy genes
`--single_copy_genes` | string | false | Comma separated list of single copy genes to normalize abundances with. Must be present in file defined by `--library`. Mutually exclusive with rpm


## `deviaTE_plot`

### Description

Final step of the pipeline. Produces a visualization of the TE family analysed by the previous script. Input is the table from `deviaTE_analyse`, output is a visualization in either PDF or EPS format. Multiple TE families can be visualized by simply combining their output tables, for instance: `cat sample1.COPIA sample2.COPIA sample3.COPIA >COPIA` and using the merged file as input to the plotting script.

### Arguments

#### Required 

Argument name | Type | default | Description
--- | --- | --- | ---
`--input` | string | none | Table of TE family information to be visualized

#### Optional 

Argument name | Type | default | Description
--- | --- | --- | ---
`--output` | string | input + .pdf/.eps | User defined name for the resulting illustration
`--out_format` | string | pdf | Format of the output. Either pdf or eps. 
`--free_yaxis`| bool | false | frees the y axis, e.g. to allow for different y axis ranges for different TE families
`--fontsize` | int | 14 | Fontsize of any text in the plot, e.g. title, axis text and axis ticks












                    
                   
