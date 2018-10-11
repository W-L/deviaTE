# Usage Manual for deviaTE

DeviaTE is a command-line program that analyzes and visualizes the diversity mobile genetic elements from sequencing data without the need for an assembled genome of the host species. Besides coverage, it displays Polymorphisms, indels, truncations and internal deletions and allows for easy comparative analyses. 

This tool is executed by running `deviaTE_main`, which is a convenient wrapper for the three tiers, `deviaTE_prep`, `deviaTE_analyse` and `deviaTE_plot` to produce a visualization of a transposon (see Fig.). These scripts are also included and can be called individually. Available arguments for all tools can be viewed by adding `-h`



![Architecture](https://github.com/W-L/deviaTE/blob/master/doc/architecture.png "Architecture of deviaTE")



## `deviaTE_main` 

### Description

This script executes the full pipeline to produce a table describing each position of the TE and a plot of that table.

### Arguments

#### Required

Argument name | Type | default | Description
`--input_fq INPUT_FQ` | String | 
`--input_bam INPUT_BAM` | String |
`--input_fq_dir` |  String |
`--input_bam_dir` | String |

[-h] [--library LIBRARY] [--read_type READ_TYPE]
                    [--min_read_len MIN_READ_LEN]
                    [--quality_threshold QUALITY_THRESHOLD]
                    [--min_alignment_len MIN_ALIGNMENT_LEN]
                    [--threads THREADS] --families FAMILIES
                    [--annotation ANNOTATION] [--no_freq_corr]
                    [--hq_threshold HQ_THRESHOLD] [--free_yaxis]
                   