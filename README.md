# deviaTE

deviaTE is a python tool for the analysis and visualization of mobile genetic element sequences.

## Dependencies

* python 3.6+
* pip for python3 (https://pip.pypa.io/en/stable/installing/)
* samtools (tested with version 1.6 and 1.9)
* bwa (tested with version 0.7.17)
* R (tested with version 3.4.3)
  * required R packages are installed automatically upon first usage

samtools and bwa need to be in your `$PATH`


## Installation

```pip3 install git+https://github.com/W-L/deviaTE/```

## Usage manual and walkthroughs

A manual and a walkthrough with examples using publicly available data can be found here:

* [Manual](https://github.com/W-L/deviaTE/blob/master/doc/MANUAL.md) 
* [Walkthrough](https://github.com/W-L/deviaTE/blob/master/doc/WALKTHROUGH.md) 


## Quick start guide

To produce a visualization of transposable element families in *Drosophila* from sequencing reads, you can use

```deviaTE --input_fq foo.fq --families TEfamily1,TEfamily2,...```

where TEfamily1 etc. are elements from the first column (ID) of the [available TE consensus sequences](https://github.com/W-L/deviaTE/blob/master/deviaTE/lib/te_table), separated by commas

DeviaTE can also be applied to multiple files in a folder using:

```deviaTE --input_fq_dir --families TEfamily1,TEfamily2,...```

Other available arguments can be seen with ```deviaTE -h``` and are documented in the [Manual](https://github.com/W-L/deviaTE/blob/master/doc/MANUAL.md) 


## Bug reports
If you find any problems please use the issue tracker on this repository


## License
deviaTE is licensed under the GPLv3 License, see LICENSE.

