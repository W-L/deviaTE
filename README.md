[![Anaconda-Server Badge version](https://anaconda.org/w-l/deviate/badges/version.svg)](https://anaconda.org/w-l/deviate)
[![Anaconda-Server Badge platforms](https://anaconda.org/w-l/deviate/badges/platforms.svg)](https://anaconda.org/w-l/deviate)
[![Anaconda-Server Badge lastupdate](https://anaconda.org/w-l/deviate/badges/latest_release_date.svg)](https://anaconda.org/w-l/deviate)
[![Anaconda-Server Badge downloads](https://anaconda.org/w-l/deviate/badges/downloads.svg)](https://anaconda.org/w-l/deviate)

# deviaTE

deviaTE is a python tool for the analysis and visualization of mobile genetic element sequences. It is available for Unix and Linux systems.

## Installation

### Using Conda

The recommended way of using deviaTE is in a conda virtual environment. This way, all dependencies should be resolved.

Create a new conda virtual environment, and install deviaTE:

```conda create deviaTE -c r -c defaults -c conda-forge -c bioconda -c w-l -n deviaTE_env```

this command loads all required anaconda channels (`r, defaults, conda-forge, bioconda, w-l`) and names the environment `deviaTE_env`

The environment can then be activated, and deviaTE with all necessary dependencies can be used inside the environment:

```source activate deviaTE_env```

After using deviaTE, the environment can be exited using:

```source deactivate```


### Using pip

DeviaTE can also be installed with the pip python package manager. However, this does not take care of all dependencies. The following tools need to be set up manually beforehand:

* python 3.6+
* pip for python3 (https://pip.pypa.io/en/stable/installing/)
* samtools (tested with version 1.6 and 1.9)
* bwa (tested with version 0.7.17)
* R (tested with version 3.4.3)
  * required R packages are installed automatically upon first usage

samtools and bwa need to be in your `$PATH`

DeviaTE can then be installed with 

```pip3 install deviaTE```


## Uninstallation

conda environment:

```conda env remove -n deviaTE_env```

pip package:

```pip3 uninstall deviaTE```


## Usage manual and walkthroughs

A manual and a walkthrough with examples using publicly available data can be found here:

* [Manual](https://github.com/W-L/deviaTE/blob/master/doc/MANUAL.md) 
* [Walkthrough](https://github.com/W-L/deviaTE/blob/master/doc/WALKTHROUGH.md) 


## Bugs, issues, questions, suggestions ...
If you find any problems, have questions or ideas for further improvement please use the issue tracker on this repository, thanks!


## License
deviaTE is licensed under the GPLv3 License

