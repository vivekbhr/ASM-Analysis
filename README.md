
# AlleleExplorer
A pipeline for allele-specific mapping and quantification of RNA-Seq and ChIP-Seq data.

## Contents

### AlleleMap-wrapper
A pipeline for allele-specific mapping of reads.

The pipeline uses bowtie2 or tophat2 to map, and lapels and suspenders tools to convert
genomic coordinates between reference genome and maternal/paternal genome.
It also performs various filtering steps.

### AlleleEx-wrapper
Based on AlleleExplorer R package for quantification of allele-specific differential expression and binding.

AlleleExplorer can be installed and used in R. After installation it can also be run from command line using the wrapper script.


## Installation

1) First Install one of the mappers :
   * For Chip-Seq : [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
   * For RNA-Seq : [tophat2](https://ccb.jhu.edu/software/tophat/index.shtml).

2) Install [suspenders](https://github.com/vivekbhr/suspenders). : This is a modified version of original suspenders package by @holtzma, with multiple bigfixes.

3) Install [lapels](https://pypi.python.org/pypi/lapels).

4) Install [modtools](https://pypi.python.org/pypi/modtools/1.0.2).

5) Clone this github repository and install the AlleleExplorer R package in the terminal.

## Usage

### Preparing Files

1) To use AlleleExplorer, you need to first create the maternal and paternal pseudogenomes using ModTools. Check out how to do it [here](https://github.com/vivekbhr/AlleleExplorer/blob/master/creating_pseudogenome.md)

2) Place the pseudogenome and the bowtie/tophat indexes in a directory named *01_rawdata/pseudogenome* 
   in your working directory.

3) Place the raw fastq files in the directory  *01_rawdata/fastq* in your working directory.

### Run AlleleExplorer

3) Run AlleleMap-wrapper (see tool help), to map the fastqs to reference genome. You will need a **config file** for this, explaining the executable locations (eg. [alleleMap.config](./alleleMap.config) ).

5) Prepare a sampleSheet for AlleleEx-wrapper, indicating the maternal/paternal mapped file location, sample names etc. (see example samplesheet).

6) Run AlleleEx-wrapper as `Rscript AlleleEx-wrapper` (see help, -h) for help.



========================
***Author: @vivekbhr***

***License: GPLv3***
