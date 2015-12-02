
# AlleleExplorer
A pipeline for allele-specific mapping and quantification of RNA-Seq and ChIP-Seq data.

## Contents

### AlleleMap-wrapper
A pipeline for allele-specific mapping of reads.

The pipeline uses bowtie2,tophat2 or hisat2 to map, and lapels and suspenders tools to convert
genomic coordinates between reference genome and maternal/paternal genome.
It also performs various filtering steps.

### AlleleEx-wrapper
Based on AlleleExplorer R package for quantification of allele-specific differential expression and binding.

AlleleExplorer can be installed and used in R. After installation it can also be run from command line using the wrapper script.


## installation

1) First Install the mapper :
    + [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
    + [tophat2](https://ccb.jhu.edu/software/tophat/index.shtml)
    + [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)
2) Install [suspenders](https://github.com/holtjma/suspenders).
3) Install [lapels](https://pypi.python.org/pypi/lapels).
4) Install [modtools](https://pypi.python.org/pypi/modtools/1.0.2).
5) Clone this github repository and install the AlleleExplorer R package in the terminal.

## Usage



========================
***Author: @vivekbhr***

***License: GPLv3***
