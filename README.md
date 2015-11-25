
# AlleleExplorer
A pipeline for allele-specific mapping and quantification of RNA-Seq and ChIP-Seq data.

## Contents

### AlleleMap-wrapper
A pipeline for allele-specific mapping of reads.

The pipeline uses bowtie2 to map, and lapels and suspenders tools to convert genomics coordinates between reference
genome and maternal/paternal genome.
It also performs various filtering steps.

### AlleleEx-wrapper
Based on AlleleExplorer R package for quantification of allele-specific differential expression and binding.

AlleleExplorer can be installed and used in R. After installation it can also be run from command line using the wrapper script.


***Author: @vivekbhr***

***License: GPLv3***
