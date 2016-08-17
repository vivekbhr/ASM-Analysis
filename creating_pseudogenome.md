## Creating pseudogenome to use with AlleleExplorer


* MOD file -- if the MOD file is not available, it must be created or downloaded
* For mouse genomes, Check whether your MOD file is available [here](http://www.csbio.unc.edu/CCstatus/index.py?run=Pseudo "MOD files supplied by Huang et al.") 
* Otherwise, you will need VCF files with the information for your strain(s) of interest and their differences to the reference genome and the reference genome FASTA file:

On terminal, type :

```
REF_FASTA=/path/to/ref.fa # fasta file for reference genome, e.g. mm9, dm3, ce10....
VCF_SNP=/path/to/snps.vcf # must contain the information about genetic variants for at least both strains of interest
VCF_INDELS=/path/to/indels.vcf # like snps.vcf, contains information about the genetic variations
```

##### information

These variables should be set regardless whether a MOD file is available already or not

On terminal, type :

```
REF_GENOME=ce10 # reference genome
MAT_STRAIN=A # maternal strain # must match entry from to vcf column
PAT_STRAIN=B # paternal strain # must match entry from vcf column
CHROMOSOMES=1,2,3,X # in case, one is interested in a couple of chromosomes only
```

##### programs

###### MOD file generation

lapels comes with a set of supplementary scripts that can be used to generate the MOD file if needed:

1. get_refmeta
2. vcf2mod
3. insilico
