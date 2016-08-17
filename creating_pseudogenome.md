## Creating pseudogenome to use with AlleleExplorer

### You will need a MOD file

MOD file is a map of coordinates between genomes. If the MOD file is not available, it must be created or downloaded

* For mouse genomes, Check whether your MOD file is available [here](http://www.csbio.unc.edu/CCstatus/index.py?run=Pseudo "MOD files supplied by Huang et al.") 
* Otherwise, you will need VCF files with the information for your strain(s) of interest and their differences to the reference genome and the reference genome FASTA file:

### Set some Environment Variables.

On terminal, type :

```
REF_FASTA=/path/to/ref.fa # fasta file for reference genome, e.g. mm9, dm3, ce10....
VCF_SNP=/path/to/snps.vcf # must contain the information about genetic variants for at least both strains of interest
VCF_INDELS=/path/to/indels.vcf # like snps.vcf, contains information about the genetic variations

REF_GENOME=dm6 # name of reference genome
MAT_STRAIN=A # maternal strain (must match entry from to vcf column)
PAT_STRAIN=B # paternal strain (must match entry from vcf column)
CHROMOSOMES=1,2,3,X # in case, one is interested in a couple of chromosomes only
```

### Installing MODtools

I recommend installing MiniConda/AnaConda before installing MODtools.

#### Installing MiniConda (terminal)

```
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh # 64 bit bash installer
bash Miniconda3-latest-Linux-x86_64.sh
```

#### Install MODtools

After you install conda, restart your terminal and type :

```
pip install modtools
```

### Running MODtools

MODtools contain three programs.

1. get_refmeta
2. vcf2mod
3. insilico


###  1. Create Meta-data for reference genome

**-->  get_refmeta**

```
get_refmeta -o ${REF_GENOME}.meta ${REF_GENOME} ${REF_FASTA}
```


### 2. generating MOD files

##### a) indexing of VCF files

```
for vcf in ${VCF_INDELS} ${VCF_SNP}
do
    tabix-0.2.6/bgzip -c ${vcf} > ${vcf}.gz
	  tabix-0.2.6/tabix -p vcf ${vcf}.gz
done
```
##### b) VCF to MOD format

**--> vcf2mod**

1 MOD file per VCF file (SNPs, INDELS) and genotype (maternal, paternal)
SNPs with bad quality (FI tag = 0) will be discarded

```
vcf2mod -c ${CHROMOSOMES} -f \
  -o ${REF_GENOME}_SNPs_${genotype}.mod \
  ${REF_GENOME} ${REF_GENOME}.meta \
	${genotype} ${VCF_SNP} 2> vcf2mod.log
```

##### c) merge MOD files for SNPs and indels per genotype

```
cat ${REF_GENOME}_SNPs_${genotype}.mod ${REF_GENOME}_indels_${genotype}.mod |\
sort -k2,2n -k3,3n | uniq > ${REF_GENOME}_indels_SNPs_${genotype}_changedChr.mod 
```

### 3. Generating pseudogenomes

**--> insilico**

CAVE: after this step, the MOD file will be gzipped (without any indication in the file name)

```
insilico \
	${REF_GENOME}_indels_SNPs_${genotype}_changedChr.mod \
	${REF_FASTA} -v -f \
	-o pseudogenome_${REF_GENOME}_${genotype}.fa > insilico.log 
```

---------------------------------------
