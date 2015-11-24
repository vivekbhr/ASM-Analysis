#!/bin/bash

### All the required variable will be stored first
# usage information
usage() { echo -e "Usage: $0 -w <working-dir> -c <config-file> -m <maternal-strain-prefix> -p <paternal-strain-prefix> -b <blacklist-regions> -n <sample-name> -t <num-processors>"; exit 1; }

# parse commandline arguments
while getopts ":w:c:m:p:b:n:t:" arg; do
    case $arg in
        w) workdir=$OPTARG
           ;;
        c) config=$OPTARG
           ;;
        m) MAT_STRAIN=$OPTARG #129S1 # refers to vcf column
           ;;
        p) PAT_STRAIN=$OPTARG  #CASTEiJ # refers to vcf column
           ;;
        b) blklist=$OPTARG  #blacklist regions to remove "/data/manke/repository/misc/annotations/blacklist_ENCODE/mm9-blacklist.bed"
           ;;
        n) sample=$OPTARG
           ;;
        t) proc=$OPTARG  # num of processors
           ;;
        \?) usage
           ;;
        :) usage
           ;;
    esac
done
shift $((OPTIND-1))

## Check if arguments are valid
if [[ -z "${workdir}" ]] || [[ -z "${MAT_STRAIN}" ]] || [[ -z "${PAT_STRAIN}" ]] || 
   [[ -z "${sample}" ]] || [[ -z "${proc}" ]] || [[ -z "${config}" ]] ; then
    usage
fi

if [[ ! -d ${workdir} ]]; then
    echo "ERROR: Input directory ${workdir} does not exist."; exit 1
elif [[ ! -d "${workdir}/01_rawdata/fastq" ]]; then
    echo "ERROR: Input fastq directory ${workdir}/01_rawdata/fastq does not exist. Please create it and place fastqs with sample names in there."; exit 1
elif [[ ! -d "${workdir}/01_rawdata/pseudogenome" ]]; then
    echo "ERROR: Input pseudogenome directory ${workdir}/01_rawdata/pseudogenome does not exist. Please create it and place fasta file, index files and mod files with maternal/paternal prefix names in it."; exit 1
elif [[ ! "${config}" ]]; then
    echo "ERROR: Config file not provided."; exit 1
fi

## Source config file
source ${config}

# Raw files
	fastq=${workdir}/01_rawdata/fastq
	pseudogen=${workdir}/01_rawdata/pseudogenome # this dir should have mat and pat strain fasta, bwt/hisat indexes, and mod files, the names should be like : <genotype>.fa, <genotype>.bt2, <genotype>.mod etc..

# Directories to make
bowtieOut=${workdir}/02_mappingToPseudogenome
refmapdir=${workdir}/03_mappingBackToRefgenome
mergedBAMs=${workdir}/04_mergeMappings
filteredBAMs=${workdir}/05_filterMergedBAMs
split=${workdir}/06_splitFilteredBAMs
#out=${workdir}/07_peakCalling_unSplitBAMs

for dir in ${bowtieOut} ${refmapdir} ${mergedBAMs} ${filteredBAMs} ${split} ${out}; do
	if [[ ! -d ${dir} ]]; then
		mkdir ${dir}; else 
		echo "directory ${dir} exists. Files with same samplename will be replaced"
	fi ; done

## ------------------------------------------------------ Map and convert
for genotype in ${MAT_STRAIN} ${PAT_STRAIN}
do 
## 01 Map
	echo "Sample : " ${sample} ". Mapping to pseudogenome : " $pseudogen/${genotype}
	${bwt} -x $pseudogen/${genotype} \
	-1 ${fastq}/${sample}_R1.fastq.gz \
	-2 ${fastq}/${sample}_R2.fastq.gz \
	-X 1000 -p ${proc} --rg-id mpi-ie --rg CN:deep_sequencing_unit --rg PL:illumina \
	| /package/samtools/samtools view -Sb - | /package/samtools/samtools sort -@ ${proc} - \
	$bowtieOut/${genotype}_${sample}
	# Index
	${samtools} index ${bowtieOut}/${genotype}_${sample}.bam
# Copy the mod file
	cp ${pseudogen}/${genotype}.mod ${bowtieOut}/${genotype}.mod

## 02 Map to Ref (Lapels)
	echo "Sample : " ${sample} ". Mapping back to reference genome."
	${lapels} -f -p ${proc} -o ${refmapdir}/${genotype}_${sample}_mapToRef.bam \
	${bowtieOut}/${genotype}.mod ${bowtieOut}/${genotype}_${sample}.bam \
	> ${genotype}_${sample}_stdout.txt 2> ${genotype}_${sample}_stderr.txt

## 03 sort 
	${samtools} sort -@ ${proc} -T ${genotype}_${sample} -O bam -n -o ${refmapdir}/${genotype}_${sample}_mapToRef.Rdsortd.bam \
	${refmapdir}/${genotype}_${sample}_mapToRef.bam

# Remove copied mod file
	rm ${bowtieOut}/${genotype}.mod*
done


## ------------------------------------------------------- Merge and Filter
# reset pythonpath
export PYTHONPATH=
## 04 merge
echo "Merging mapping : "${sample} ". Files :" ${MAT_STRAIN}_${sample}_mapToRef.Rdsortd.bam ${PAT_STRAIN}_${sample}_mapToRef.Rdsortd.bam

${suspenders} --pileup -p ${proc} -c ${mergedBAMs}/pileupImage_${sample}.png \
	${refmapdir}/${sample}_suspMerged.bam \
	${refmapdir}/${MAT_STRAIN}_${sample}_mapToRef.Rdsortd.bam \
	${refmapdir}/${PAT_STRAIN}_${sample}_mapToRef.Rdsortd.bam

mv ${refmapdir}/${sample}_suspMerged.bam ${mergedBAMs}/


## Filtering for random and blacklisted regions in the dir (need to make blklist optional)
echo "Filtering and sorting. Sample : ${sample} " 
${samtools} sort -@ ${proc} -T ${sample} ${mergedBAMs}/${sample}_suspMerged.bam -O sam | ${samfilt} --filter_out_from_BED ${blklist} --random \
--chrM --lowqual | ${samtools} sort -@ ${proc} -T ${sample} -O bam -o ${filteredBAMs}/${sample}_filt.bam - # needs sorting here
# index
${samtools} index ${filteredBAMs}/${sample}_filt.bam

## AllelicFilter.py : Filtering merged BAMS
echo "Splitting by alleles. Sample : ${sample} "
${allelefilt} --BAMfile ${filteredBAMs}/${sample}_filt.bam --outfile1 ${split}/${sample}_129S1_sep.bam --outfile2 ${split}/${sample}_CASTEiJ_sep.bam --outfile3 ${split}/${sample}_cannotTell.bam --removeMultiMapped --coordinateSorting 

## DONE

echo " Done..!! Now prepare samplesheet for AlleleEx-wrapper. "
