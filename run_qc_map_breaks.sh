#!/bin/bash

# Name: Basic bash pipeline - convert to a workflow engine later
# Author: Graham Rose
# Date: 14-06-2018
# Version 0.1
# Description: Breakpoint mapping analysis
#----------------------------------------------------

########## Load modules
module load fastqc # defaults to v.X
module load trimmomatic/0.36
module load bwa/0.7.15
module load crest/1.0 # requires Bio:DB:Sam module, and an aligner eg.
module load blat/35
module load samtools/1.2 # dependencie for CREST? v.0.1.17 maybe needed
module load bioperl/1.6.0
module load picard/2.6.0
module load breakdancer/1.2

########## Static vars ##########
DBS="/home/malikian/bp_mapping/dbs"
outputDir="results"


######### Functions/Modules ##########

# FastQC. Temp opts: --nogroup
func_fastqc() {
	echo "Running FastQC..."
	echo "FastQC on single-node using 4 Threads"
        fastqc -k 8 $1 -t 4
	fastqc -k 8 $2 -t 4
}

# Trimmomatic QC. Not clip fasta needs adapter seq. used. No idea on library for MiSeq. Vars created for later functions.
func_trimmomatic() {
	echo "Running Trimmomatic..."
	echo "PE mode, 4 threads"
	
	trimmomatic PE -phred33 -threads 4 -trimlog $fastq_base_name".trimmomatic.log" $1 $2 $fastq_r1_base"_PE.fastq" $fastq_r1_base"_SE.fastq" $fastq_r2_base"_PE.fastq" $fastq_r2_base"_SE.fastq" ILLUMINACLIP:$DBS/clip.fa:2:30:10 LEADING:20 TRAILING:20 MINLEN:40 	

	# New vars for trimmed PE fastq files
	fastq_r1_trimmed=$fastq_r1_base"_PE.fastq"
	fastq_r2_trimmed=$fastq_r2_base"_PE.fastq"

}

# Mapping and parse BAM. Genome fasta and bwa index in DBS variable declared at top
func_bwa() {
	echo "Running bwa mem..."
        echo "PE, 8 threads"
	echo "Ouput SAM filename $fastq_base_name"

	# Create syminks to genome index
	echo "Creating symlinks to indexed hg38 genome at /home/malikian/bp_mapping/dbs"
	ln -s /home/malikian/bp_mapping/dbs/hg38.fa.* .

	# Run bwa mem
	bwa mem -t 8 hg38.fa $1 $2 > $fastq_base_name".sam"

	# SAM>BAM
	echo "Convert to BAM..."
	picard SortSam VALIDATION_STRINGENCY=LENIENT INPUT=$fastq_base_name".sam" OUTPUT=$fastq_base_name".bam" SORT_ORDER=coordinate
	samtools index $fastq_base_name".bam"
}

func_mapping_stats() {

	# Basic mapping QC
        echo "Mapping QC stats..."
        samtools flagstat $fastq_base_name".bam" > $fastq_base_name".bam.flagstat.txt"
        samtools idxstats $fastq_base_name".bam" > $fastq_base_name".bam.idxstats.txt"
        samtools stats $fastq_base_name".bam" > $fastq_base_name".bam.stats.txt"
}


# CREST. As extractSclip.pl
func_crest() {
	#
	echo "Running CREST..."
	
	echo "Can't locate Bio/DB/Sam.pm"	

}

# BREAKDANCER. As bam2cfg.pl>breakdancer-max
func_breakdancer() {
	echo "Detecting inter-chromosomal translocations..."

	# Setup results dir, * add ignore is present
	if [ -d "$outputDir" ]; then
  		echo "Results folder already exists"
	else
		echo "Creating output results folder"
		mkdir $outputDir
	fi

	# bam2cfg.pl -g -h 120925_M00368_0055_A000000000-A1EY4_CGATGT_A_S1.bam > *.cfg
        echo "Building configuration file for BreakDancer..."
	bam2cfg.pl -g -h $fastq_base_name".bam" > $fastq_base_name".bam.cfg"

	# breakdancer-max
	echo "Running BreakDancer..."
	echo "Output in path:" $outputDir"/"$fastq_base_name".breakdancer.ctx"
	breakdancer_max $fastq_base_name".bam.cfg" > $outputDir"/"$fastq_base_name".breakdancer.ctx"

}

 

########## Main logic ########## 

# Temp hard set instead of passing comd line args, skip control of args with temp logic skip.
r1=120925_M00368_0055_A000000000-A1EY4_CGATGT_A_S1_L001_R1_001.fastq.gz
r2=120925_M00368_0055_A000000000-A1EY4_CGATGT_A_S1_L001_R2_001.fastq.gz

# Remove commentsif using cmd line args. Remove above hardset filenames
#Check input for 2 filein
#if [ $# -ne 2 ]; then 
#	echo "Enter 2 fastq files"
#	exit 1
#else
	# Set fastq names
	fastq_r1=$r1	#$1
	fastq_r2=$r2	#$2

	fastq_r1_base=`echo $fastq_r1 | sed 's/.fastq.gz//'` # core r1 filename
	fastq_r2_base=`echo $fastq_r2 | sed 's/.fastq.gz//'` # core r2 filename
	fastq_base_name=`echo $fastq_r1 | sed 's/_L001_R.*//'` # core filename

	# 1.Run fastqQC pre trim
	#func_fastqc $fastq_r1 $fastq_r2
	
	# 2.Run trimming. Creates trimmed vars.
	func_trimmomatic $fastq_r1 $fastq_r2
	
	# 3.Run fastQC post trim. Note variables assigned in trimming function.
	#func_fastqc $fastq_r1_trimmed $fastq_r2_trimmed

	# 4.Run mapping and QC
	#func_bwa $fastq_r1_trimmed $fastq_r2_trimmed
	func_mapping_stats
	
	# 5.Run CREST
	#func_crest

	# 6.Run BREAKDANCER
	func_breakdancer	

	echo "Done"

# Remove comment if using cmd args loop str
#fi


########## End ##########
