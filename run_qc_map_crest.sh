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
module load crest/1.0 
module load blat/35
module load samtools/1.2
module load bioperl/1.6.0
module load picard/2.6.0

########## Static vars ##########
DBS="/home/malikian/bp_mapping/dbs"


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

#
# CREST
#func_crest() {
#	# Crest
#
#}
 

########## Main logic ########## 

# Temp hard set instead of passing comd line args, skip control of args with temp logic skip.
r1=120925_M00368_0055_A000000000-A1EY4_CGATGT_A_S1_L001_R1_001.fastq.gz
r2=120925_M00368_0055_A000000000-A1EY4_CGATGT_A_S1_L001_R2_001.fastq.gz

# Uncomment the below args handling if removing filenames hardset above.
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
	func_fastqc $fastq_r1 $fastq_r2
	
	# 2.Run trimming
	func_trimmomatic $fastq_r1 $fastq_r2
	
	# 3.Run fastQC post trim. Note variables assigned in trimming function.
	func_fastqc $fastq_r1_trimmed $fastq_r2_trimmed

	# 4.Run mapping
	func_bwa $fastq_r1_trimmed $fastq_r2_trimmed
	
	# 5.Run CREST
	#func_crest

	echo "Done"

# Remove comment if using cmd args
#fi


########## End ##########

