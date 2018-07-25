#!/bin/bash

# Name: Basic bash pipeline - convert to a workflow engine later
# Author: Graham Rose
# Date: 14-06-2018
# Version 0.1
# Description: Breakpoint mapping analysis
#----------------------------------------------------

########## Static vars ##########
outputDir="results"
refGenome="hg38.fa" # indexed previously
refGenomePath="/home/malikian/bp_mapping/dbs"


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
	
	trimmomatic PE -phred33 -threads 4 -trimlog $fastq_base_name".trimmomatic.log" $1 $2 $fastq_r1_base"_PE.fastq" $fastq_r1_base"_SE.fastq" $fastq_r2_base"_PE.fastq" $fastq_r2_base"_SE.fastq" ILLUMINACLIP:$refGenomePath/clip.fa:2:30:10 LEADING:20 TRAILING:20 MINLEN:40 	

	# New vars for trimmed PE fastq files
	fastq_r1_trimmed=$fastq_r1_base"_PE.fastq"
	fastq_r2_trimmed=$fastq_r2_base"_PE.fastq"

}

# Mapping and parse BAM. Genome fasta and bwa index in refGenomePath variable declared at top
func_bwa_mem() {
	echo "Running bwa mem..."
        echo "PE, 8 threads"
	echo "Ouput SAM filename $fastq_base_name"

	# Create syminks to genome index
	echo "Creating symlinks to indexed hg38 genome at" $refGenomePath
	ln -s $refGenomePath"/"$refGenome* .

	# Run bwa mem
	bwa mem -t 8 $refGenome $1 $2 > $fastq_base_name".sam"

	# SAM>BAM
	echo "Convert to BAM..."
	picard SortSam VALIDATION_STRINGENCY=LENIENT INPUT=$fastq_base_name".sam" OUTPUT=$fastq_base_name".bam" SORT_ORDER=coordinate
	
	samtools index $fastq_base_name".bam"
}

# Mapping and parse BAM. Genome fasta and bwa index in refGenomePath variable declared at top
func_bwa_sampe() {
        echo "Running bwa sampe..."
        echo "PE, 8 threads"
        echo "Ouput SAM filename $fastq_base_name"

        # Create syminks to genome index
        echo "Creating symlinks to indexed hg38 genome at" $refGenomePath
        ln -s $refGenomePath"/"$refGenome.* .

	# Run index
	bwa aln $refGenome $1 > $fastq_base_name"_R1.sai"
	bwa aln $refGenome $2 > $fastq_base_name"_R2.sai"

        # Run bwa sampe
	bwa sampe $refGenome $fastq_base_name"_R1.sai" $fastq_base_name"_R2.sai" $1 $2 > $fastq_base_name".sam"

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


# CREST. Using anaconda package
func_crest() {

	echo "Running CREST - Clipping REveals STructure..."
	echo "Using anaconda3/personal for Bio/DB/Sam.pm module dependency, and custom modules at:"
	echo "		home/malikian/bp_mapping/CREST_modules/"

	echo "Get soft-clipping positions..."
	echo "Extracting chromosomes:"

	arr=( chr8 chr9 chr21 chr22 chr23 )

	for element in "${arr[@]}"
	do	
		echo $element

		perl -I ~/anaconda3/lib/perl5/site_perl/5.22.0/ ~/anaconda3/bin/extractSClip.pl -i $fastq_base_name".bam" --ref_genome $refGenome -r $element
	done
	
	# Starting blat server
	echo "Starting blat server using gfServer from blat suite"	
	#gfServer start localhost 8666 $refGenome".2bit"

	# Concatenate split chrs
	cat "$fastq_base_name.bam".*".cover" > $fastq_base_name".bam.cover"
	
	echo "SV detection..."
	perl -I ~/anaconda3/lib/perl5/site_perl/5.22.0/ ~/anaconda3/bin/CREST.pl -f $fastq_base_name".bam.cover" -d $fastq_base_name".bam" --ref_genome $refGenome -t hg38.fa.2bit --blatserver localhost --blatport 8111

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
	func_fastqc $fastq_r1 $fastq_r2
	
	# 2.Run trimming. Creates trimmed vars.
	func_trimmomatic $fastq_r1 $fastq_r2
	
	# 3.Run fastQC post trim. Note variables assigned in trimming function.
	func_fastqc $fastq_r1_trimmed $fastq_r2_trimmed

	# 4.Run mapping and QC
	func_bwa_sampe $fastq_r1_trimmed $fastq_r2_trimmed
	#func_mapping_stats
	
	# 5.Run CREST
	func_crest

	# 6.Run BREAKDANCER
	#func_breakdancer	

	echo "Done"

# Remove comment if using cmd args loop str
#fi


########## End ##########

