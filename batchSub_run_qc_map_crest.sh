#!/bin/bash
#PBS -lwalltime=12:00:00
#PBS -l select=1:ncpus=8:mem=16gb
# Email when job ends 
#PBS -m e
#
# Job name (default is name of pbs script file)
#PBS -N crest_analysis
########## End of PBS scripting

# Usage: Batch script to run a serial job under SGE.

# Run info
echo running on `hostname`
echo Node file is  $PBS_NODEFILE
echo The assigned nodes are:
cat $PBS_NODEFILE

echo --------------------------------------

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
########## End of load modules


#Run
cd $PBS_O_WORKDIR
bash run_qc_map_breaks.sh


########## Done ##########



