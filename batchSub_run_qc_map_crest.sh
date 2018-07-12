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
module load picard/2.6.0
module load breakdancer/1.2
module load anaconda3/personal # for CRESt and dependencies
########## End of load modules


#Run
cd $PBS_O_WORKDIR
bash run_qc_map_breaks.sh


########## Done ##########



