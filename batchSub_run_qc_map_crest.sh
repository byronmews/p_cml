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

#Run
cd $PBS_O_WORKDIR
bash run_qc_map_crest.sh


########## Done ##########



