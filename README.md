# p_cml

# Dependencies on CTX1

CREST module on PBS does not include teh Perl module Bio:DB:Sam, or the assembler CAP3

Install locally using anaconda instead, then CREST conda package which installs crest and all dependencies (blat, cap3, bioperl, samtools) 

# Load conda module first and set up personal installation (if not previously)
module load anaconda3/personal
anaconda-setup

# Install CREST conda package
conda install crest 

Note, might also need to use samtools pre v1.0 (eg.v0.1.17)?
