# p_cml

## Dependencies on CTX1 node @ HPC

CREST module on PBS does not include teh Perl module Bio:DB:Sam, or the assembler CAP3

Install locally using anaconda instead, then CREST conda package which installs crest and all dependencies (blat, cap3, bioperl, samtools) 

Load conda module first and set up personal installation (if not previously)
```
module load anaconda3/personal
anaconda-setup
```

Install CREST conda package
```
conda install crest
```

Install perl-bio-db-sam dependency
```
conda install -c bioconda perl-bio-db-sam
```

gfServer setup and use module on HPC. Manage server.
```
gfServer start localhost 8111 hg38.fa.2bit
```


Many other dependencies:
	requires ref in 2bit format(fatotwobit)
	gfServer
	CREST own perl Tree module(within install)

Make sure when running CREST.pl all modules are within @INC path
 

Note, might also need to use samtools pre v1.0 (eg.v0.1.17)?
