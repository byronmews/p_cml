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
gfServer start localhost 50000 /the/full/path/to/hg38.2bit
```


Many other dependencies:
	requires ref in 2bit format(fatotwobit)
	gfServer
	CREST uses its own perl Tree module(within install dir)

Make sure when running CREST.pl all modules are within @INC path
 

Note, path to CREST.pl must be the same full path as used to start gfServer
