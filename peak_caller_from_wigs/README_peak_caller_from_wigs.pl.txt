OVERVIEW: peak_caller_from_wigs.pl is designed to facilitate the use
	of custom input parameters to accommodate different goals as well as
	distinct characteristics of different organisms, such as viruses,
	which have higher densities of promoters. 

	Important NOTE: You must input one or more of only positive strand or one or more
	of only negative strand wiggle files. After identifying peaks for the positive strand
	and the negative strand, the files can then be concatenated.

	Inputs:
		-w	one or more wig files
			(input wigs must be of same strand).
		-mw	maximum peak width 
		-fva	fraction of maximum peak value for adjacent signals
		-mspd	minimum value for a single position depth signal to seed a peak
		-s	strand of wig files analyzed
		-CAGE	shifts output 1bp downstream for CAGE data (e.g. STAR aligner
			"--outWigType read1_5p" option outputs 1bp upstream
			from start of read) 

	Option:
		-h help
    
    

	Usage: perl peak_caller_from_wigs.pl [INPUTS] -w <file1.wig,file2.wig,file3.wig,...> -mw <max peak width> -fva <minimum fraction of main peak for adjacent signals to be included in peak> -mspd <minimum single position depth to seed peak> -s <strand of wig files (+ or -)> -CAGE <y or n (default is n)>

	Example: perl peak_caller_from_wigs.pl -w PATH/file1.wig,PATH/file2.wig,file3.wig -mw 8 -fva 0.2 -mspd 10 -s + -CAGE y


