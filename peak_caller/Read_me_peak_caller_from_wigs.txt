	OVERVIEW: peak_caller_from_wigs.pl is designed to facilitate the use of custom input parameters to accommodate different goals as well as distinct characteristics of different organisms, such as viruses, which have higher densities of promoters. 

	Inputs:
		-w	option facilitates the input of multiple wig files for an analysis (input  			wigs must be of same strand).
		-mw	max peak width 
		-fva	fraction of max peak value for adjacent signals
		-mspd	minimum value for a single position depth signal to seed a peak
		-s	strand of wig files analyzed

	Option:
		-h help
    
    

	Usage: perl peak_caller_from_wigs.pl [INPUTS] -w <file1.wig,file2.wig,file3.wig,...> -mw <max peak width> -fva <minimum fraction of main peak for adjacent signals to be included in peak> -mspd <minimum single position depth to seed peak> -s <strand of wig files (+ or -)>

	Example: perl peak_caller_from_wigs.pl -w PATH/file1.wig,PATH/file2.wig,PATH/file3.wig -mw 8 -fva 0.2 -mspd 200 -s + 
