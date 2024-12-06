OVERVIEW: LR_validate.pl is designed to validate mapped long read data using 1) splice junction
	data from short read alignments, 2) CAGE 5' end peaks, 3) ONT read 3' end peaks, and
	4) an input mapped long read file in bed12 format. 
    
    

	Inputs:
		-5Pp	CAGE peak file (.bed10 format CAGE peak output from peak_caller_from_wigs.pl).
		-mcde	Minimum CAGE peak depth
		-mcdi	Maximum distance between CAGE peak and start of ONT read
		-3Pp	3' peak file (bed6 format)
		-3Pde	Minimum 3' peak read depth
		-3Pdi	Maximum 3' distance between ONT end and 3' peak
		-minSJ	Minimum number of short read splice junctions to validate spliced ONT reads
		-SJt	Splice junction tab file from STAR alignment of short reads
		-f	genome fasta file
		-LR	Input long read bed12 file
    
    

	Options:
		-ATG	known ATG starts site file
		-h	help
    
    

	Usage: perl LR_validate.pl [INPUTS] -5Pp <5' CAGE peak file (bed10 format from peak_caller_from_wigs.pl)> -mcde <Minimum CAGE peak depth> -mcdi <Maximum distance between CAGE peak and start of ONT read> -3Pp <3' peak file> -3Pde <Minimum read depth of 3' peak clusters> -3Pdi <Maximum 3' distance between ONT end and 3' peaks> -minSJ <Minimum number of validating short read splice junctions detected> -SJt <Splice junction tab file from STAR alignment of short reads> -ATG <known ATG start site file (bed6) (default = no known starts)> -f <genome fasta file> -LR <Input long read bed12 file>
    
    

	Example: perl validate_ONT.pl -5Pp /PATH/5Ppeakfile.bed -mcde 10 -mcdi 2 -3Pp /PATH/3Ppeakfile.bed -3Pde 10 -3Pdi 10 -minSJ 2 -SJt /PATH/SJfile.tab -f /PATH/genome.fa -LR /PATH/LRfile.bed
