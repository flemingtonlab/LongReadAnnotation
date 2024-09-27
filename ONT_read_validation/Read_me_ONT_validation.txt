OVERVIEW: ONT_validation.pl is designed to validate ONT reads using 1) splice junction data from short read alignments, 2) CAGE 5' end peaks, 3) ONT read 3' end peaks, and an input long read file in bed12 format. 
    
    

	Inputs:
		-d	Indicate whether dependency on BALF2, vPIC, or OriLyt is desired.
		-mcde	Minimum CAGE peak depth
		-mcdi	Maximum distance between CAGE peak and start of ONT read
		-3Pp	3' peak file (bed6 format)
		-3Pde	Minimum 3' peak read depth
		-3Pdi	Maximum 3' distance between ONT end and 3' peak
		-minSJ	Minimum number of short read splice junctions to validate spliced ONT reads
		-SJt	Splice junction tab file from STAR alignment of short reads
		-ONT	Input long read bed12 file
    
    

	Option:
		-h help
    
    

	Usage: perl ONT_validation.pl [INPUTS] -d <BALF2 or vPIC or OriLYt> -mcde <Minimum CAGE peak depth> -mcdi <Maximum distance between CAGE peak and start of ONT read> -3Pp <3' peak file> -3Pde <Minimum read depth of 3' peak clusters> -3Pdi <Maximum 3' distance between ONT end and 3' peaks> -minSJ <Minimum number of validating short read splice junctions detected> -SJt <Splice junction tab file from STAR alignment of short reads> -ONT <Input long read bed12 file>
    
    

	Example: perl ONT_validation.pl -d BALF2 -mcde 1000 -mcdi 2 -3Pp 3Ppeakfile.bed -3Pde 400 -3Pdi 10 -minSJ 200 -SJt SJfile.tab -ONT ONTfile.bed

