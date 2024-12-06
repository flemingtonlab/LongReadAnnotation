OVERVIEW: LR_3prime_end_summation.pl takes an input bed12
	long read file and sums the number of 3' ends at each genomic position.
	Output is in wiggle format.

	Inputs:
		-bed	bed12 full length long read file

	Option:
		-h help
    
    

	Usage: perl LR_3prime_end_summation.pl [INPUTS] -bed <bed12 full length long read file>

	Example: perl LR_3prime_end_summation.pl -bed PATH/LongRead_file.bed
