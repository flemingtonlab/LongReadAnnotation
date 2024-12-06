OVERVIEW: get_known_ORF_start_sites.pl takes an input bed12 transcript
	annotation file and pulls out all known annotated ATG start sites for inputting into
	LR_validate.pl if desired.

	Inputs:
		-bed	bed12 annotation file

	Option:
		-h help
    
    

	Usage: perl get_known_ORF_start_sites.pl [INPUTS] -bed <bed12 annotation file>

	Example: perl get_known_ORF_start_sites.pl -bed PATH/annotation.bed
