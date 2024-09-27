	OVERVIEW: convert_5P_CAGE_wig_files_to_first_base_of_read.pl shifts CAGE wiggle file 
	outputs generated from STAR alignments using the "read1 5p" option to the first base of the read 
	alignment. This is necessary to properly link resulting CAGE peak files to ONT 5' ends. 

    
	SPECIAL NOTE FOR NEGATIVE STRAND WIGGLE FILES:
		wiggle files representing the negative strand must have negative values. If this is not
		the case, simply convert all signal values to a negative value before inputting into this script. 

    

	Inputs:
		wiggle files separated by spaces will convert each individually
    
    
	Option:
		-h help
    
    
	Usage: convert_5P_CAGE_wig_files_to_first_base_of_read.pl [INPUTS] file1.wig file2.wig file3.wig ...
    
    
	Example: perl convert_5P_CAGE_wig_files_to_first_base_of_read.pl PATH/file1.wig PATH/file2.wig PATH/file3.wig
