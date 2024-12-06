OVERVIEW: invert_wigs.pl takes one or more input wig files
	and converts the coverage value to -n.
	Can input multiple second strand wig files at the same time.

	Inputs:
		-wig	one or more input wig files

	Option:
		-h help
    
    

	Usage: perl invert_wigs.pl [INPUTS] -wig <wiggle file(s)>

	Example: perl invert_wigs.pl -wig PATH/file1_neg_strand.wig,PATH/file2_neg_strand.wig,PATH/file3_neg_strand.wig
