# LongReadAnnotation 
![Validation schematic](images/2_long_read_validation_schematic.jpg)

# Citation

Erik K Flemington

Dinh Truong Nguyen

# Requirements
Perl 5

# Dependencies
File::Basename

# Installation
```
git clone https://github.com/flemingtonlab/LongReadAnnotation.git
```


# File formats

| File    | Description     |
|:---------------|:---------------:|
|Wiggle  | Standard wiggle format |
| BED12 | Standard BED12 format  |
| BED6 | Standard BED6 format  |
| Genome fasta | Standard genome fasta (can be wrapped or unwrapped)|

# LongReadAnnotation Pipeline
1) identify start site clusters
	
	Call 5' peaks from wig files (use “-CAGE y” option which shifts output 1bp downstream for CAGE data because STAR aligner CAGE output option (—outWigType read1_5p) outputs signals 1bp upstream from start of read) 
	- perl /PATH/peak_caller_from_wigs.pl -w PATH/MC1_Unique.str1.out_chr1.wig,/PATH/MC2_Unique.str1.out_chr1.wig,/PATH/MC4_Unique.str1.out_chr1.wig -mw 8 -fva 0.2 -mspd 10 -s + -CAGE y

		Negative strand input wigs values must be negative! (If coverage values for negative strand wigs are positive, use invert_wigs.pl to change sign)
	- perl /PATH/peak_caller_from_wigs.pl -w PATH/MC1_Unique.str2.out.negative_values_chr1.wig,/PATH/MC2_Unique.str2.out.negative_values_chr1.wig,/PATH/MC4_Unique.str2.out.negative_values_chr1.wig -mw 8 -fva 0.2 -mspd 10 -s - -CAGE y


	Concatenate positive and negative strand bed files
	- cat /PATH/CAGE_peaks_positive_strand.bed /PATH/CAGE_peaks_negative_strand.bed > /PATH/CAGE_peaks_positive_plus_negative_strand.bed

	Sort concatenated positive plus negative strand bed files (not essential but is a good configuration for the file)
	- sort -V -k 1,1 -k 2,2n -k 3,3n /PATH/CAGE_peaks_positive_plus_negative_strand.bed > /PATH/CAGE_peaks_positive_plus_negative_strand_sorted.bed
