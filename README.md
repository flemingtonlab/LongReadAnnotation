# LongReadAnnotation 
![Validation schematic](images/2_long_read_validation_schematic.jpg)

## Project Description
LongReadAnnotation is a tool for annotating and validating long reads from sequencing data. It includes pipelines for identifying start and end site clusters and validating the long reads.

## Citation

Erik K Flemington, Dinh Truong Nguyen

## Requirements
- Perl 5

## Dependencies
- File::Basename

##Installation
1. Clone the repository:
```sh
git clone https://github.com/flemingtonlab/LongReadAnnotation.git
```
2. Install the dependencies (ensure Perl and required modules are installed).


## File Formats

| File          | Description                                      |
|:--------------|:------------------------------------------------:|
| Wiggle        | Standard wiggle format                           |
| BED12         | Standard BED12 format                            |
| BED6          | Standard BED6 format                             |
| Genome fasta  | Standard genome fasta (can be wrapped or unwrapped) 
|

## LongReadAnnotation Pipeline
### 1. Identify Start Site Clusters
1. Call 5' peaks from wig files:
   ```sh
   perl /PATH/peak_caller_from_wigs.pl -w PATH/MC1_Unique.str1.out_chr1.wig,PATH/MC2_Unique.str1.out_chr1.wig,PATH/MC4_Unique.str1.out_chr1.wig -mw 8 -fva 0.2 -mspd 10 -s + -CAGE y
   ```
2. Ensure negative strand input wigs values are negative! (If coverage values for negative strand wigs are positive, use invert_wigs.pl to change sign):
   ```sh
   perl /PATH/peak_caller_from_wigs.pl -w PATH/MC1_Unique.str2.out.negative_values_chr1.wig,PATH/MC2_Unique.str2.out.negative_values_chr1.wig,PATH/MC4_Unique.str2.out.negative_values_chr1.wig -mw 8 -fva 0.2 -mspd 10 -s -
   ```
3. Concatenate and sort bed files:
   ```sh
   cat /PATH/CAGE_peaks_positive_strand.bed /PATH/CAGE_peaks_negative_strand.bed > /PATH/CAGE_peaks_positive_plus_negative_strand.bed
   sort -V -k 1,1 -k 2,2n -k 3,3n /PATH/CAGE_peaks_positive_plus_negative_strand.bed > /PATH/CAGE_peaks_positive_plus_negative_strand_sorted.bed
   ```

### 2. Identify 3' Site Clusters
1. Summarize 3' LR end coverage:
   ```sh
   perl /PATH/LR_3prime_end_summation.pl -bed /PATH/LR.bed
   ```
2. Identify and sort 3' end clusters:
   ```sh
   perl /PATH/peak_caller_from_wigs.pl -w PATH/MC1_3p_chr1.wig,PATH/MC2_3p_chr1.wig,PATH/MC4_3p_chr1.wig -mw 8 -fva 0.2 -mspd 10 -s +
   perl /PATH/peak_caller_from_wigs.pl -w PATH/MC1_3p_negative_values_chr1.wig,PATH/MC2_3p_negative_values_chr1.wig,PATH/MC4_3p_negative_values_chr1.wig -mw 8 -fva 0.2 -mspd 10 -s -
   cat /PATH/3p_peaks_positive_strand.bed /PATH/3p_peaks_negative_strand.bed > /PATH/3p_peaks_positive_plus_negative_strand.bed
   sort -V -k 1,1 -k 2,2n -k 3,3n /PATH/3p_peaks_positive_plus_negative_strand.bed > /PATH/3p_peaks_positive_plus_negative_strand_sorted.bed
   ```

### 3. Validate Long Reads
1. Perform long read validation:
   ```sh
   perl /PATH/LR_validate.pl -5Pp /PATH/MU_5P_CAGE_peaks_chr1.bed -mcde 10 -mcdi 2 -3Pp /PATH/MU_3P_peaks_chr1.bed -3Pde 10 -3Pdi 10 -minSJ 1 -SJt /PATH/MU-SJ.out.tab -f /PATH/hg38_chr1_first_por...
   ```
2. Validate with previously identified ATG start sites:
   ```sh
   perl /PATH/LR_validate.pl -5Pp /PATH/MU_5P_CAGE_peaks_chr1.bed -mcde 10 -mcdi 2 -3Pp /PATH/MU_3P_peaks_chr1.bed -3Pde 10 -3Pdi 10 -minSJ 1 -SJt /PATH/MU-SJ.out.tab -f /PATH/hg38_chr1_first_portion.fa  -LR /PATH/MU_LR_fullLength.merged_1million.bed -ATG /PATH/LR_validate/test_data/hg38_chr1_known_ORF_start_sites.bed
   ```
