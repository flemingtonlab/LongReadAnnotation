#!/usr/bin/perl

use warnings;
use strict;  

my $file;

sub program_info {
    print "\n\tOVERVIEW: LR_3prime_end_summation.pl takes an input bed12\n\tlong read file and sums the number of 3' ends at each genomic position.\n\tOutput is in wiggle format.\n\n\tInputs:\n\t\t-bed\tbed12 full length long read file\n\n\tOption:\n\t\t\-h help
    
    \n\n\tUsage: perl LR_3prime_end_summation.pl [INPUTS] -bed <bed12 full length long read file>\n\n\tExample: perl LR_3prime_end_summation.pl -bed PATH/LongRead_file.bed\n\n";
    exit;
}

sub options {
    if (scalar @ARGV == 0) {
        program_info;
        exit;
    }    
    for (my $i=0; $i < scalar(@ARGV); $i++) {
        if ($ARGV[$i] eq "\-bed") {
            $file = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-h") {
            program_info;
            exit;
        }
    }
}

sub qc {
    if (not defined($file)) {
        print "\nLong read file not defined!\n";
        program_info;
        exit;
    }
}

sub ONT_3prime_end_summaries {
    my $cleaned_file_name = $file;
    $cleaned_file_name =~ s/\.bed$//g;
    my $file_ends = $cleaned_file_name."_ends";
    my $file_ends_sorted = $cleaned_file_name."_ends_sorted";
    my $file_ends_sorted_summed_pos = $cleaned_file_name."_ONT_3p_end_distribution_pos.wig";
    my $file_ends_sorted_summed_neg = $cleaned_file_name."_ONT_3p_end_distribution_neg.wig";
    print "\nOutputting 3' ends in bed format...\n\n";
    open (INF, "<$file") or die "couldn't open input file";
    open (OUT, ">$file_ends") or die "couldn't open output file";
    while (my $line = <INF>) {
	    chomp($line);
        my @split_line = split("\t", $line);
        if($split_line[5] eq "\+") {
            print OUT $split_line[0], "\t", $split_line[2]-1, "\t", $split_line[2]-1, "\tblank\t1\t", $split_line[5], "\n";
        }
        elsif($split_line[5] eq "\-") {
            print OUT $split_line[0], "\t", $split_line[1], "\t", $split_line[1], "\tblank\t1\t", $split_line[5], "\n";
        }
    }
    close(INF);
    close(OUT);

    print "Sorting bed file...\n\n";
    `sort -V -k 1,1 -k 2,2n -k 3,3n -k 6,6 $file_ends > $file_ends_sorted`;

    print "Summing counts for each coordinate...\n\n";
    open (INF, "<$file_ends_sorted") or die "couldn't open input file";
    open (OUT1, ">$file_ends_sorted_summed_pos") or die "couldn't open output file";
    open (OUT2, ">$file_ends_sorted_summed_neg") or die "couldn't open output file";
    my $previous_line_chr = "null";
    my $previous_line_coord = 0;
    my $previous_strand = "null";
    my $sum = 0;
    while (my $line = <INF>) {
	    chomp($line);
        my @split_line = split("\t", $line);
        if($previous_line_chr ne $split_line[0]) {
            
            if ($previous_strand eq "+") {
                print OUT1 $previous_line_coord+1, "\t", $sum, "\n";
                if (eof(INF)) {
                    print OUT1 $split_line[1]+1, "\t", $split_line[4], "\n";
                }
            }
            elsif ($previous_strand eq "-") {
                print OUT2 $previous_line_coord+1, "\t", -$sum, "\n";
                if (eof(INF)) {
                    print OUT2 $split_line[1]+1, "\t", -$split_line[4], "\n";
                }
            }
            print OUT1 "variableStep chrom=", $split_line[0], "\n";
            print OUT2 "variableStep chrom=", $split_line[0], "\n";
            $previous_line_chr = $split_line[0];
            $previous_line_coord = $split_line[1];
            $previous_strand = $split_line[5];
            $sum = 1;
        }
        elsif($split_line[0] eq $previous_line_chr and $split_line[1] == $previous_line_coord and $split_line[5] eq $previous_strand) {
            $sum = $sum + $split_line[4];
            if (eof(INF)) {
                if ($previous_strand eq "+") {
                    print OUT1 $previous_line_coord+1, "\t", $sum, "\n";
                }
                elsif ($previous_strand eq "-") {
                    print OUT2 $previous_line_coord+1, "\t", -$sum, "\n";
                }
            }
        }
        else {
            if ($previous_strand eq "+") {
                print OUT1 $previous_line_coord+1, "\t", $sum, "\n";
                if (eof(INF)) {
                    print OUT1 $split_line[1]+1, "\t", $split_line[4], "\n";
                }
            }
            elsif ($previous_strand eq "-") {
                print OUT2 $previous_line_coord+1, "\t", -$sum, "\n";
                if (eof(INF)) {
                    print OUT2 $split_line[1]+1, "\t", -$split_line[4], "\n";
                }
            }
            $previous_line_chr = $split_line[0];
            $previous_line_coord = $split_line[1];
            $previous_strand = $split_line[5];
            $sum = 1;
        }
    }
    `rm $file_ends`;
    `rm $file_ends_sorted`
}
options;
qc;
ONT_3prime_end_summaries;