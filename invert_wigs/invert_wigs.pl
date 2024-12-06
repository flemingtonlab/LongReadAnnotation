#!/usr/bin/perl

use warnings;
use strict;  

my @files;
sub program_info {
    print "\n\tOVERVIEW: invert_wigs.pl takes one or more input wig files\n\tand converts the coverage value to -n.\n\tCan input multiple second strand wig files at the same time.\n\n\tInputs:\n\t\t-wig\tone or more input wig files\n\n\tOption:\n\t\t\-h help
    
    \n\n\tUsage: perl invert_wigs.pl [INPUTS] -wig <wiggle file(s)>\n\n\tExample: perl invert_wigs.pl -wig PATH/file1_neg_strand.wig,PATH/file2_neg_strand.wig,PATH/file3_neg_strand.wig\n\n";
    exit;
}

sub options {
    if (scalar @ARGV == 0) {
        program_info;
        exit;
    }    
    for (my $i=0; $i < scalar(@ARGV); $i++) {
        if ($ARGV[$i] eq "\-wig") {
            my $input = $ARGV[$i+1];
            @files = split("\,", $input);
        }
        elsif ($ARGV[$i] eq "\-h") {
            program_info;
            exit;
        }
    }
}

sub qc {
    if (scalar(@files) == 0) {
        print "\ninput wiggle files not defined!\n";
        program_info;
        exit;
    }
}

sub invert_wigs {
    foreach my $file(@files) {
        my $output_file = $file;
        $output_file =~ s/out.wig/out./g;
        $output_file = $output_file."negative_values.wig";
        print $output_file, "\n";
        open (INF, "<$file") or die "couldn't open input file";
        open (OUT, ">$output_file") or die "couldn't open output file";
        while (my $line = <INF>) {
	        chomp($line);
            my @split_line = split("\t", $line);
            if($line =~ m/^variable/g) {
                print OUT $line, "\n";
            }
            else {
                print OUT $split_line[0], "\t", -$split_line[1], "\n";
            }
        }
        close(INF);
        close(OUT);
    }
}
options;
qc;
invert_wigs;