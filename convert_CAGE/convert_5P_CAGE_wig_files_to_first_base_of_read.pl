#!/usr/bin/perl

use warnings;
use strict;  
use File::Basename;

# NOTE: Requires that negative strand wig files have negative values!!!!!

sub program_info {
    print "\n\tOVERVIEW: convert_5P_CAGE_wig_files_to_first_base_of_read.pl shifts CAGE wiggle file \n\toutputs generated from STAR alignments using the \"read1 5p\" option to the first base of the read \n\talignment. This is necessary to properly link resulting CAGE peak files to ONT 5' ends. 

    \n\tSPECIAL NOTE FOR NEGATIVE STRAND WIGGLE FILES:\n\t\twiggle files representing the negative strand must have negative values. If this is not\n\t\tthe case, simply convert all signal values to a negative value before inputting into this script. 

    \n\n\tInputs:\n\t\twiggle files separated by spaces will convert each individually
    
    \n\tOption:\n\t\t\-h help
    
    \n\tUsage: convert_5P_CAGE_wig_files_to_first_base_of_read.pl [INPUTS] file1.wig file2.wig file3.wig ...
    
    \n\tExample: perl convert_5P_CAGE_wig_files_to_first_base_of_read.pl PATH/file1.wig PATH/file2.wig PATH/file3.wig\n\n";
    exit;
}

my @files = @ARGV;



sub options {
    if (scalar @ARGV == 0) {
        print "\n\n\tNO INPUT!\n\n";
        program_info;
        exit;
    }    
    for (my $i=0; $i < scalar @ARGV; $i++) {
        if ($ARGV[$i] eq "\-h") {
            program_info;
            exit;
        }
    }
}

sub convert_wigs {
    print "\nConverting wiggle files...\n\n"; 

    foreach my $file(@files) {
        my $output_file = $file;
        $output_file =~ s/\.wig//g;
        $output_file = $output_file."_CAGE_wig_converted_to_first_base_of_read\.wig";
        my $strand = "null";
        open(INF, "<$file") or die "couldn't open wig file";
        while(my $line = <INF>) {
            chomp($line);
            my @split_line = split("\t", $line);
            if($. == 3) {
                last;
            }
            elsif($. == 2) {
                if($split_line[1] > 0) {
                    $strand = "\+";
                }
                elsif($split_line[1] < 0) {
                    $strand = "\-";
                }
            }
        }
        close(INF);
        if($strand eq "\+") {
            open(INF, "<$file") or die "couldn't open wig file";
            open(OUT, ">$output_file") or die "couldn't open wig file";
            while(my $line = <INF>) {
                chomp($line);
                my @split_line = split("\t", $line);
                if($line =~ m/chrom/g) {
                    print OUT $line, "\n";
                }
                else {
                    print OUT $split_line[0]+1, "\t", $split_line[1], "\n";;
                }
            }
            close(INF);
            close(OUT);
        }
        elsif($strand eq "\-") {
            open(INF, "<$file") or die "couldn't open wig file";
            open(OUT, ">$output_file") or die "couldn't open wig file";
            while(my $line = <INF>) {
                chomp($line);
                my @split_line = split("\t", $line);
                if($line =~ m/chrom/g) {
                    print OUT $line, "\n";
                }
                else {
                    print OUT $split_line[0]-1, "\t", $split_line[1], "\n";;
                }
            }
            close(INF);
            close(OUT);
        }
        else {
            print "\n\nERROR!! Strand not determined!\n\n";
        }
    }
}
options;
convert_wigs;