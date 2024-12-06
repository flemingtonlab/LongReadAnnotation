#!/usr/bin/perl

use warnings;
use strict;  
use File::Basename;

my $annotation_file;
my $SS_count = 0;

sub program_info {
    print "\n\tOVERVIEW: get_known_ORF_start_sites.pl takes an input bed12 transcript\n\tannotation file and pulls out all known annotated ATG start sites for inputting into\n\tLR_validate.pl if desired.\n\n\tInputs:\n\t\t-bed\tbed12 annotation file\n\n\tOption:\n\t\t\-h help
    
    \n\n\tUsage: perl get_known_ORF_start_sites.pl [INPUTS] -bed <bed12 annotation file>\n\n\tExample: perl get_known_ORF_start_sites.pl -bed PATH/annotation.bed\n\n";
    exit;
}

sub options {
    if (scalar @ARGV == 0) {
        program_info;
        exit;
    }    
    for (my $i=0; $i < scalar(@ARGV); $i++) {
        if ($ARGV[$i] eq "\-bed") {
            $annotation_file = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-h") {
            program_info;
            exit;
        }
    }
}

sub qc {
    if (not defined($annotation_file)) {
        print "\nAnnotation file not defined!\n";
        program_info;
        exit;
    }
}

sub get_ATGs {
    open (INF, "<$annotation_file") or die "couldn't open input file";
    open (OUT, ">$annotation_file.temp") or die "couldn't open output file";

    while (my $line = <INF>) {
	    chomp($line);
        if($line !~ m/^chr/g) {
            next;
        }
        my @split_line = split("\t", $line);
        if($split_line[7] - $split_line[6] > 10) {
            if(eof) {
                if($split_line[5] eq "+" and $split_line[1] ne $split_line[6]) {
                    $SS_count++;
                    print OUT $split_line[0], "\t", $split_line[6], "\t", $split_line[6], "\tStart_site_", $SS_count, "\t\.\t",$split_line[5];
                }
                elsif($split_line[5] eq "-" and $split_line[2] ne $split_line[7]) {
                    $SS_count++;
                    print OUT $split_line[0], "\t", $split_line[7], "\t", $split_line[7], "\tStart_site_", $SS_count, "\t\.\t", $split_line[5];
                }
            }
            elsif($split_line[5] eq "+" and $split_line[1] ne $split_line[6]) {
                $SS_count++;
                print OUT $split_line[0], "\t", $split_line[6], "\t", $split_line[6], "\tStart_site_", $SS_count, "\t\.\t",$split_line[5], "\n";;
            }
            elsif($split_line[5] eq "-" and $split_line[2] ne $split_line[7]) {
                $SS_count++;
                print OUT $split_line[0], "\t", $split_line[7], "\t", $split_line[7], "\tStart_site_", $SS_count, "\t\.\t", $split_line[5], "\n";;
            }
        }
    }
    close(INF);
    close(OUT);

    my $temp_file = $annotation_file.".temp";
    my $sorted_file = $annotation_file.".known_ORF_start_sites.bed";
    `sort -V -k 1,1 -k 2,2n -k 3,3n $temp_file > $sorted_file`;
    `rm $temp_file`;
}

options;
qc;
get_ATGs;