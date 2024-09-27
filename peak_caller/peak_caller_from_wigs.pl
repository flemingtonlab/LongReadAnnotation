#!/usr/bin/perl

use warnings;
use strict;  
use File::Basename;

my @wig_files = ();
my $max_peak_width;
my $fraction_value_for_adjacent_peaks;
my $min_single_position_depth;
my $strand;

sub program_info {
    print "\n\tOVERVIEW: peak_caller_from_wigs.pl is designed to facilitate the use of custom input parameters to accommodate different \n\tgoals as well as distinct characteristics of different organisms, such as viruses, which have higher densities of promoters. \n\n\tInputs:\n\t\t-w\toption facilitates the input of multiple wig files for an analysis (input wigs must be of same strand).\n\t\t-mw\tmax peak width \n\t\t-fva\tfraction of max peak value for adjacent signals\n\t\t-mspd\tminimum value for a single position depth signal to seed a peak\n\t\t-s\tstrand of wig files analyzed\n\n\tOption:\n\t\t\-h help
    
    \n\n\tUsage: perl peak_caller_from_wigs.pl [INPUTS] -w <file1.wig,file2.wig,file3.wig,...> -mw <max peak width> -fva <minimum fraction of main peak for adjacent signals to be included in peak> -mspd <minimum single position depth to seed peak> -s <strand of wig files (+ or -)>\n\n\tExample: perl peak_caller_from_wigs.pl -w PATH/file1.wig,PATH/file2.wig,PATH/file3.wig -mw 8 -fva 0.2 -mspd 200 -s + \n\n";
    exit;
}

sub options {
    if (scalar @ARGV == 0) {
        program_info;
        exit;
    }
    for (my $i=0; $i < scalar @ARGV; $i++) {
        if ($ARGV[$i] eq "\-w") {
            @wig_files = split("\,", $ARGV[$i+1]);
        }
        elsif ($ARGV[$i] eq "\-mw") {
            $max_peak_width = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-fva") {
            $fraction_value_for_adjacent_peaks = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-mspd") {
            $min_single_position_depth = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-s") {
            $strand = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-h") {
            program_info;
            exit;
        }
    }
}

sub qc {
    if (scalar(@wig_files) == 0) {
        print "\nwiggle files not defined!\n";
        program_info;
        exit;
    }
    elsif (not defined($max_peak_width)) {
        print "\nMax peak width not defined!\n";
        program_info;
        exit;
    }
    elsif (not defined($fraction_value_for_adjacent_peaks)) {
        print "\nMinimum fraction of main peak amplitude to be included in defining peak width not defined!\n";
        program_info;
        exit;
    }
    elsif (not defined($min_single_position_depth)) {
        print "\nMinimum single position depth to seed peak not defined!\n";
        program_info;
        exit;
    }
    elsif (not defined($strand)) {
        print "\nStrand is not defined!\n";
        program_info;
        exit;
    }
}

my @all_positions = ();
my $all_positions_length = 0;

my @first_run_peaks = ();
my $first_run_peaks_length = 0;

my @second_run_peaks = ();
my $second_run_peaks_length = 0;

sub merge_wigs {
    print "\nMerging wiggle files...\n\n"; 
    my @merged_positions = ();
    foreach my $file(@wig_files) {
        my $chr_1;
        open(INF, "<$file") or die "couldn't open wig file";
        while(my $line = <INF>) {
            chomp($line);
            if ($line =~ m/chr/g) {
                my @split_line = split("\=", $line);
                $chr_1 = $split_line[1];
            }
            else {
                my @split_line = split("\t", $line);
                my $positions_bed_format = $chr_1."\t".$split_line[0]."\t".$split_line[0]."\t\.\t".$split_line[1]."\t".$strand;
                push(@merged_positions, $positions_bed_format);
            }
        }
        close(INF);
    }
    my $out_file = $wig_files[0]."merged";
    open(OUT, ">$out_file") or die "couldn't open wig file";
    print OUT join("\n", @merged_positions);
    close(OUT);
    my $sorted_merged = $wig_files[0]."sorted";
    `sort -V -k 1,1 -k 2,2n -k 3,3n $out_file > $sorted_merged`;
    `rm $out_file`;

    open(INF, "<$sorted_merged") or die "couldn't open sorted merged file";
    my $chr;
    my $coord1;
    my $total_reads;
    while(my $line = <INF>) {
        chomp($line);
        my @split_line = split("\t", $line);

        if($. == 1) {
            $chr = $split_line[0];
            $coord1 = $split_line[1];
            $total_reads = $split_line[4]
        }
        elsif($coord1 != $split_line[1] || $chr ne $split_line[0]) {
            my $out_line = $chr."\t".$coord1."\t".$coord1."\t\.\t".$total_reads."\t".$strand;
            push(@all_positions, $out_line);
            $chr = $split_line[0];
            $coord1 = $split_line[1];
            $total_reads = $split_line[4];
            if(eof(INF)) {
                my $out_line = $split_line[0]."\t".$split_line[1]."\t".$split_line[1]."\t\.\t".$split_line[4]."\t".$strand;
                push(@all_positions, $out_line);
            }
        }
        elsif($chr eq $split_line[0] and $coord1 == $split_line[1]) {
            $total_reads = $split_line[4] + $total_reads;
            if(eof(INF)) {
                my $out_line = $chr."\t".$coord1."\t".$coord1."\t\.\t".$total_reads."\t".$strand;
                push(@all_positions, $out_line);
            }
        }
    }
    close(INF);
    `rm $sorted_merged`;
    $all_positions_length = @all_positions;
}

sub determine_positions_with_min_single_position_depth {
    print "Extracting all genome positions with greater than minimum number of reads...\n\n"; 

    foreach my $line(@all_positions) {
        my @split_line = split("\t", $line);
        if(abs($split_line[4]) > $min_single_position_depth) {
            my $peak_line = $split_line[0]."\t".$split_line[1]."\t".$split_line[2]."\tPeak\t".$split_line[4]."\t".$strand;
            push(@first_run_peaks, $peak_line);
        }
    }
    $first_run_peaks_length = @first_run_peaks;
}

sub filter_first_run_peaks {
    print "Identifying max peaks within max distance...\n\n"; 

    for(my $i=0; $i < $first_run_peaks_length; $i++) {
        my @split_i = split("\t", $first_run_peaks[$i]);
        my $omit_count = 0;
        for(my $j=0; $j < $first_run_peaks_length; $j++) {
            my @split_j = split("\t", $first_run_peaks[$j]);
            if($split_i[0] eq $split_j[0]) {
                if($split_j[1] >= $split_i[1]-$max_peak_width and $split_j[1] <= $split_i[1]+$max_peak_width) {
                    if(abs($split_j[4]) > abs($split_i[4])) {
                        $omit_count++;
                    }
                }
            }
        }
        if($omit_count == 0) {
            push(@second_run_peaks, $first_run_peaks[$i]);
        }
    }
    $second_run_peaks_length = @second_run_peaks;
}

sub identify_clusters {
    print "Identify clusters...\n\n"; 
    my $outdirectory = dirname($wig_files[0]);

    my $output_file = $outdirectory."\/PEAKS_max_width_".$max_peak_width.".fraction_max_peaks_".$fraction_value_for_adjacent_peaks.".min_single_pos_depth_".$min_single_position_depth.".strand_".$strand.".bed";
    open(OUT, ">$output_file") or die "couldn't open output file";
    my $peak_count = 0;
    for(my $i=0; $i < $second_run_peaks_length; $i++) {
        $peak_count++;
        my @split_i = split("\t", $second_run_peaks[$i]);
        my @add_j_positions_loop_value = ();
        for(my $j=0; $j < $all_positions_length; $j++) {
            my @split_j = split("\t", $all_positions[$j]);
            
            if($split_i[0] eq $split_j[0]) {
                if($split_j[1] >= $split_i[1]-$max_peak_width and $split_j[1] <= $split_i[1]+$max_peak_width) {
                    if(abs($split_j[4]) >= abs($split_i[4])*$fraction_value_for_adjacent_peaks) {
                        push(@add_j_positions_loop_value, $j);
                    }
                }
                elsif($split_j[1] > $split_i[1]+$max_peak_width) {
                    last;
                }
            }
        }

        my $add_j_positions_loop_value_length = @add_j_positions_loop_value;
        my $loop_start = $add_j_positions_loop_value[0];
        my $loop_end = $add_j_positions_loop_value[-1];

        my $chr;
        my $start_coord;
        my $end_coord;
        my $total_counts = 0;

        my $max_amp = 0;
        my $max_amp_chr;
        my $max_amp_coord;

        for(my $j=$loop_start; $j <= $loop_end; $j++) {
            my @split_j = split("\t", $all_positions[$j]);
            if($j == $loop_start) {
                $chr = $split_j[0];
                $start_coord = $split_j[1]-1;
                $end_coord = $split_j[2];
                $total_counts = abs($split_j[4])+$total_counts;
                $strand = $split_j[5];

                $max_amp = abs($split_j[4]);
                $max_amp_chr = $chr;
                $max_amp_coord = $split_j[1];
            }
            elsif($j == $loop_end) {
                $end_coord = $split_j[2];
                $total_counts = abs($split_j[4])+$total_counts;
                if(abs($split_j[4]) > $max_amp) {
                    $max_amp = abs($split_j[4]);
                    $max_amp_coord = $split_j[1];
                }
            }
            else {
                $total_counts = abs($split_j[4])+$total_counts;
                if(abs($split_j[4]) > $max_amp) {
                    $max_amp = abs($split_j[4]);
                    $max_amp_coord = $split_j[1];
                }
            }
        }
        print OUT $chr, "\t", $start_coord, "\t", $end_coord, "\t",  $chr, "\;", $max_amp_coord-1, "\;", $max_amp_coord, "\;", "peak_", $peak_count, "_strand_", $strand, "\;", $max_amp, "\;", $strand, "\t", $total_counts, "\t", $strand, "\n";
    }
    close(OUT);
}
options;
qc;

merge_wigs;
determine_positions_with_min_single_position_depth;
filter_first_run_peaks;
identify_clusters;
