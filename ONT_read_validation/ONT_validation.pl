#!/usr/bin/perl

use warnings;
use strict;  
use File::Basename;

my $base_dirpath = dirname(__FILE__);
my $CAGE_cons_file;
my $min_CAGE_depth;
my $max_CAGE_distance;
my $threeP_file;
my $min_3P_depth;
my $max_3P_distance;
my $min_SJ_reads;
my $short_read_SJ_tab_file;
my $ONT_file;
my $genome_fasta = $base_dirpath."/accessory_files/EBV_fa/chrEBV_Akata_inverted_2.fa";
my $known_EBV_atg_file = $base_dirpath."/accessory_files/known_EBV_ATG/chrEBV_Akata_inverted_refined_genes_plus_features_annotation_cleaned_Eric_mod_revised.bed.known_EBV_ORF_start_sites.bed";

my @CAGE_peaks = ();
my $CAGE_peaks_length = 0;

my @threeP_peaks = ();
my $threeP_peaks_length = 0;

my @CAGE_validated = ();
my $CAGE_validated_length = 0;

my @CAGE_and_3P_validated = ();
my $CAGE_and_3P_validated_length = 0;

my @CAGE_and_3P_validated_AAUAAA = ();
my $CAGE_and_3P_validated_AAUAAA_length = 0;

my @genome_fa = ();
my $genome_fa_length = 0;

my @splice_junc = ();
my $splice_junc_length = 0;

my @CAGE_and_3P_validated_AAUAAA_SJ_val = ();
my $CAGE_and_3P_validated_AAUAAA_SJ_val_length = 0;

my @pos_EBV_atg_array = ();
my @neg_EBV_atg_array = ();

my $basename;
my $out_directory;
my $summary_file;
my $out_file;
my $negative_out_file;

sub program_info {
    print "\n\tOVERVIEW: ONT_validation.pl is designed to validate ONT reads using 1) splice junction data from short read alignments, 2) CAGE 5' end peaks, 3) ONT read 3' end peaks, and an input long read file in bed12 format. 
    
    \n\n\tInputs:\n\t\t-d\tIndicate whether dependency on BALF2, vPIC, or OriLyt is desired.\n\t\t-mcde\tMinimum CAGE peak depth\n\t\t-mcdi\tMaximum distance between CAGE peak and start of ONT read\n\t\t-3Pp\t3' peak file (bed6 format)\n\t\t-3Pde\tMinimum 3' peak read depth\n\t\t-3Pdi\tMaximum 3' distance between ONT end and 3' peak\n\t\t-minSJ\tMinimum number of short read splice junctions to validate spliced ONT reads\n\t\t-SJt\tSplice junction tab file from STAR alignment of short reads\n\t\t-ONT\tInput long read bed12 file
    
    \n\n\tOption:\n\t\t\-h help
    
    \n\n\tUsage: perl ONT_validation.pl [INPUTS] -d <BALF2 or vPIC or OriLYt> -mcde <Minimum CAGE peak depth> -mcdi <Maximum distance between CAGE peak and start of ONT read> -3Pp <3' peak file> -3Pde <Minimum read depth of 3' peak clusters> -3Pdi <Maximum 3' distance between ONT end and 3' peaks> -minSJ <Minimum number of validating short read splice junctions detected> -SJt <Splice junction tab file from STAR alignment of short reads> -ONT <Input long read bed12 file>
    
    \n\n\tExample: perl ONT_validation.pl -d BALF2 -mcde 1000 -mcdi 2 -3Pp 3Ppeakfile.bed -3Pde 400 -3Pdi 10 -minSJ 200 -SJt SJfile.tab -ONT ONTfile.bed\n\n";
    exit;
}

sub options {
    
    if (scalar @ARGV == 0) {
        program_info;
        exit;
    }    
    for (my $i=0; $i < scalar @ARGV; $i++) {
        if ($ARGV[$i] eq "\-d") {
            if($ARGV[$i+1] =~ m/BALF2|balf2/g) {
                $CAGE_cons_file = $base_dirpath."/accessory_files/CAGE_peaks/Akata_BCR_plus_Mutu_Zta_combined_CAGE_PEAKS_max_width_8.fraction_max_peaks_0.2.min_single_pos_depth_1000.strand_-_and_+_EBV_only.compiled.BALF2.bed";
            }
            elsif($ARGV[$i+1] =~ /vPIC|vpic|VPIC/g) {
                $CAGE_cons_file = $base_dirpath."/accessory_files/CAGE_peaks/Akata_BCR_plus_Mutu_Zta_combined_CAGE_PEAKS_max_width_8.fraction_max_peaks_0.2.min_single_pos_depth_1000.strand_-_and_+_EBV_only.compiled.vPIC.bed";
            }
            elsif($ARGV[$i+1] =~ /OriLyt|ORILYT|orilyt/g) {
                $CAGE_cons_file = $base_dirpath."/accessory_files/CAGE_peaks/Akata_BCR_plus_Mutu_Zta_combined_CAGE_PEAKS_max_width_8.fraction_max_peaks_0.2.min_single_pos_depth_1000.strand_-_and_+_EBV_only.compiled.orilyt.bed";
            }
        }
        elsif ($ARGV[$i] eq "\-mcde") {
            $min_CAGE_depth = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-mcdi") {
            $max_CAGE_distance = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-3Pp") {
            $threeP_file = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-3Pde") {
            $min_3P_depth = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-3Pdi") {
            $max_3P_distance = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-minSJ") {
            $min_SJ_reads = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-SJt") {
            $short_read_SJ_tab_file = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-ONT") {
            $ONT_file = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-h") {
            program_info;
            exit;
        }
    }
}

sub qc {
    if (not defined($CAGE_cons_file)) {
        print "\nBALF2, vPIC, or OriLyt value not defined!\n";
        program_info;
        exit;
    }
    elsif (not defined($min_CAGE_depth)) {
        print "\nMinimum CAGE depth not defined!\n";
        program_info;
        exit;
    }
    elsif (not defined($max_CAGE_distance)) {
        print "\nMaximum CAGE distance not defined!\n";
        program_info;
        exit;
    }
    elsif (not defined($threeP_file)) {
        print "\n3 prime peak file not defined!\n";
        program_info;
        exit;
    }
    elsif (not defined($min_3P_depth)) {
        print "\nMinimum 3P peak cluster coverage depth not defined!\n";
        program_info;
        exit;
    }
    elsif (not defined($max_3P_distance)) {
        print "\nMaximum distance between 3' ONT read and 3' cluster not defined!\n";
        program_info;
        exit;
    }
    elsif (not defined($min_SJ_reads)) {
        print "\nMinimum number of validating splice junction reads not defined!\n";
        program_info;
        exit;
    }
    elsif (not defined($short_read_SJ_tab_file)) {
        print "\nShort read splice junction (SJ) tab file not defined!\n";
        program_info;
        exit;
    }
    elsif (not defined($ONT_file)) {
        print "\nInput long read bed file not defined!\n";
        program_info;
        exit;
    }
}

sub setup {
    $basename = basename($ONT_file);
    $basename =~ s/\.bed$//g;
    $out_directory = dirname($ONT_file)."\/OUT_Directory_".$basename;
    `mkdir $out_directory`;
    $summary_file = $out_directory."\/".$basename."_summary_analysis_stats.txt";
    $out_file = $out_directory."\/VALIDATED_".$basename."\.bed";

    open(OUT2, ">$summary_file") or die "couldn't open summary file";
    print OUT2 "CAGE consensus file: ", basename($CAGE_cons_file), "\nMin CAGE depth: ", $min_CAGE_depth, "\nMax distance from ONT start to CAGE peak: ", $max_CAGE_distance, "\nThree prime end peak file: ", basename($threeP_file), "\nMinimum 3Prime peak depth: ", $min_3P_depth, "\nMax distance from ONT end to 3P peak: ", $max_3P_distance, "\nMinimum splice junction reads: ", $min_SJ_reads, "\nGenome fasta file: ", basename($genome_fasta), "\nShort read SJ tab file: ", basename($short_read_SJ_tab_file), "\nKnown EBV start site file: ", basename($known_EBV_atg_file), "\nONT input file: ", basename($ONT_file), "\n";

    $negative_out_file = $out_directory."\/".$basename."_negative_splice_junction_reads.bed";
    open(OUT3, ">$negative_out_file") or die "couldn't open summary file";
}

sub fa_unwrapper {
    print "\nProcessing genome fasta file...\n\n"; #Unwrapping and putting into array
    open(INF, "<$genome_fasta") or die "couldn't open genome fasta file";
    my @line_array = ();
    while(my $line = <INF>) {
        chomp($line);
        if ($. == 1) {
            my @split_line = split(" ", $line);
            my $chromosome = $split_line[0];
            $chromosome =~ s/\>//g;
            push(@genome_fa, $chromosome);
        }
        elsif ($line =~ m/^\>chr/) {
            push(@genome_fa, join("", @line_array));
            my @split_line = split(" ", $line);
            my $chromosome = $split_line[0];
            $chromosome =~ s/\>//g;
            push(@genome_fa, $chromosome);
            @line_array = ();
        }
        elsif (eof(INF)) {
            $line =~ tr/a-z/A-Z/;
            push(@line_array, $line);
            push(@genome_fa, join("", @line_array));
        }
        else { 
            $line =~ tr/a-z/A-Z/;
            push(@line_array, $line);
        }
    }
    close(INF);
    $genome_fa_length = @genome_fa;
}

sub process_CAGE_cons_file {
    print "Processing CAGE consensus file...\n\n";
    open (INF, "<$CAGE_cons_file") or die "couldn't open input file";

    while (my $line = <INF>) {
	    chomp($line);
        $line =~ s/\"//g;
        my @split_line = split("\t", $line);
        if($split_line[4] >= $min_CAGE_depth) {
            push(@CAGE_peaks, $line);
        }
    }
    close(INF);
    $CAGE_peaks_length = @CAGE_peaks;
}

sub process_3P_file {
    print "Processing 3 prime end file...\n\n";
    open (INF, "<$threeP_file") or die "couldn't open input file";

    while (my $line = <INF>) {
	    chomp($line);
        my @split_line = split("\t", $line);
        if(abs($split_line[4]) >= $min_3P_depth) {
            push(@threeP_peaks, $line);
        }
    }
    close(INF);
    $threeP_peaks_length = @threeP_peaks;
}

sub process_splicing_file {
    print "Processing splicing tab file...\n\n"; #
    open(INF, "<$short_read_SJ_tab_file") or die "couldn't open genome fasta file";
    while(my $line = <INF>) {
        chomp($line);
        my @split_line = split("\t", $line);
        if($split_line[6] > $min_SJ_reads){
            if ($split_line[3] == 1) {
                my $bed_line = join("\t", @split_line[0..2])."\tjunc\t".$split_line[6]."\t\+\n";
                push(@splice_junc, $bed_line);
            }
            elsif ($split_line[3] == 2) {
                my $bed_line = join("\t", @split_line[0..2])."\tjunc\t".$split_line[6]."\t\-\n";
                push(@splice_junc, $bed_line);
            }
        }
    }
    close(INF);
    $splice_junc_length = @splice_junc;
}

sub ONT_CAGE_overlap {
    print "Determining 5' end CAGE overlap with ONT reads (will take a few minutes)...\n\n";
    open (INF, "<$ONT_file") or die "couldn't open input file";
    my $count = 0;
    my $hit_count = 0;
    while (my $line = <INF>) {
        $count++;
	    chomp($line);
        my @split_line = split("\t", $line);
        
        for (my $i=0; $i < $CAGE_peaks_length; $i++) {
            my @split_CAGE_line = split("\t", $CAGE_peaks[$i]);
            my @split_CAGE_info = split("\;", $split_CAGE_line[3]);
            if($split_CAGE_line[0] eq $split_line[0]) {
                if($split_CAGE_line[5] eq "+" and $split_line[5] eq "+") {
                    if($split_line[1] >= $split_CAGE_line[1]-$max_CAGE_distance+1 and  $split_line[1] <= $split_CAGE_line[2]+$max_CAGE_distance+1) {
                        my $line_plus_CAGE_value = join("\t", @split_line[0..3])."\t".$split_line[4]."_CAGE_val-".abs($split_CAGE_line[4])."-coord-".$split_CAGE_info[1]."\t".join("\t", @split_line[5..7])."\t".$split_CAGE_line[8]."\t".join("\t", @split_line[9..11]);
                        push(@CAGE_validated, $line_plus_CAGE_value);
                        $hit_count++;
                    }
                }
                elsif($split_CAGE_line[5] eq "-" and $split_line[5] eq "-") {
                    if($split_line[2] >= $split_CAGE_line[1]-$max_CAGE_distance-1 and  $split_line[2] <= $split_CAGE_line[2]+$max_CAGE_distance-1) {
                        my $line_plus_CAGE_value = join("\t", @split_line[0..3])."\t".$split_line[4]."_CAGE_val-".abs($split_CAGE_line[4])."-coord-".$split_CAGE_info[2]."\t".join("\t", @split_line[5..7])."\t".$split_CAGE_line[8]."\t".join("\t", @split_line[9..11]);
                        push(@CAGE_validated, $line_plus_CAGE_value);
                        $hit_count++;
                    }
                }
            }
        }
    }
    $CAGE_validated_length = @CAGE_validated;
    close(INF);
    print OUT2 "\nFraction of ONT reads with CAGE peak validation = ", $hit_count/$count, "\n";
}

sub ONT_3_prime_overlap {
    print "Determining 3' end overlap with CAGE validated ONT reads...\n\n";
    my $count = 0;
    my $hit_count = 0;
    for(my $j=0; $j < $CAGE_validated_length; $j++) {
        my @split_CAGE_validated = split("\t", $CAGE_validated[$j]);
        $count++;

        for (my $i=0; $i < $threeP_peaks_length; $i++) {
            my @split_threeP_peaks_line = split("\t", $threeP_peaks[$i]);
            my @split_threeP_peaks_info = split("\;", $split_threeP_peaks_line[3]);
            if($split_CAGE_validated[0] eq $split_threeP_peaks_line[0]) {
                if($split_CAGE_validated[5] eq "+" and $split_threeP_peaks_line[5] eq "+") {
                    if($split_CAGE_validated[2] >= $split_threeP_peaks_line[1]-$max_3P_distance and $split_CAGE_validated[2] <= $split_threeP_peaks_line[2]+$max_3P_distance) {
                        my $push_CAGE_validated = join("\t", @split_CAGE_validated[0..4])."-3P_val-".abs($split_threeP_peaks_line[4])."-coord-".$split_threeP_peaks_info[2]."\t".join("\t", @split_CAGE_validated[5..11]);
                        push(@CAGE_and_3P_validated, $push_CAGE_validated);
                        $hit_count++;
                    }
                }
                elsif($split_CAGE_validated[5] eq "-" and $split_threeP_peaks_line[5] eq "-") {
                    if($split_CAGE_validated[1] >= $split_threeP_peaks_line[1]-$max_3P_distance and  $split_CAGE_validated[1] <= $split_threeP_peaks_line[2]+$max_3P_distance) {
                        my $push_CAGE_validated = join("\t", @split_CAGE_validated[0..4])."-3P_val-".abs($split_threeP_peaks_line[4])."-coord-".$split_threeP_peaks_info[1]."\t".join("\t", @split_CAGE_validated[5..11]);
                        push(@CAGE_and_3P_validated, $push_CAGE_validated);
                        $hit_count++;
                    }
                }
            }
        }
    }
    $CAGE_and_3P_validated_length = @CAGE_and_3P_validated;
    print OUT2 "\nFraction of CAGE validated ONT reads with 3P validation = ", $hit_count/$count, "\n";
}

sub ONT_AAUAAA_motif {
    print "Assessing presence of AAUAAA motifs in 5' and 3' validated ONT reads...\n\n";
    my $count = 0;
    my $hit_count = 0;
    my $hit_AUUAAA_count = 0;
    my $hit_AAUACA_count = 0;
    my $hit_GAUAAA_count = 0;
    for (my $i=0; $i < $CAGE_and_3P_validated_length; $i++) {
        $count++;
        my @split_transcript = split("\t", $CAGE_and_3P_validated[$i]);

        my $threeP_seq;
        for (my $j=0; $j < $genome_fa_length; $j++) {
            if($split_transcript[0] eq $genome_fa[$j]) {
                my $chr_seq = $genome_fa[$j+1];
                if($split_transcript[5] eq "+") {
                    $threeP_seq = substr($chr_seq, $split_transcript[2]-35, 30);
                    if($threeP_seq =~ m/AATAAA/) {
                        my $line = join("\t", @split_transcript[0..11])."\tAAUAAA"; 
                        push(@CAGE_and_3P_validated_AAUAAA, $line); 
                        $hit_count++;
                    }
                    elsif($threeP_seq =~ m/ATTAAA/) {
                        my $line = join("\t", @split_transcript[0..11])."\tAUUAAA"; 
                        push(@CAGE_and_3P_validated_AAUAAA, $line); 
                        $hit_AUUAAA_count++;
                    }
                    elsif($threeP_seq =~ m/AATACA/) {
                        my $line = join("\t", @split_transcript[0..11])."\tAAUACA"; 
                        push(@CAGE_and_3P_validated_AAUAAA, $line); 
                        $hit_AAUACA_count++;
                    }
                    elsif($threeP_seq =~ m/GATAAA/) {
                        my $line = join("\t", @split_transcript[0..11])."\tGAUAAA"; 
                        push(@CAGE_and_3P_validated_AAUAAA, $line); 
                        $hit_GAUAAA_count++;
                    }
                    else {
                        my $line = join("\t", @split_transcript[0..11])."\tNone";
                    }
                }
                elsif($split_transcript[5] eq "-") {
                    $threeP_seq = substr($chr_seq, $split_transcript[1]+5, 30);
                    if($threeP_seq =~ m/TTTATT/) {
                        my $line = join("\t", @split_transcript[0..11])."\tAAUAAA"; 
                        push(@CAGE_and_3P_validated_AAUAAA, $line); 
                        $hit_count++;
                    }
                    elsif($threeP_seq =~ m/TTTAAT/) {
                        my $line = join("\t", @split_transcript[0..11])."\tAUUAAA"; 
                        push(@CAGE_and_3P_validated_AAUAAA, $line); 
                        $hit_AUUAAA_count++;
                    }
                    elsif($threeP_seq =~ m/TGTATT/) {
                        my $line = join("\t", @split_transcript[0..11])."\tAAUACA"; 
                        push(@CAGE_and_3P_validated_AAUAAA, $line); 
                        $hit_AAUACA_count++;
                    }
                    elsif($threeP_seq =~ m/TTTATC/) {
                        my $line = join("\t", @split_transcript[0..11])."\tGAUAAA"; 
                        push(@CAGE_and_3P_validated_AAUAAA, $line); 
                        $hit_GAUAAA_count++;
                    }
                    else {
                        my $line = join("\t", @split_transcript[0..11])."\tNone";
                    }
                }
            }

        }
    }
    $CAGE_and_3P_validated_AAUAAA_length = @CAGE_and_3P_validated_AAUAAA;
    print OUT2 "\nFraction of 3P and CAGE validated ONT reads with upstream AAUAAA = ", $hit_count/$count, "\n";
    print OUT2 "Fraction of 3P and CAGE validated ONT reads with upstream AUUAAA = ", $hit_AUUAAA_count/$count, "\n";
    print OUT2 "Fraction of 3P and CAGE validated ONT reads with upstream AAUACA = ", $hit_AAUACA_count/$count, "\n";
    print OUT2 "Fraction of 3P and CAGE validated ONT reads with upstream GAUAAA = ", $hit_GAUAAA_count/$count, "\n";
}

sub validate_SJs {
    print "Validating SJ in spliced ONT reads using short read data (STAR .tab file)...\n\n";
    for(my $j=0; $j < $CAGE_and_3P_validated_AAUAAA_length; $j++) {
        my @split_CAGE_and_3P_validated_AAUAAA_line = split("\t", $CAGE_and_3P_validated_AAUAAA[$j]);

        if($split_CAGE_and_3P_validated_AAUAAA_line[9] == 1) {
            my $line = $CAGE_and_3P_validated_AAUAAA[$j]."\tno splice";
            push(@CAGE_and_3P_validated_AAUAAA_SJ_val, $line);
        }
        else {
            my @hit_array = ();
            my @size_array = split("\,", $split_CAGE_and_3P_validated_AAUAAA_line[10]);
            my @start_array = split("\,", $split_CAGE_and_3P_validated_AAUAAA_line[11]);
            
            for (my $e=0; $e < $split_CAGE_and_3P_validated_AAUAAA_line[9]-1; $e++) {
                my $intron_start = $split_CAGE_and_3P_validated_AAUAAA_line[1]+$start_array[$e]+$size_array[$e]+1;
                my $intron_end = $split_CAGE_and_3P_validated_AAUAAA_line[1]+$start_array[$e+1];
                
                my $intron_hits = 0;
                for (my $i=0; $i < $splice_junc_length; $i++) {
                    my @splice_junc_line = split("\t", $splice_junc[$i]);
                    if($split_CAGE_and_3P_validated_AAUAAA_line[0] eq $splice_junc_line[0]) {
                        if($intron_start == $splice_junc_line[1] and $intron_end == $splice_junc_line[2]) {
                            $intron_hits++;
                        }
                    }
                }
                if($intron_hits == 0) {
                    push(@hit_array, "false");
                }
                else {
                    push(@hit_array, "true");
                }
            }
            my $false = 0;
            for(my $e=0; $e < scalar(@hit_array); $e++) {
                if($hit_array[$e] eq "false") {
                    $false++;
                }
            }

            if($false == 0) {
                my $line = $CAGE_and_3P_validated_AAUAAA[$j]."\t".join("\,", @hit_array);
                push(@CAGE_and_3P_validated_AAUAAA_SJ_val, $line);
            }
            else {
                my $line = $CAGE_and_3P_validated_AAUAAA[$j]."\t".join("\,", @hit_array);
                print OUT3 $line, "\n";
            }
        }
    }
    my $CAGE_and_3P_validated_AAUAAA_SJ_val_length = @CAGE_and_3P_validated_AAUAAA_SJ_val;
}

sub print_CAGE_3P_validated_ONT_reads {
    print "Printing CAGE validated, 3' end validated ONT reads...\n\n";

    open (OUT, ">$out_file") or die "couldn't open output file";
    print OUT join("\n", @CAGE_and_3P_validated_AAUAAA_SJ_val);
    close(OUT);
}

sub consolidate_reads {
    print "Consolidating reads and generating isoform annotation...\n\n";
    my $transcript_out_file = $out_directory."\/VALIDATED_transcripts_".$basename."\.bed";
    open (INF, "<$out_file") or die "couldn't open input file";
    open (OUT, ">$transcript_out_file") or die "couldn't open output file";
    while (my $line = <INF>) {
	    chomp($line);
        my @split_line = split("\t", $line);
        my @split_CAGE_3P_data = split("\-", $split_line[4]);
        my $offset1;
        my $offset2;
        my $exon_number = $split_line[9];
        my @split_exon_sizes = split("\,", $split_line[10]);
        my @revised_exon_sizes = ();
        my @split_exon_starts = split("\,", $split_line[11]);
        my @revised_exon_starts = ();
        if($split_line[5] eq "+") {
            $offset1 = $split_CAGE_3P_data[3] - $split_line[1];
            $offset2 = $split_line[2] - $split_CAGE_3P_data[7];
            for (my $i=0; $i<$exon_number; $i++) {
                if($i == 0) {
                    push(@revised_exon_starts, $split_exon_starts[$i]);
                    if($exon_number == 1) {
                        my $new_exon_size = $split_exon_sizes[$i] - ($offset1 + $offset2);
                        push(@revised_exon_sizes, $new_exon_size);
                    }
                    else {
                        my $new_exon_size = $split_exon_sizes[$i] - $offset1;
                        push(@revised_exon_sizes, $new_exon_size);
                    }
                }
                elsif($i == $exon_number-1) {
                    my $new_exon_size = $split_exon_sizes[$i] - $offset2;
                    push(@revised_exon_sizes, $new_exon_size);
                    my $new_exon_start = $split_exon_starts[$i] - $offset1;
                    push(@revised_exon_starts, $new_exon_start);
                }
                else {
                    push(@revised_exon_sizes, $split_exon_sizes[$i]);
                    my $new_exon_start = $split_exon_starts[$i] - $offset1;
                    push(@revised_exon_starts, $new_exon_start);
                }
            }
           
            print OUT $split_line[0], "\t", $split_CAGE_3P_data[3], "\t", $split_CAGE_3P_data[7], "\t", join("\t", @split_line[3..5]), "\t",  $split_CAGE_3P_data[3], "\t", $split_CAGE_3P_data[7], "\t", join("\t", @split_line[8..9]), "\t", join("\,", @revised_exon_sizes), "\t", join("\,", @revised_exon_starts), "\t", join("\t", @split_line[12..13]), "\t", $split_CAGE_3P_data[1], "\t", $split_CAGE_3P_data[5], "\n";
        }
        elsif($split_line[5] eq "-") {
            $offset1 = $split_CAGE_3P_data[7] - $split_line[1];
            $offset2 = $split_line[2] - $split_CAGE_3P_data[3];
            for (my $i=0; $i<$exon_number; $i++) {
                if($i == 0) {
                    push(@revised_exon_starts, $split_exon_starts[$i]);
                    if($exon_number == 1) {
                        my $new_exon_size = $split_exon_sizes[$i] - ($offset1 + $offset2);
                        push(@revised_exon_sizes, $new_exon_size);
                    }
                    else {
                        my $new_exon_size = $split_exon_sizes[$i] - $offset1;
                        push(@revised_exon_sizes, $new_exon_size);
                    }
                }
                elsif($i == $exon_number-1) {
                    my $new_exon_size = $split_exon_sizes[$i] - $offset2;
                    push(@revised_exon_sizes, $new_exon_size);
                    my $new_exon_start = $split_exon_starts[$i] - $offset1;
                    push(@revised_exon_starts, $new_exon_start);
                }
                else {
                    push(@revised_exon_sizes, $split_exon_sizes[$i]);
                    my $new_exon_start = $split_exon_starts[$i] - $offset1;
                    push(@revised_exon_starts, $new_exon_start);
                }
            }
            print OUT $split_line[0], "\t", $split_CAGE_3P_data[7], "\t", $split_CAGE_3P_data[3], "\t", join("\t", @split_line[3..5]), "\t", $split_CAGE_3P_data[7], "\t", $split_CAGE_3P_data[3], "\t", join("\t", @split_line[8..9]), "\t", join("\,", @revised_exon_sizes), "\t", join("\,", @revised_exon_starts), "\t", join("\t", @split_line[12..13]), "\t", $split_CAGE_3P_data[1], "\t", $split_CAGE_3P_data[5], "\n";
        }
    }
    close(INF);
    close(OUT);

    my $sorted_transcript_out_file = $out_directory."\/VALIDATED_transcripts_".$basename.".sorted\.bed";
    `sort -V -k 1,1 -k 2,2n -k 3,3n -k 11,11 -k 12,12 $transcript_out_file > $sorted_transcript_out_file`;

    my $single_entry_transcript_out_file = $out_directory."\/VALIDATED_transcripts_".$basename.".single_entry\.bed";
    open (INF, "<$sorted_transcript_out_file") or die "couldn't open input file";
    open (OUT, ">$single_entry_transcript_out_file") or die "couldn't open output file";
    my $count = 0;
    my $prev_line;
    my $prev_chr = "null";
    my $prev_start = "0";
    my $prev_end = "0";
    my $prev_sizes = "null";
    my $prev_starts = "null";
    while (my $line = <INF>) {
	    chomp($line);
        my @split_line = split("\t", $line);

        if($. == 1) {
            $prev_chr = $split_line[0];
            $prev_start = $split_line[1];
            $prev_end = $split_line[2];
            $prev_sizes = $split_line[10];
            $prev_starts = $split_line[11];
            $prev_line = $line;
            $count++;
        }

        elsif($split_line[0] eq $prev_chr and $split_line[1] == $prev_start and $split_line[2] == $prev_end and $prev_sizes eq $split_line[10] and $split_line[11] eq $prev_starts) {
            $count++;
            if(eof(INF)) {
                my ($thick_start, $thick_end) = find_reading_frames($prev_line);
                my @split_prev_line = split("\t", $prev_line);
                print OUT join("\t", @split_prev_line[0..5]), "\t", $thick_start, "\t", $thick_end, "\t", join("\t", @split_prev_line[8..15]), "\t", $count, "\n";
            }
        }
        elsif(eof(INF)) {
            my ($thick_start, $thick_end) = find_reading_frames($prev_line);
            my @split_prev_line = split("\t", $prev_line);
            print OUT join("\t", @split_prev_line[0..5]), "\t", $thick_start, "\t", $thick_end, "\t", join("\t", @split_prev_line[8..15]), "\t", $count, "\t1\n";
        }
        else {
            my ($thick_start, $thick_end) = find_reading_frames($prev_line);
            my @split_prev_line = split("\t", $prev_line);
            print OUT join("\t", @split_prev_line[0..5]), "\t", $thick_start, "\t", $thick_end, "\t", join("\t", @split_prev_line[8..15]), "\t", $count, "\n";
            $prev_chr = $split_line[0];
            $prev_start = $split_line[1];
            $prev_end = $split_line[2];
            $prev_sizes = $split_line[10];
            $prev_starts = $split_line[11];
            $prev_line = $line;
            $count = 1;
        }
    }
    close(INF);
    close(OUT);
}

sub process_known_atg_file {
    print "Processing CAGE consensus file...\n\n";
    open (INF, "<$known_EBV_atg_file") or die "couldn't open input file";

    while (my $line = <INF>) {
	    chomp($line);
        my @split_line = split("\t", $line);
        if($split_line[1] eq "+") {
            push(@pos_EBV_atg_array, $split_line[0]);
        }
        elsif($split_line[1] eq "-") {
            push(@neg_EBV_atg_array, $split_line[0]);
        }
    }
    close(INF);
}

sub find_reading_frames {
    my $line = $_[0];
    my @split_line = split("\t", $line);
    my @split_sizes = split("\,", $split_line[10]);
    my @split_starts = split("\,", $split_line[11]);
    my $block_number = $split_line[9];
    my $chr = $split_line[0];
    my @ATG_positions = ();
    my @TERM_positions = ();

    my $thick_start = $split_line[1];
    my $thick_end = $split_line[2];

    my $known_atg_coord = "null";
    my $known_atg_transcript_coord = "null";
    for (my $x = 0; $x < scalar(@genome_fa); $x++) {
        if($chr eq $genome_fa[$x]) {
            my $chr_seq = $genome_fa[$x+1];
            my $transcript_seq = "";
            my @positions_array = ();
            
            my $transcript_position_count = 0;
            my $for_loop_count = 0;
            for (my $b = 0; $b < $block_number; $b++) {
                my $block_start = $split_line[1]+$split_starts[$b];
                my $block_end = $split_line[1]+$split_starts[$b]+$split_sizes[$b];

                if($split_line[5] eq "+") {
                    
                    for(my $atg=0; $atg < scalar(@pos_EBV_atg_array); $atg++) {
                        if($for_loop_count == 0 and $pos_EBV_atg_array[$atg] >= $block_start and $pos_EBV_atg_array[$atg] <= $block_end) {
                            $known_atg_coord = $pos_EBV_atg_array[$atg];
                            $for_loop_count++;
                        }
                    }
                }
                elsif($split_line[5] eq "-") {
                    for(my $atg=0; $atg < scalar(@neg_EBV_atg_array); $atg++) {
                        if($neg_EBV_atg_array[$atg] >= $block_start and $neg_EBV_atg_array[$atg] <= $block_end) {
                            $known_atg_coord = $neg_EBV_atg_array[$atg];
                        }
                    }
                }

                my $block_seq = substr($chr_seq, $block_start, $split_sizes[$b]);
                $transcript_seq = $transcript_seq.$block_seq;

                for (my $pa = $block_start; $pa < $block_end; $pa++) {
                    $transcript_position_count++;
                    push(@positions_array, $pa);
                    if($known_atg_coord ne "null") {
                        if($pa == $known_atg_coord) {
                            if ($split_line[5] eq "+") {
                                $known_atg_transcript_coord = ($transcript_position_count-1)."\;".($transcript_position_count+1);
                            }
                            elsif ($split_line[5] eq "-") {
                                $known_atg_transcript_coord = ($transcript_position_count-4)."\;".($transcript_position_count-2);
                            }
                        }
                    }
                }
            }
                
            if($split_line[5] eq "+") {
                if($known_atg_coord eq "null") {
                    @ATG_positions = match_all_positions("CATG|ATGG", $transcript_seq, "plus");
                }
                else {
                    push(@ATG_positions, $known_atg_transcript_coord);
                }
                
                @TERM_positions = match_all_positions("TAA|TGA|TAG", $transcript_seq, "null");

                my $hit = 0;
                for (my $i=0; $i < scalar(@ATG_positions); $i++) {
                    my @split_start = split("\;", $ATG_positions[$i]);
                    my $ATG_position = $split_start[0];
                    for (my $j=0; $j < scalar(@TERM_positions); $j++) {
                        my @split_TERM = split("\;", $TERM_positions[$j]);
                        my $end_position = $split_TERM[0];
                        my $length = $end_position - $ATG_position;

                        if ($length % 3 == 0) {
                            if($length >= 300 || ($known_atg_coord ne "null" and $length > 0)) {
                                $hit++;
                                $thick_start = $positions_array[$ATG_position];
                                $thick_end = $positions_array[$end_position]+3;
                                return($thick_start, $thick_end);
                                $i = scalar(@ATG_positions);
                            }
                            if($length > 0) {
                                last;
                            }
                        }
                    }
                }
                if($hit == 0) {
                    return($thick_start, $thick_start);
                }
            }
            elsif($split_line[5] eq "-") {
                if($known_atg_coord eq "null") {
                    @ATG_positions = match_all_positions("CCAT|CATG", $transcript_seq, "minus");
                }
                else {
                    push(@ATG_positions, $known_atg_transcript_coord);
                }
                @TERM_positions = match_all_positions("TTA|TCA|CTA", $transcript_seq, "null");

                my $hit = 0;
                for (my $i=scalar(@ATG_positions)-1; $i >= 0; $i--) {
                    my @split_start = split("\;", $ATG_positions[$i]);
                    my $ATG_position = $split_start[1]+1;
                    for (my $j=scalar(@TERM_positions)-1; $j >= 0; $j--) {
                        my @split_TERM = split("\;", $TERM_positions[$j]);
                        my $end_position = $split_TERM[1]+1;
                        my $length = $ATG_position - $end_position;
                        if ($length % 3 == 0) {
                            if($length >= 300 || ($known_atg_coord ne "null" and $length > 0)) {
                                $hit++;
                                $thick_start = $positions_array[$end_position]-3;
                                $thick_end = $positions_array[$ATG_position];
                                return($thick_start, $thick_end);
                                $i = -1;
                            }
                            if($length > 0) {
                                last;
                            }
                        }
                    }
                }
                if($hit == 0) {
                    return($thick_start, $thick_start);
                }
            }
            
        }
        last;
    }
}

sub match_all_positions {
    my ($term, $sequence, $strand) = @_;
    my @output;
    while ($sequence =~ m/$term/g) {
        if($strand eq "plus") {
            if(substr($sequence, $-[0], ($+[0]) - ($-[0])) eq "CATG") {
                push (@output, join("\;", $-[0]+1, $+[0]-1));
            }
            elsif(substr($sequence, $-[0], ($+[0]) - ($-[0])) eq "ATGG") {
                push (@output, join("\;", $-[0], $+[0]-2));
            }
        }
        elsif($strand eq "minus") {
            if(substr($sequence, $-[0], ($+[0]) - ($-[0])) eq "CCAT") {
                push (@output, join("\;", $-[0]+1, $+[0]-1));

            }
            elsif(substr($sequence, $-[0], ($+[0]) - ($-[0])) eq "CATG") {
                push (@output, join("\;", $-[0], $+[0]-2));
            }
        }
        else {
            push (@output, join("\;", $-[0], $+[0]-1));
        }
    }
    return @output;
}

sub deduplicate_from_validated_single_entry_file {
    print "Deduplicating validated single entry file...\n\n";

    my $file = $out_directory."\/VALIDATED_transcripts_".$basename.".single_entry\.bed";
    my $sorted_file = $out_directory."\/VALIDATED_transcripts_".$basename.".single_entry_sorted";
    my $output_file = $out_directory."\/VALIDATED_transcripts_".$basename.".single_entry_deduplicated.bed";

    `sort -k 4,4 $file > $sorted_file`;

    open(INF, "<$sorted_file") or die "couldn't open genome fasta file";
    open(OUT, ">$output_file") or die "couldn't open summary file";
    my $prev_coord1 = 0;
    my $prev_coord2 = 0;
    my $prev_transcript_ID = "null";
    my $prev_CAGE_depth = 0;
    my $prev_line = "null";
    while(my $line = <INF>) {
        chomp($line);
        my @split_line = split("\t", $line);
        if($. == 1) {
            $prev_coord1 = $split_line[1];
            $prev_coord2 = $split_line[2];
            $prev_transcript_ID = $split_line[3];
            $prev_CAGE_depth = $split_line[14];
            $prev_line = $line;
        }
        elsif(eof) {
            if($prev_transcript_ID ne $split_line[3]) {
                print OUT $prev_line, "\n";
                print OUT $line, "\n";
            }
            else {
                if ($split_line[5] eq "+") {
                    if($split_line[2] >= $prev_coord2 and $split_line[14] >= $prev_CAGE_depth) {
                        print OUT $line, "\n";
                    }
                    else {
                        print OUT $prev_line, "\n";
                    }
                }
                elsif ($split_line[5] eq "-") {
                    if($split_line[1] <= $prev_coord2 and $split_line[14] >= $prev_CAGE_depth) {
                        print OUT $line, "\n";
                    }
                    else {
                        print OUT $prev_line, "\n";
                    }
                }
            }
        }
        elsif($prev_transcript_ID ne $split_line[3]) {
            print OUT $prev_line, "\n";
            $prev_coord1 = $split_line[1];
            $prev_coord2 = $split_line[2];
            $prev_transcript_ID = $split_line[3];
            $prev_CAGE_depth = $split_line[14];
            $prev_line = $line;
        }
        elsif ($split_line[5] eq "+") {
            if($split_line[2] >= $prev_coord2 and $split_line[14] >= $prev_CAGE_depth) {
                $prev_coord1 = $split_line[1];
                $prev_coord2 = $split_line[2];
                $prev_transcript_ID = $split_line[3];
                $prev_CAGE_depth = $split_line[14];
                $prev_line = $line;
            }
        }
        elsif ($split_line[5] eq "-") {
            if($split_line[1] <= $prev_coord1 and $split_line[14] >= $prev_CAGE_depth) {
                $prev_coord1 = $split_line[1];
                $prev_coord2 = $split_line[2];
                $prev_transcript_ID = $split_line[3];
                $prev_CAGE_depth = $split_line[14];
                $prev_line = $line;
            }
        }
    }
    close(INF);
    close(OUT);
    `rm $sorted_file`;
}

sub final_clean_merge_competing_dual_3p_ends {
    print "Merging competing dual 3p ends...\n\n";

    my $file = $out_directory."\/VALIDATED_transcripts_".$basename.".single_entry_deduplicated.bed";
    my $output_file = $file;
    $output_file =~ s/\.single_entry_deduplicated\.bed//g;
    $output_file =~ s/VALIDATED_//g;
    $output_file = $output_file."_confirmed.bed";

    my @array = ();
    open (INF, "<$file") or die "couldn't open input file";
    open (OUT, ">$output_file") or die "couldn't open output file";
    while (my $line = <INF>) {
	    chomp($line);
        push(@array, $line);
    }
    close(INF);

    my @new_array = ();
    my @skip_elements = ();
    for (my $i = 0; $i<scalar(@array); $i++) {
        my @split_i = split("\t", $array[$i]);
        my $hit = 0;
        my $ONT_count = $split_i[16];
        for (my $j = 0; $j<scalar(@array); $j++) {
            my @split_j = split("\t", $array[$j]);
            next if ($i==$j);
            if($split_i[5] eq "+" and $split_j[5] eq "+" and $split_i[1] == $split_j[1] and $split_i[2] > $split_j[2] and $split_i[2]-37 < $split_j[2] and $split_i[9] == $split_j[9]) {
                my @split_i_starts = split("\,", $split_i[11]);
                my @split_i_sizes = split("\,", $split_i[10]);
                my @split_j_starts = split("\,", $split_j[11]);
                my @split_j_sizes = split("\,", $split_j[10]);
            
                for (my $k=0; $k < $split_i[9]; $k++) {
                    if($split_i[9]-1 == $k) {
                        if($split_i_sizes[$k]-($split_i[2]-$split_j[2]) == $split_j_sizes[$k]) {
                            $hit++;
                        }
                    }
                    elsif($split_i_sizes[$k] == $split_j_sizes[$k]) {
                        $hit++;
                    }
                }
                if($hit == $split_i[9]) {
                    $ONT_count = $ONT_count + $split_j[16];
                    push(@skip_elements, $j);
                }
            }
            elsif($split_i[5] eq "-" and $split_j[5] eq "-" and $split_i[2] == $split_j[2] and $split_i[1] < $split_j[1] and $split_i[1]+37 > $split_j[1] and $split_i[9] == $split_j[9]) {
                my @split_i_starts = split("\,", $split_i[11]);
                my @split_i_sizes = split("\,", $split_i[10]);
                my @split_j_starts = split("\,", $split_j[11]);
                my @split_j_sizes = split("\,", $split_j[10]);
            
                for (my $k=0; $k < $split_i[9]; $k++) {
                    if($k == 0) {
                        if($split_i_sizes[$k]-($split_j[1]-$split_i[1]) == $split_j_sizes[$k]) {
                            $hit++;
                        }
                    }
                    elsif($split_i_sizes[$k] == $split_j_sizes[$k]) {
                        $hit++;
                    }
                }
                if($hit == $split_i[9]) {
                    $ONT_count = $ONT_count + $split_j[16];
                    push(@skip_elements, $j);
                }
            }
        }

        my $new_element = join("\t", @split_i[0..15])."\t".$ONT_count;
        push(@new_array, $array[$i]);

    }

    for(my $a=0; $a<scalar(@new_array); $a++) {
        my $hit = 0;
        for(my $b=0; $b<scalar(@skip_elements); $b++) {
            if($a == $skip_elements[$b]) {
                $hit++;
            }
        }
        if($hit == 0) {
            print OUT $new_array[$a], "\n";
        }
    }
    close(OUT);
    my $rm1 = $out_directory."\/VALIDATED\*";
    `rm $rm1`;
    my $rm2 = $out_directory."\/\*negative_splice_junction_reads\*";
    `rm $rm2`;
}

options;
qc;
setup;
fa_unwrapper;
process_CAGE_cons_file;
process_3P_file;
process_splicing_file;
process_known_atg_file;
ONT_CAGE_overlap;
ONT_3_prime_overlap;
ONT_AAUAAA_motif;
validate_SJs;
print_CAGE_3P_validated_ONT_reads;
close(OUT2);

close(OUT3);

consolidate_reads;
deduplicate_from_validated_single_entry_file;
final_clean_merge_competing_dual_3p_ends;
print "DONE!!!\n\n";