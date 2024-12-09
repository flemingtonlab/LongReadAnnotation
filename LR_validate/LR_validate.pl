#!/usr/bin/perl
use warnings;
use strict;  
use File::Basename;
#use DateTime;

my $base_dirpath = dirname(__FILE__);
my $CAGE_cons_file;
my $min_CAGE_depth;
my $max_CAGE_distance;
my $threeP_file;
my $threeP_file_pos;
my $threeP_file_neg;
my $threeP_file_pos_sorted;
my $threeP_file_neg_sorted;
my $min_3P_depth;
my $max_3P_distance;
my $min_SJ_reads;
my $short_read_SJ_tab_file;
my $ONT_file;
my $genome_fasta;
my $known_atg_file = "none";

my @CAGE_peaks_pos = ();
my @CAGE_peaks_neg = ();

my @threeP_peaks_pos = ();
my @threeP_peaks_neg = ();

my @ONT_reads_pos;
my @ONT_reads_neg;

my @CAGE_validated_pos = ();
my @CAGE_validated_neg = ();

my @CAGE_and_3P_validated = ();

my @CAGE_and_3P_validated_AAUAAA_pos = ();
my @CAGE_and_3P_validated_AAUAAA_neg = ();

my @genome_fa = ();

my @splice_junc_pos = ();
my @splice_junc_neg = ();

my @CAGE_and_3P_validated_AAUAAA_SJ_val_pos = ();
my @CAGE_and_3P_validated_AAUAAA_SJ_val_neg = ();

my @pos_atg_array = ();
my @neg_atg_array = ();

my $prev_reading_frame_chr = "null";

my $pos_atg_array_chr_start = 0;
my $neg_atg_array_chr_start = 0;

my $basename;
my $out_directory;
my $summary_file;
my $out_file;
my $negative_out_file;

sub program_info {
    print "\n\tOVERVIEW: LR_validate.pl is designed to validate mapped long read data using 1) splice junction\n\tdata from short read alignments, 2) CAGE 5' end peaks, 3) LR read 3' end peaks, and\n\t4) an input mapped long read file in bed12 format. 
    
    \n\n\tInputs:\n\t\t-5Pp\tCAGE peak file (.bed10 format CAGE peak output from peak_caller_from_wigs.pl).\n\t\t-mcde\tMinimum CAGE peak depth\n\t\t-mcdi\tMaximum distance between CAGE peak and start of LR read\n\t\t-3Pp\t3' peak file (bed6 format)\n\t\t-3Pde\tMinimum 3' peak read depth\n\t\t-3Pdi\tMaximum 3' distance between LR end and 3' peak\n\t\t-minSJ\tMinimum number of short read splice junctions to validate spliced LRs\n\t\t-SJt\tSplice junction tab file from STAR alignment of short reads\n\t\t-f\tgenome fasta file\n\t\t-LR\tInput long read bed12 file
    
    \n\n\tOptions:\n\t\t\-ATG\tknown ATG starts site file\n\t\t\-h\thelp
    
    \n\n\tUsage: perl LR_validate.pl [INPUTS] -5Pp <5' CAGE peak file (bed10 format from peak_caller_from_wigs.pl)> -mcde <Minimum CAGE peak depth> -mcdi <Maximum distance between CAGE peak and start of LR read> -3Pp <3' peak file> -3Pde <Minimum read depth of 3' peak clusters> -3Pdi <Maximum 3' distance between LR end and 3' peaks> -minSJ <Minimum number of validating short read splice junctions detected> -SJt <Splice junction tab file from STAR alignment of short reads> -ATG <known ATG start site file (bed6) (default = no known starts)> -f <genome fasta file> -LR <Input long read bed12 file>
    
    \n\n\tExample: perl LR_validate.pl -5Pp /PATH/5Ppeakfile.bed -mcde 10 -mcdi 2 -3Pp /PATH/3Ppeakfile.bed -3Pde 10 -3Pdi 10 -minSJ 2 -SJt /PATH/SJfile.tab -f /PATH/genome.fa -LR /PATH/LRfile.bed\n\n";
    exit;
}

sub options {
    if (scalar @ARGV == 0) {
        program_info;
        exit;
    }    
    for (my $i=0; $i < scalar @ARGV; $i++) {
        if ($ARGV[$i] eq "\-5Pp") {
            $CAGE_cons_file = $ARGV[$i+1];
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
        elsif ($ARGV[$i] eq "\-ATG") {
            $known_atg_file = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-f") {
            $genome_fasta = $ARGV[$i+1];
        }
        elsif ($ARGV[$i] eq "\-LR") {
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
        print "\nCAGE file not defined!\n";
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
        print "\nMaximum distance between 3' LR read and 3' cluster not defined!\n";
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
    elsif (not defined($genome_fasta)) {
        print "\nGenome fasta file not defined!\n";
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
    print OUT2 "CAGE consensus file: ", basename($CAGE_cons_file), "\nMin CAGE depth: ", $min_CAGE_depth, "\nMax distance from ONT start to CAGE peak: ", $max_CAGE_distance, "\nThree prime end peak file: ", basename($threeP_file), "\nMinimum 3Prime peak depth: ", $min_3P_depth, "\nMax distance from LR end to 3P peak: ", $max_3P_distance, "\nMinimum splice junction reads: ", $min_SJ_reads, "\nGenome fasta file: ", basename($genome_fasta), "\nShort read SJ tab file: ", basename($short_read_SJ_tab_file), "\nKnown start site file: ", basename($known_atg_file), "\nLR input file: ", basename($ONT_file), "\n";

    $negative_out_file = $out_directory."\/".$basename."_negative_splice_junction_reads.bed";
    open(OUT3, ">$negative_out_file") or die "couldn't open summary file";
}

sub fa_unwrapper {
    #my $datetime = DateTime->now;   
    print "\nProcessing genome fasta file...\n\n"; #(",$datetime, ")\n\n";
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
}

sub process_CAGE_cons_file {
    #my $datetime = DateTime->now; 
    print "Processing CAGE file...\n\n"; #(",$datetime, ")\n\n";
    my $CAGE_cons_file_pos = $CAGE_cons_file."_OUT_POS.bed";
    my $CAGE_cons_file_neg = $CAGE_cons_file.".OUT_NEG.bed";

    open(INF, "<$CAGE_cons_file") or die "Couldn't open input file '$CAGE_cons_file': $!";
    open(OUT_POS, ">$CAGE_cons_file_pos") or die "Couldn't open output file '$CAGE_cons_file_pos': $!";
    open(OUT_NEG, ">$CAGE_cons_file_neg") or die "Couldn't open output file '$CAGE_cons_file_neg': $!";

    while (my $line = <INF>) {
        chomp($line);
        $line =~ s/\"//g; 
        my @split_line = split("\t", $line);

        if ($split_line[5] eq "+") {
           print OUT_POS "$line\n";
        }
        elsif ($split_line[5] eq "-") {
            print OUT_NEG "$line\n";
        }
    }
    close(INF);
    close(OUT_POS);
    close(OUT_NEG);

    my $CAGE_cons_file_pos_sorted = $CAGE_cons_file_pos;
    $CAGE_cons_file_pos_sorted =~ s/\.bed$/_sortedforvalidation.bed/;
    `sort -V -k 1,1 -k 2,2n -k 3,3n  $CAGE_cons_file_pos > $CAGE_cons_file_pos_sorted`;
    `rm $CAGE_cons_file_pos`;
    $CAGE_cons_file_pos = $CAGE_cons_file_pos_sorted;

    my $CAGE_cons_file_neg_sorted = $CAGE_cons_file_neg;
    $CAGE_cons_file_neg_sorted =~ s/\.bed$/_sortedforvalidation.bed/;
    `sort -V -k 1,1 -k 2,2n -k 3,3n   $CAGE_cons_file_neg > $CAGE_cons_file_neg_sorted`;
    `rm $CAGE_cons_file_neg`;
    $CAGE_cons_file_neg = $CAGE_cons_file_neg_sorted;

    open(INF, "<$CAGE_cons_file_pos") or die "Couldn't open input file '$CAGE_cons_file_pos': $!";
    while (my $line = <INF>) {
        chomp($line);
        push(@CAGE_peaks_pos, $line)
    }
    close(INF);

    open(INF,"<$CAGE_cons_file_neg") or die "Couldn't open input file '$CAGE_cons_file_neg': $!";
    while (my $line = <INF>) {
        chomp($line);
        push(@CAGE_peaks_neg, $line);
    }
    close(INF);

    if (!@CAGE_peaks_pos) {
        warn "Warning: @CAGE_peaks_pos is empty or undefined.";
    }if (!@CAGE_peaks_neg) {
        warn "Warning: @CAGE_peaks_neg is empty or undefined.";
    }
    `rm $CAGE_cons_file_pos`;
    `rm $CAGE_cons_file_neg`;
}

sub process_3P_file {
    #my $datetime = DateTime->now; 
    print "Processing 3P file...\n\n"; #(",$datetime, ")\n\n";
    my @threeP_reads_pos;
    my @threeP_reads_neg;
    $threeP_file_pos = "$threeP_file.OUT_POS.bed";
    $threeP_file_neg = "$threeP_file.OUT_NEG.bed";
    
    open(INF, "<$threeP_file") or die "Couldn't open input file '$threeP_file': $!";
    open(OUT_POS, ">$threeP_file_pos") or die "Couldn't open output file '$threeP_file_pos': $!";
    open(OUT_NEG, ">$threeP_file_neg") or die "Couldn't open output file '$threeP_file_neg': $!";

    while (my $line = <INF>) {
        chomp($line);
        $line =~ s/\"//g; 
        my @split_line = split("\t", $line);
        if ($split_line[5] eq "+") {
           print OUT_POS "$line\n";
        }
        elsif ($split_line[5] eq "-") {
            print OUT_NEG "$line\n";
        }
    }
    close(INF);
    close(OUT_POS);
    close(OUT_NEG);

    $threeP_file_pos_sorted = $threeP_file_pos;
    $threeP_file_pos_sorted =~ s/\.bed$/_sortedforvalidation.bed/;
    `sort -V -k 1,1 -k 2,2n -k 3,3n  $threeP_file_pos > $threeP_file_pos_sorted`;
    `rm $threeP_file_pos`;
    $threeP_file_pos = $threeP_file_pos_sorted;

    $threeP_file_neg_sorted = $threeP_file_neg;
    $threeP_file_neg_sorted =~ s/\.bed$/_sortedforvalidation.bed/;
    `sort -V -k 1,1 -k 2,2n -k 3,3n $threeP_file_neg > $threeP_file_neg_sorted`;
    `rm $threeP_file_neg`;
    $threeP_file_neg = $threeP_file_neg_sorted;

    open(INF, "<$threeP_file_pos") or die "Couldn't open input file '$threeP_file_pos': $!";
    while (my $line = <INF>) {
        chomp($line);
        push(@threeP_peaks_pos, $line)
    }
    close(INF);

    open(INF,"<$threeP_file_neg") or die "Couldn't open input file '$threeP_file_neg': $!";
    while (my $line = <INF>) {
        chomp($line);
        push(@threeP_peaks_neg, $line);
    }
    close(INF);

    if (!@threeP_peaks_pos) {
        warn "Warning: @threeP_peaks_pos is empty or undefined.";
    }if (!@threeP_peaks_neg) {
        warn "Warning: @threeP_peaks_neg is empty or undefined.";
    }
    `rm $threeP_file_pos`;
    `rm $threeP_file_neg`;
}

sub process_SJTab_file {
    #my $datetime = DateTime->now; 
    print "Processing SJ file...\n\n"; #(",$datetime, ")\n\n";
    (my $short_read_SJ_tab_file_pos = $short_read_SJ_tab_file) =~ s/\.tab$//;
    $short_read_SJ_tab_file_pos .= "_pos.tab";
    (my $short_read_SJ_tab_file_neg = $short_read_SJ_tab_file) =~ s/\.tab$//;
    $short_read_SJ_tab_file_neg .= "_neg.tab";
    my @splice_junc;

    open(INF, "<$short_read_SJ_tab_file") or die "Couldn't open input file '$short_read_SJ_tab_file': $!";
    open(OUT_POS, ">$short_read_SJ_tab_file_pos") or die "Couldn't open output file '$short_read_SJ_tab_file_pos': $!";
    open(OUT_NEG, ">$short_read_SJ_tab_file_neg") or die "Couldn't open output file '$short_read_SJ_tab_file_neg': $!";

    while (my $line = <INF>) {
        chomp($line);
        $line =~ s/\"//g;
        my @split_line = split("\t", $line);
        if ($split_line[3] eq "1") {
           print OUT_POS "$line\n";
        }
        elsif ($split_line[3] eq "2") {
            print OUT_NEG "$line\n";
        }
    }
    close(INF);
    close(OUT_POS);
    close(OUT_NEG);

    my $short_read_SJ_tab_file_pos_sorted = $short_read_SJ_tab_file_pos;
    $short_read_SJ_tab_file_pos_sorted =~ s/(\.\w+)?$/_sortedforvalidation.tab/;
    `sort -V -k 1,1 -k 2,2n -k 3,3n $short_read_SJ_tab_file_pos > $short_read_SJ_tab_file_pos_sorted`;
    `rm $short_read_SJ_tab_file_pos`;
    $short_read_SJ_tab_file_pos = $short_read_SJ_tab_file_pos_sorted;
    my $short_read_SJ_tab_file_neg_sorted = $short_read_SJ_tab_file_neg;
    $short_read_SJ_tab_file_neg_sorted =~ s/(\.\w+)?$/_sortedforvalidation.tab/;
    `sort -V -k 1,1 -k 3,3n -k 2,2n $short_read_SJ_tab_file_neg > $short_read_SJ_tab_file_neg_sorted`;
    `rm $short_read_SJ_tab_file_neg`;
    $short_read_SJ_tab_file_neg = $short_read_SJ_tab_file_neg_sorted;

    open(INF, "<$short_read_SJ_tab_file_pos") or die "couldn't open genome fasta file";
    while(my $line = <INF>) {
        chomp($line);
        my @split_line = split("\t", $line);
        if($split_line[6] > $min_SJ_reads){
            my $bed_line = join("\t", @split_line[0..2])."\tjunc\t".$split_line[6]."\t\+\n";
            push(@splice_junc_pos, $bed_line);
        }
    }
    close(INF);

    open(INF, "<$short_read_SJ_tab_file_neg") or die "couldn't open genome fasta file";
    while(my $line = <INF>) {
        chomp($line);
        my @split_line = split("\t", $line);
        if($split_line[6] > $min_SJ_reads){
            my $bed_line = join("\t", @split_line[0..2])."\tjunc\t".$split_line[6]."\t\-\n";
            push(@splice_junc_neg, $bed_line);
        }
    }
    close(INF);
    `rm $short_read_SJ_tab_file_pos`;
    `rm $short_read_SJ_tab_file_neg`;
}

sub process_ONT_file {
    #my $datetime = DateTime->now; 
    print "Processing long read file...\n\n"; #(",$datetime, ")\n\n";
    my $ont_file_pos = $ONT_file."_LR_POS.bed";
    my $ont_file_neg = $ONT_file."_LR_NEG.bed";
    my $ontfile_pos_sorted;
    my $ontfile_neg_sorted;
   
    open(INF, "<$ONT_file") or die "Couldn't open input file '$ONT_file': $!";
    open(OUT_POS, ">$ont_file_pos") or die "Couldn't open output file '$ont_file_pos': $!";
    open(OUT_NEG, ">$ont_file_neg") or die "Couldn't open output file '$ont_file_neg': $!";

    while (my $line = <INF>) {
        chomp($line);
        $line =~ s/\"//g;
        my @split_line = split("\t", $line);
        if ($split_line[5] eq "+") {
           print OUT_POS "$line\n";
        }
        elsif ($split_line[5] eq "-") {
            print OUT_NEG "$line\n";
        }
    }
    close(INF);
    close(OUT_POS);
    close(OUT_NEG);

    $ontfile_pos_sorted = $ont_file_pos;
    $ontfile_pos_sorted =~ s/\.bed$/_sortedforvalidation.bed/;
    `sort -V -k 1,1 -k 2,2n -k 3,3n  $ont_file_pos > $ontfile_pos_sorted`;
    `rm $ont_file_pos`;
    $ont_file_pos = $ontfile_pos_sorted;

    $ontfile_neg_sorted = $ont_file_neg;
    $ontfile_neg_sorted =~ s/\.bed$/_sortedforvalidation.bed/;
    `sort -V -k 1,1 -k 3,3n -k 2,2n   $ont_file_neg > $ontfile_neg_sorted`;
    `rm $ont_file_neg`;
    $ont_file_neg = $ontfile_neg_sorted;

    open(INF, "<$ont_file_pos") or die "Couldn't open input file '$ont_file_pos': $!";
    while (my $line = <INF>) {
        chomp($line);
        push(@ONT_reads_pos, $line)
    }
    close(INF);

    open(INF,"<$ont_file_neg") or die "Couldn't open input file '$ont_file_neg': $!";
    while (my $line = <INF>) {
        chomp($line);
        push(@ONT_reads_neg, $line);
    }
    close(INF);

    if (!@ONT_reads_pos) {
        warn "Warning: @ONT_reads_pos is empty or undefined.";
    }
    if (!@ONT_reads_neg) {
        warn "Warning: @ONT_reads_neg is empty or undefined.";
    }
    `rm $ont_file_pos`;
    `rm $ont_file_neg`;
}

sub ONT_CAGE_overlap {
    my $count = 0;
    my $hit_count = 0;
    my $k = 0;
    my $j = 0;

    #my $datetime = DateTime->now; 
    print "Determining 5' end CAGE overlap with LRs (may take a few minutes)...\n\n"; #(",$datetime, ")\n\n";
    my $pos_output_file = "$ONT_file.CAGE_validated_pos.bed";
    open(OUT, ">$pos_output_file") or die or die "Couldn't open output file '$pos_output_file': $!";
    for (; $j < scalar(@ONT_reads_pos); $j++) {
        $count++;
        my @split_line = split("\t", $ONT_reads_pos[$j]);   
        for (my $i=$k; $i < scalar(@CAGE_peaks_pos); $i++) {
            my @split_CAGE_line = split("\t", $CAGE_peaks_pos[$i]);
            my @split_CAGE_info = split("\;", $split_CAGE_line[3]);
            if($split_CAGE_line[0] eq $split_line[0]) {
                if($split_line[1] >= $split_CAGE_line[1]-$max_CAGE_distance+1 and $split_line[1] <= $split_CAGE_line[2]+$max_CAGE_distance+1) {
                    my $line_plus_CAGE_value = join("\t", @split_line[0..3])."\t".$split_line[4]."_CAGE_val-".abs($split_CAGE_line[4])."-coord-".$split_CAGE_info[1]."\t".join("\t", @split_line[5..7])."\t".$split_CAGE_line[8]."\t".join("\t", @split_line[9..11]);
                    print OUT "$line_plus_CAGE_value\n";
                    $hit_count++;
                    if($i>2) {
                        $k = $i-2;
                    }
                    last;
                }
                elsif($split_CAGE_line[1] > $split_line[1]+100) {
                    last;
                }
            }
        }
    }
    close OUT;

    my $neg_output_file = "$ONT_file.CAGE_validated_neg.bed";
    open(OUT, ">$neg_output_file") or die or die "Couldn't open output file '$neg_output_file': $!";

    $j = 0;
    $k = 0;
    for (; $j < scalar(@ONT_reads_neg); $j++) {
        $count++; 
        my @split_line = split("\t", $ONT_reads_neg[$j]);
        for (my $i=$k; $i < scalar(@CAGE_peaks_neg); $i++) {
            my @split_CAGE_line = split("\t", $CAGE_peaks_neg[$i]);
            my @split_CAGE_info = split("\;", $split_CAGE_line[3]);
            if($split_CAGE_line[0] eq $split_line[0]) {
                if($split_line[2] >= $split_CAGE_line[1]-$max_CAGE_distance-1 and $split_line[2] <= $split_CAGE_line[2]+$max_CAGE_distance-1) {
                    my $line_plus_CAGE_value = join("\t", @split_line[0..3])."\t".$split_line[4]."_CAGE_val-".abs($split_CAGE_line[4])."-coord-".$split_CAGE_info[2]."\t".join("\t", @split_line[5..7])."\t".$split_CAGE_line[8]."\t".join("\t", @split_line[9..11]);
                    print OUT "$line_plus_CAGE_value\n";
                    $hit_count++;
                    if($i>2) {
                        $k = $i-2;
                    }
                    last;
                }
                elsif($split_CAGE_line[1] > $split_line[2]+100) {
                    last;
                }
            }
        }
    }
    close OUT;
    print OUT2 "\nFraction of LRs with CAGE peak validation = ", $hit_count/$count, "\n";

    my $pos_output_file_sorted = $pos_output_file;
    $pos_output_file_sorted =~ s/\.bed$/_sortedforvalidation.bed/;
    `sort -V -k 1,1 -k 3,3n -k 2,2n $pos_output_file > $pos_output_file_sorted`;
    `rm $pos_output_file`;
    $pos_output_file = $pos_output_file_sorted;

    my $neg_output_file_sorted = $neg_output_file;
    $neg_output_file_sorted =~ s/\.bed$/_sortedforvalidation.bed/;
    `sort -V -k 1,1 -k 2,2n -k 3,3n $neg_output_file > $neg_output_file_sorted`;
    `rm $neg_output_file`;
    $neg_output_file = $neg_output_file_sorted;

    open(INF, "<$pos_output_file") or die "Couldn't open input file '$pos_output_file': $!";
    while (my $line = <INF>) {
        chomp($line);
        push(@CAGE_validated_pos, $line)
    }
    close(INF);

    open(INF,"<$neg_output_file") or die "Couldn't open input file '$neg_output_file': $!";
    while (my $line = <INF>) {
        chomp($line);
        push(@CAGE_validated_neg, $line);
    }
    close(INF);

    if (!@CAGE_validated_pos) {
        warn "Warning: @CAGE_validated_pos is empty or undefined.";
    }if (!@CAGE_validated_neg) {
        warn "Warning: @CAGE_validated_neg is empty or undefined.";
    }
    `rm $pos_output_file`;
    `rm $neg_output_file`;
}

sub ONT_3_prime_overlap {
    #my $datetime = DateTime->now; 
    print "Determining 3' end overlap with CAGE validated LRs...\n\n"; #(",$datetime, ")\n\n";
    my $count = 0;
    my $hit_count = 0;
    my $k = 0;
    my $j = 0;

    for (; $j < scalar(@CAGE_validated_pos); $j++) {
        $count++;
        my @split_line = split("\t", $CAGE_validated_pos[$j]);
        for (my $i=$k; $i < scalar(@threeP_peaks_pos); $i++) {
            my @split_threeP_peaks_line = split("\t", $threeP_peaks_pos[$i]);
            my @split_threeP_peaks_info = split("\;", $split_threeP_peaks_line[3]); 
            if ($split_threeP_peaks_line[0] eq $split_line[0]) { 
                if ($split_line[2] >= $split_threeP_peaks_line[1] - $max_3P_distance and $split_line[2] <= $split_threeP_peaks_line[2] + $max_3P_distance) {
                    my $push_CAGE_validated = join("\t", @split_line[0..4]) . "-3P_val-" . abs($split_threeP_peaks_line[4]) . "-coord-" . $split_threeP_peaks_info[2] . "\t" . join("\t", @split_line[5..11]);
                    push(@CAGE_and_3P_validated, $push_CAGE_validated);
                    $hit_count++;
                    if ($i > 2) {
                        $k = $i - 2;
                    }
                    last;
                }
                elsif ($split_threeP_peaks_line[1] > $split_line[2] + 100) {
                    last;
                }
            }
        }
    }
    $k = 0;
    $j = 0;

    for (; $j < scalar(@CAGE_validated_neg); $j++) {
        $count++;
        my @split_line = split("\t", $CAGE_validated_neg[$j]);
        for (my $i = $k; $i < scalar(@threeP_peaks_neg); $i++) {
            my @split_threeP_peaks_line = split("\t", $threeP_peaks_neg[$i]);
            my @split_threeP_peaks_info = split(";", $split_threeP_peaks_line[3]);
            if ($split_threeP_peaks_line[0] eq $split_line[0]) {
                if ($split_line[1] >= $split_threeP_peaks_line[1] - $max_3P_distance and $split_line[1] <= $split_threeP_peaks_line[2] + $max_3P_distance) {
                    my $push_CAGE_validated = join("\t", @split_line[0..4]) . "-3P_val-" . abs($split_threeP_peaks_line[4]) . "-coord-" . $split_threeP_peaks_info[1] . "\t" . join("\t", @split_line[5..11]);
                    push(@CAGE_and_3P_validated, $push_CAGE_validated);
                    $hit_count++;
                    if ($i > 2) {
                        $k = $i - 2;
                    }
                    last;
                }
                elsif ($split_threeP_peaks_line[1] > $split_line[1] + 100) {
                    last;
                }
            }
        }
    }
    print OUT2 "\nFraction of CAGE validated LRs with 3P validation = ", $hit_count / $count, "\n";
}

sub ONT_AAUAAA_motif {
    #my $datetime = DateTime->now; 
    print "Assessing presence of AAUAAA motifs in 5' and 3' validated LRs...\n\n"; #(",$datetime, ")\n\n";
    my $count = 0;
    my $hit_count = 0;
    my $hit_AUUAAA_count = 0;
    my $hit_AAUACA_count = 0;
    my $hit_GAUAAA_count = 0;

    my $output_file_pos = "$ONT_file.CAGE_and_3P_validated_AAUAAA_pos.bed";
    my $output_file_neg = "$ONT_file.CAGE_and_3P_validated_AAUAAA_neg.bed";

    open(my $OUT_POS, '>', $output_file_pos) or die "Couldn't open output file '$output_file_pos': $!";
    open(my $OUT_NEG, '>', $output_file_neg) or die "Couldn't open output file '$output_file_neg': $!";

    for (my $i = 0; $i < scalar(@CAGE_and_3P_validated); $i++) {
        $count++;
        my @split_transcript = split("\t", $CAGE_and_3P_validated[$i]);
        my $threeP_seq;
        for (my $j = 0; $j < scalar(@genome_fa); $j++) {
            if ($split_transcript[0] eq $genome_fa[$j]) {
                my $chr_seq = $genome_fa[$j+1];
                if ($split_transcript[5] eq "+") {
                    $threeP_seq = substr($chr_seq, $split_transcript[2] - 35, 30);
                    my $line;
                    if ($threeP_seq =~ m/AATAAA/) {
                        $line = join("\t", @split_transcript[0..11]) . "\tAAUAAA";
                        print $OUT_POS "$line\n";
                        $hit_count++;
                    }
                    elsif ($threeP_seq =~ m/ATTAAA/) {
                        $line = join("\t", @split_transcript[0..11]) . "\tAUUAAA";
                        print $OUT_POS "$line\n";
                        $hit_AUUAAA_count++;
                    }
                    elsif ($threeP_seq =~ m/AATACA/) {
                        $line = join("\t", @split_transcript[0..11]) . "\tAAUACA";
                        print $OUT_POS "$line\n";
                        $hit_AAUACA_count++;
                    }
                    elsif ($threeP_seq =~ m/GATAAA/) {
                        $line = join("\t", @split_transcript[0..11]) . "\tGAUAAA";
                        print $OUT_POS "$line\n";
                        $hit_GAUAAA_count++;
                    }
                    else {
                        $line = join("\t", @split_transcript[0..11]) . "\tNone";
                    }
                }
                elsif ($split_transcript[5] eq "-") {
                    $threeP_seq = substr($chr_seq, $split_transcript[1] + 5, 30);
                    my $line;
                    if ($threeP_seq =~ m/TTTATT/) {
                        $line = join("\t", @split_transcript[0..11]) . "\tAAUAAA";
                        print $OUT_NEG "$line\n";
                        $hit_count++;
                    }
                    elsif ($threeP_seq =~ m/TTTAAT/) {
                        $line = join("\t", @split_transcript[0..11]) . "\tAUUAAA";
                        print $OUT_NEG "$line\n";
                        $hit_AUUAAA_count++;
                    }
                    elsif ($threeP_seq =~ m/TGTATT/) {
                        $line = join("\t", @split_transcript[0..11]) . "\tAAUACA";
                        print $OUT_NEG "$line\n";
                        $hit_AAUACA_count++;
                    }
                    elsif ($threeP_seq =~ m/TTTATC/) {
                        $line = join("\t", @split_transcript[0..11]) . "\tGAUAAA";
                        print $OUT_NEG "$line\n";
                        $hit_GAUAAA_count++;
                    }
                    else {
                        $line = join("\t", @split_transcript[0..11]) . "\tNone";
                    }
                }
            }
        }
    }
    close($OUT_POS);
    close($OUT_NEG);

    print OUT2 "\nFraction of 3P and CAGE validated LRs with upstream AAUAAA = ", $hit_count / $count, "\n";
    print OUT2 "Fraction of 3P and CAGE validated LRs with upstream AUUAAA = ", $hit_AUUAAA_count / $count, "\n";
    print OUT2 "Fraction of 3P and CAGE validated LRs with upstream AAUACA = ", $hit_AAUACA_count / $count, "\n";
    print OUT2 "Fraction of 3P and CAGE validated LRs with upstream GAUAAA = ", $hit_GAUAAA_count / $count, "\n";

    my $output_file_pos_sorted = $output_file_pos;
    $output_file_pos_sorted =~ s/\.bed$/_sortedforvalidation.bed/;
    `sort -V -k 1,1 -k 2,2n -k 3,3n $output_file_pos > $output_file_pos_sorted`;
    `rm $output_file_pos`;
    $output_file_pos = $output_file_pos_sorted;

    my $output_file_neg_sorted = $output_file_neg;
    $output_file_neg_sorted =~ s/\.bed$/_sortedforvalidation.bed/;
    `sort -V -k 1,1 -k 3,3n -k 2,2n $output_file_neg > $output_file_neg_sorted`;
    `rm $output_file_neg`;
    $output_file_neg = $output_file_neg_sorted; 
    
    open(INF, "<$output_file_pos") or die "Couldn't open input file '$output_file_pos': $!";
    while (my $line = <INF>) {
        chomp($line);
        push(@CAGE_and_3P_validated_AAUAAA_pos, $line)
    }
    close(INF);

    open(INF,"<$output_file_neg") or die "Couldn't open input file '$output_file_neg': $!";
    while (my $line = <INF>) {
        chomp($line);
        push(@CAGE_and_3P_validated_AAUAAA_neg, $line);
    }
    close(INF);

    if (!@CAGE_and_3P_validated_AAUAAA_pos) {
        warn "Warning: @CAGE_and_3P_validated_AAUAAA_pos is empty or undefined.";
    }
    if (!@CAGE_and_3P_validated_AAUAAA_neg) {
        warn "Warning: @CAGE_and_3P_validated_AAUAAA_neg is empty or undefined.";
    }
    `rm $output_file_pos`;
    `rm $output_file_neg`;
}

sub validate_SJs {
    #my $datetime = DateTime->now; 
    print "Validating SJs in spliced LRs using short read data (STAR .tab file)...\n\n"; #(",$datetime, ")\n\n";
    my $i = 0;
    my $k = 0;
    for (my $j = 0; $j < scalar(@CAGE_and_3P_validated_AAUAAA_pos); $j++) {
        my @split_CAGE_and_3P_validated_AAUAAA_line_pos = split("\t", $CAGE_and_3P_validated_AAUAAA_pos[$j]);
        if ($split_CAGE_and_3P_validated_AAUAAA_line_pos[9] == 1) {
            my $line = $CAGE_and_3P_validated_AAUAAA_pos[$j] . "\tno splice";
            push(@CAGE_and_3P_validated_AAUAAA_SJ_val_pos, $line);
        }
        else {
            my @hit_array_pos = ();
            my @size_array_pos = split(",", $split_CAGE_and_3P_validated_AAUAAA_line_pos[10]);
            my @start_array_pos = split(",", $split_CAGE_and_3P_validated_AAUAAA_line_pos[11]);
            for (my $e = 0; $e < $split_CAGE_and_3P_validated_AAUAAA_line_pos[9] - 1; $e++) {
                my $intron_start_pos = $split_CAGE_and_3P_validated_AAUAAA_line_pos[1] + $start_array_pos[$e] + $size_array_pos[$e] + 1;
                my $intron_end_pos = $split_CAGE_and_3P_validated_AAUAAA_line_pos[1] + $start_array_pos[$e + 1];
                my $intron_hits_pos = 0;
                for ($i = $k; $i < scalar(@splice_junc_pos); $i++) {
                    my @splice_junc_line_pos = split("\t", $splice_junc_pos[$i]);
                    if ($split_CAGE_and_3P_validated_AAUAAA_line_pos[0] eq $splice_junc_line_pos[0]) {
                        if ($intron_start_pos == $splice_junc_line_pos[1] and $intron_end_pos == $splice_junc_line_pos[2]) {
                            $intron_hits_pos++;
                            if ($i > 2) {
                                $k = $i - 2;
                            }
                            last;
                        }
                        elsif ($splice_junc_line_pos[1] > $intron_end_pos + 100) {
                            last;
                        }
                    }
                }    
                if ($intron_hits_pos == 0) {
                    push(@hit_array_pos, "false");
                }
                else {
                    push(@hit_array_pos, "true");
                }
            }
            my $false = 0;
            for (my $e = 0; $e < scalar(@hit_array_pos); $e++) {
                if ($hit_array_pos[$e] eq "false") {
                    $false++;
                }
            }
            if ($false == 0) {
                my $line = $CAGE_and_3P_validated_AAUAAA_pos[$j] . "\t" . join(",", @hit_array_pos);
                push(@CAGE_and_3P_validated_AAUAAA_SJ_val_pos, $line);
            }
            else {
                my $line = $CAGE_and_3P_validated_AAUAAA_pos[$j] . "\t" . join(",", @hit_array_pos);
                print OUT3 $line, "\n";
            }
        } 
        for(my $z = $k; $k>=0; $z--) {
            my @backup_splice_junc_line_pos = split("\t", $splice_junc_pos[$z]);
            if($split_CAGE_and_3P_validated_AAUAAA_line_pos[0] ne $backup_splice_junc_line_pos[0] || ($split_CAGE_and_3P_validated_AAUAAA_line_pos[0] eq $backup_splice_junc_line_pos[0] and $split_CAGE_and_3P_validated_AAUAAA_line_pos[1] > $backup_splice_junc_line_pos[1])) {
                $k = $z;
                last;
            }
        }  
    }    
    $i = 0;
    $k = 0;
    for (my $j = 0; $j < scalar(@CAGE_and_3P_validated_AAUAAA_neg); $j++) {
        my @split_CAGE_and_3P_validated_AAUAAA_line_neg = split("\t", $CAGE_and_3P_validated_AAUAAA_neg[$j]);
        if ($split_CAGE_and_3P_validated_AAUAAA_line_neg[9] == 1) {
            my $line = $CAGE_and_3P_validated_AAUAAA_neg[$j] . "\tno splice";
            push(@CAGE_and_3P_validated_AAUAAA_SJ_val_neg, $line);
        }
        else {
            my @hit_array_neg = ();
            my @size_array_neg = split(",", $split_CAGE_and_3P_validated_AAUAAA_line_neg[10]);
            my @start_array_neg = split(",", $split_CAGE_and_3P_validated_AAUAAA_line_neg[11]);
            for (my $e = 0; $e < $split_CAGE_and_3P_validated_AAUAAA_line_neg[9] - 1; $e++) {
                my $intron_start_neg = $split_CAGE_and_3P_validated_AAUAAA_line_neg[1] + $start_array_neg[$e] + $size_array_neg[$e] + 1;
                my $intron_end_neg = $split_CAGE_and_3P_validated_AAUAAA_line_neg[1] + $start_array_neg[$e + 1];
                my $intron_hits_neg = 0;
                for (my $i = $k; $i < scalar(@splice_junc_neg); $i++) {
                    my @splice_junc_line_neg = split("\t", $splice_junc_neg[$i]);
                    if ($split_CAGE_and_3P_validated_AAUAAA_line_neg[0] eq $splice_junc_line_neg[0]) {
                        if ($intron_start_neg == $splice_junc_line_neg[1] and $intron_end_neg == $splice_junc_line_neg[2]) {
                            $intron_hits_neg++;
                            if ($i > 2) {
                                $k = $i - 2;
                            }
                            last;
                        }
                        elsif ($splice_junc_line_neg[1] > $intron_end_neg + 100) {
                            last;
                        }
                    }
                }    
                if ($intron_hits_neg == 0) {
                    push(@hit_array_neg, "false");
                }
                else {
                    push(@hit_array_neg, "true");
                }
            }
            my $false = 0;
            for (my $e = 0; $e < scalar(@hit_array_neg); $e++) {
                if ($hit_array_neg[$e] eq "false") {
                    $false++;
                }
            }
            if ($false == 0) {
                my $line = $CAGE_and_3P_validated_AAUAAA_neg[$j] . "\t" . join(",", @hit_array_neg);
                push(@CAGE_and_3P_validated_AAUAAA_SJ_val_neg, $line);
            }
            else {
                my $line = $CAGE_and_3P_validated_AAUAAA_neg[$j] . "\t" . join(",", @hit_array_neg);
                print OUT3 $line, "\n";
            }
        }  
        for(my $z = $k; $k>=0; $z--) {
            my @backup_splice_junc_line_neg = split("\t", $splice_junc_neg[$z]);
            if($split_CAGE_and_3P_validated_AAUAAA_line_neg[0] ne $backup_splice_junc_line_neg[0] || ($split_CAGE_and_3P_validated_AAUAAA_line_neg[0] eq $backup_splice_junc_line_neg[0] and $split_CAGE_and_3P_validated_AAUAAA_line_neg[1] > $backup_splice_junc_line_neg[1])) {
                $k = $z;
                last;
            }
        }  
    }    
}

sub print_CAGE_3P_validated_ONT_reads {
    #my $datetime = DateTime->now; 
    print "Printing validated LRs...\n\n"; #(",$datetime, ")\n\n";
    open (OUT, ">$out_file") or die "couldn't open output file";
    print OUT join("\n", @CAGE_and_3P_validated_AAUAAA_SJ_val_pos);
    print OUT join("\n", @CAGE_and_3P_validated_AAUAAA_SJ_val_neg);
    close(OUT);
}

sub consolidate_reads {
    #my $datetime = DateTime->now; 
    print "Consolidating reads and generating isoform annotation (may take a few minutes)...\n\n"; #(",$datetime, ")\n\n";
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
            if(not defined $thick_end) {
                $thick_end = $split_prev_line[2];
            }
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
    #my $datetime = DateTime->now; 
    print "Processing known ATG file...\n\n"; #(",$datetime, ")\n\n";
    my $pos_count = -1;
    my $neg_count = -1;
    if($known_atg_file eq "none") {
        my $pos = "chrnull\t0\t0\t.\t.\t+";
        push(@pos_atg_array, $pos);
        my $neg = "chrnull\t0\t0\t.\t.\t-";
        push(@neg_atg_array, $neg);
    }
    else {
        my $sorted_known_atg_file = $known_atg_file.".sorted\.bed";
        `sort -V -k 1,1 -k 2,2n -k 3,3n  $known_atg_file > $sorted_known_atg_file`;
        open (INF, "<$sorted_known_atg_file") or die "couldn't open input file";

        while (my $line = <INF>) {
	        chomp($line);
            my @split_line = split("\t", $line);
            if($split_line[5] eq "+") {
                push(@pos_atg_array, $line);
            }
            elsif($split_line[5] eq "-") {
                push(@neg_atg_array, $line);
            }
        }
        my $last_pos = "chrnull\t0\t0\t.\t.\t+";
        push(@pos_atg_array, $last_pos);
        my $last_neg = "chrnull\t0\t0\t.\t.\t-";
        push(@neg_atg_array, $last_neg);
        close(INF);
        `rm $sorted_known_atg_file`;
    }
}

sub find_reading_frames {
    my $line = $_[0];
    my @split_line = split("\t", $line);
    my @split_sizes = split("\,", $split_line[10]);
    my @split_starts = split("\,", $split_line[11]);
    my $block_number = $split_line[9];
    my $chr = $split_line[0];
    if($prev_reading_frame_chr ne $chr) {
        for(my $p=0; $p<scalar(@pos_atg_array); $p++) {
            my @split_pos = split("\t", $pos_atg_array[$p]);
            if($split_pos[0] eq $chr) {
                $pos_atg_array_chr_start = $p;
                last;
            }
        }
        for(my $n=0; $n<scalar(@neg_atg_array); $n++) {
            my @split_neg = split("\t", $neg_atg_array[$n]);
            if($split_neg[0] eq $chr) {
                $neg_atg_array_chr_start = $n;
                last;
            }
        }
    }
    $prev_reading_frame_chr = $chr;
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
                    for(my $atg=$pos_atg_array_chr_start; $atg < scalar(@pos_atg_array); $atg++) {
                        my @split_pos_atg_array = split("\t", $pos_atg_array[$atg]);
                        if($split_pos_atg_array[2] > $split_line[2]) {
                            last;
                        }
                        if($for_loop_count == 0 and $chr eq $split_pos_atg_array[0] and $split_pos_atg_array[1] > $block_start and $split_pos_atg_array[1] <= $block_end) {
                            $known_atg_coord = $split_pos_atg_array[1];
                            $for_loop_count++;
                        }
                    }
                }
                elsif($split_line[5] eq "-") {
                    for(my $atg=$neg_atg_array_chr_start; $atg < scalar(@neg_atg_array); $atg++) {
                        my @split_neg_atg_array = split("\t", $neg_atg_array[$atg]);
                        if($split_neg_atg_array[2] > $split_line[2]) {
                            last;
                        }
                        if($chr eq $split_neg_atg_array[0] and $split_neg_atg_array[1] >= $block_start and $split_neg_atg_array[1] < $block_end) {
                            $known_atg_coord = $split_neg_atg_array[1];
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
                        my $length;
                        if($ATG_position eq "null") {
                            $length = 0;
                        }
                        else {
                            $length = $end_position - $ATG_position;
                        }
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
                    if(not defined $split_start[1]) {
                        print $ATG_positions[$i], "\n";
                        print $known_atg_coord, "\n";
                        print $known_atg_transcript_coord, "\n";
                        print $line, "\n";
                        print join("\t", @ATG_positions), "\n";
                        print $transcript_seq, "\n";
                    }

                    my $ATG_position = $split_start[1]+1;
                    
                    for (my $j=scalar(@TERM_positions)-1; $j >= 0; $j--) {
                        my @split_TERM = split("\;", $TERM_positions[$j]);
                        my $end_position = $split_TERM[1]+1;
                        my $length;
                        if($ATG_position eq "null") {
                            $length = 0;
                        }
                        else {
                            $length = $ATG_position - $end_position;
                        }
                        if ($length % 3 == 0) {
                            if($length >= 300 || ($known_atg_coord ne "null" and $length > 0)) {
                                
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
            last;
        }
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
            if(substr($sequence, $-[0], ($+[0]) - ($-[0])) eq "CATG") {
                push (@output, join("\;", $-[0], $+[0]-2));
            }
            elsif(substr($sequence, $-[0], ($+[0]) - ($-[0])) eq "CCAT") {
                push (@output, join("\;", $-[0]+1, $+[0]-1));

            }
            
        }
        else {
            push (@output, join("\;", $-[0], $+[0]-1));
        }
    }
    return @output;
}

sub deduplicate_from_validated_single_entry_file {
    #my $datetime = DateTime->now; 
    print "Deduplicating validated single entry file...\n\n"; #(",$datetime, ")\n\n";
    my $file = $out_directory."\/VALIDATED_transcripts_".$basename.".single_entry\.bed";
    my $sorted_file = $out_directory."\/VALIDATED_transcripts_".$basename.".single_entry_sorted";
    my $output_file_pos = $out_directory."\/VALIDATED_transcripts_".$basename.".single_entry_deduplicated_pos.bed";
    my $output_file_neg = $out_directory."\/VALIDATED_transcripts_".$basename.".single_entry_deduplicated_neg.bed";

    `sort -k 4,4 $file > $sorted_file`;

    open(INF, "<$sorted_file") or die "couldn't open genome fasta file";
    open(OUTp, ">$output_file_pos") or die "couldn't open summary file";
    open(OUTn, ">$output_file_neg") or die "couldn't open summary file";
    my $prev_coord1 = 0;
    my $prev_coord2 = 0;
    my $prev_transcript_ID = "null";
    my $prev_CAGE_depth = 0;
    my $prev_strand = "null";
    my $prev_line = "null";
    while(my $line = <INF>) {
        chomp($line);
        my @split_line = split("\t", $line);
        if($. == 1) {
            $prev_coord1 = $split_line[1];
            $prev_coord2 = $split_line[2];
            $prev_transcript_ID = $split_line[3];
            $prev_CAGE_depth = $split_line[14];
            $prev_strand = $split_line[5];
            $prev_line = $line;
        }
        elsif(eof) {
            if($prev_transcript_ID ne $split_line[3]) {
                if($prev_strand eq "+") {
                    print OUTp $prev_line, "\n";
                }
                elsif($prev_strand eq "-") {
                    print OUTn $prev_line, "\n";
                }
                if($split_line[5] eq "+") {
                    print OUTp $line, "\n";
                }
                elsif($split_line[5] eq "-") {
                    print OUTn $line, "\n";
                }
            }
            else {
                if ($split_line[5] eq "+") {
                    if($split_line[2] >= $prev_coord2 and $split_line[14] >= $prev_CAGE_depth) {
                        print OUTp $line, "\n";
                    }
                    else {
                        print OUTp $prev_line, "\n";
                    }
                }
                elsif ($split_line[5] eq "-") {
                    if($split_line[1] <= $prev_coord2 and $split_line[14] >= $prev_CAGE_depth) {
                        print OUTn $line, "\n";
                    }
                    else {
                        print OUTn $prev_line, "\n";
                    }
                }
            }
        }
        elsif($prev_transcript_ID ne $split_line[3]) {
            if($prev_strand eq "+") {
                print OUTp $prev_line, "\n";
            }
            elsif($prev_strand eq "-") {
                print OUTn $prev_line, "\n";
            }
            $prev_coord1 = $split_line[1];
            $prev_coord2 = $split_line[2];
            $prev_transcript_ID = $split_line[3];
            $prev_strand = $split_line[5];
            $prev_CAGE_depth = $split_line[14];
            $prev_line = $line;
        }
        elsif ($split_line[5] eq "+") {
            if($split_line[2] >= $prev_coord2 and $split_line[14] >= $prev_CAGE_depth) {
                $prev_coord1 = $split_line[1];
                $prev_coord2 = $split_line[2];
                $prev_transcript_ID = $split_line[3];
                $prev_strand = $split_line[5];
                $prev_CAGE_depth = $split_line[14];
                $prev_line = $line;
            }
        }
        elsif ($split_line[5] eq "-") {
            if($split_line[1] <= $prev_coord1 and $split_line[14] >= $prev_CAGE_depth) {
                $prev_coord1 = $split_line[1];
                $prev_coord2 = $split_line[2];
                $prev_transcript_ID = $split_line[3];
                $prev_strand = $split_line[5];
                $prev_CAGE_depth = $split_line[14];
                $prev_line = $line;
            }
        }
    }
    close(INF);
    close(OUTp);
    close(OUTn);
    my $output_file_pos_sorted = $out_directory."\/VALIDATED_transcripts_".$basename.".single_entry_deduplicated_pos_sorted.bed";
    my $output_file_neg_sorted = $out_directory."\/VALIDATED_transcripts_".$basename.".single_entry_deduplicated_neg_sorted.bed";

    `sort -V -k 1,1 -k 2,2n -k 3,3n $output_file_pos > $output_file_pos_sorted`;
    `sort -V -k 1,1 -k 3,3n -k 2,2n $output_file_neg > $output_file_neg_sorted`;
    `rm $sorted_file`;
    `rm $output_file_pos`;
    `rm $output_file_neg`;
}

sub final_clean_merge_competing_dual_3p_ends {
    #my $datetime = DateTime->now; 
    print "Merging competing dual 3p ends...\n\n"; #(",$datetime, ")\n\n";
    my $file_pos = $out_directory."\/VALIDATED_transcripts_".$basename.".single_entry_deduplicated_pos_sorted.bed";
    my $output_file = $file_pos;
    $output_file =~ s/\.single_entry_deduplicated_pos_sorted\.bed//g;
    $output_file =~ s/VALIDATED_//g;
    $output_file = $output_file."_confirmed.bed";

    my $file_neg = $out_directory."\/VALIDATED_transcripts_".$basename.".single_entry_deduplicated_neg_sorted.bed";
    
    my @array_pos = ();
    open (INFp, "<$file_pos") or die "couldn't open input file";
    while (my $line = <INFp>) {
	    chomp($line);
        push(@array_pos, $line);
    }
    close(INFp);

    open (OUT, ">$output_file") or die "couldn't open output file";
    my @new_array_pos = ();
    my @skip_elements_pos = ();
    my $z = 0;
    for (my $i = 0; $i<scalar(@array_pos); $i++) {
        my @split_i = split("\t", $array_pos[$i]);
        my $ONT_count = $split_i[16];
        my $five_P_hit = 0;
        for (my $j = $z; $j<scalar(@array_pos); $j++) {
            my @split_j = split("\t", $array_pos[$j]);
            my $hit = 0;
            if ($i==$j) {
                $five_P_hit++;
                if($five_P_hit == 1) {
                    if($j > 2) {
                        $z = $j-2;
                    }
                }
                next;
            }
            elsif($split_i[1] == $split_j[1]) {
                $five_P_hit++;
                if($split_i[2] > $split_j[2] and $split_i[2]-37 < $split_j[2] and $split_i[9] == $split_j[9]) {
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
                }
                if($hit == $split_i[9]) {
                    $ONT_count = $ONT_count + $split_j[16];
                    push(@skip_elements_pos, $j);
                }
            }
            if($five_P_hit == 1) {
                if($j > 2) {
                    $z = $j-2;
                }
            }
            if($split_i[1] < $split_j[1]) {
                last;
            }
        }
        my $element = join("\t", @split_i[0..15])."\t".$ONT_count;
        push(@new_array_pos, $element);
    }

    for(my $a=0; $a<scalar(@new_array_pos); $a++) {
        my $hit = 0;
        for(my $b=0; $b<scalar(@skip_elements_pos); $b++) {
            if($a == $skip_elements_pos[$b]) {
                $hit++;
            }
        }
        if($hit == 0) {
            print OUT $new_array_pos[$a], "\n";
        }
    }

    my @array_neg = ();
    open (INFn, "<$file_neg") or die "couldn't open input file";
    while (my $line = <INFn>) {
	    chomp($line);
        push(@array_neg, $line);
    }
    close(INFn);

    my @new_array_neg = ();
    my @skip_elements_neg = ();
    $z = 0;
    for (my $i = 0; $i<scalar(@array_neg); $i++) {
        my @split_i = split("\t", $array_neg[$i]);
        my $ONT_count = $split_i[16];
        my $five_P_hit = 0;
        for (my $j = $z; $j<scalar(@array_neg); $j++) {
            my @split_j = split("\t", $array_neg[$j]);
            my $hit = 0;
            if ($i==$j) {
                $five_P_hit++;
                if($five_P_hit == 1) {
                    if($j > 2) {
                        $z = $j-2;
                    }
                }
                next;
            }
            elsif($split_i[2] == $split_j[2]) {
                $five_P_hit++;
                if($split_i[1] < $split_j[1] and $split_i[1]+37 > $split_j[1] and $split_i[9] == $split_j[9]) {
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
                        push(@skip_elements_neg, $j);
                    }
                }
            }
            if($five_P_hit == 1) {
                if ($j > 2) {
                    $z = $j-2;
                }
            }
            if($split_i[2] < $split_j[2]) {
                last;
            }
        }
        my $element = join("\t", @split_i[0..15])."\t".$ONT_count;
        push(@new_array_neg, $element);
    }

    for(my $a=0; $a<scalar(@new_array_neg); $a++) {
        my $hit = 0;
        for(my $b=0; $b<scalar(@skip_elements_neg); $b++) {
            if($a == $skip_elements_neg[$b]) {
                $hit++;
            }
        }
        if($hit == 0) {
            print OUT $new_array_neg[$a], "\n";
        }
    }
    close(OUT);

    my $output_file_sorted = $output_file."_confirmed_sorted.bed";
    `sort -V -k 1,1 -k 2,2n -k3,3n $output_file > $output_file_sorted`;
    `mv $output_file_sorted $output_file`;

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
process_SJTab_file;
process_ONT_file;
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
#my $datetime = DateTime->now; 
print "DONE!!!\n\n"; #(",$datetime, ")\n\n";
