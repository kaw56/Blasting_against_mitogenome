#!usr/bin/perl
# blast_annotated_genome.pl
# version 2 for blasting against mitochondrial genome
# this time using a database of mitochondrial annotations
use warnings;
use strict;

use Bio::SearchIO;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlastPlus;

# code reused from version 1.5 (annotate hits)

##########
# set up #
##########

die 'usage: annotate_hits.pl <query> <query format> <database for comparison>' 
    unless (@ARGV == 3);

# filenames    
my $sequence_filename   = shift;
my $format              = shift;
my $database_filename   = shift;
my $report_filename     = "$sequence_filename-RESULTS.txt";

# check whether the file to search against exists in this directory...
die "Check that the file to be searched against is in this directory\n" 
    unless (-e $database_filename);

# open output for printing results...
open (my $report, ">", $report_filename) 
        or die "Can't open $report_filename for writing, $!";
        
# create seqIO object...
my $seq_in  = Bio::SeqIO->new(
                            -format => $format,
                            -file   => $sequence_filename,
    );

# create a blast factory, set up database as subject...
my $factory = Bio::Tools::Run::StandAloneBlastPlus->new(
                            -db_data => $database_filename,
                            -create  => 1,
    ); 

#########
# blast # 
#########


# run a blast
# if a hit save the report
# also keep a list of the names of hits as a summary results output
