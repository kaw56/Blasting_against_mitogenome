#!usr/bin/perl
# ep_mito_blast.pl
#version 1 for blasting against mitochondrial genome
use warnings;
use strict;

use Bio::SearchIO;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlastPlus;

die "usage: ep_mito_blast <fasta database> <output>\n" unless (@ARGV == 2);

my $input_file  = shift;
my $output      = shift;
my $database    = 'partial_Ep_mitogenome';

# create a fasta seq object...
my $seq_in  = Bio::SeqIO->new(
                            -format => 'fasta',
                            -file   => $input_file,
    );

# check whether the file to search against exists...
die "Check that the database/file to be searched against exists!\n" 
                                                    unless (-e $database);

# create a blast factory, set up ep mito as subject...
my $factory = Bio::Tools::Run::StandAloneBlastPlus->new(
                            -db_data => $database,
                            -create  => 1,
    ); 

# open output for printing results...
open (my $outfile, ">", $output) or die "Can't open $output for writing, $!";

# blast each sequence in the input against ep mito...
my $seq;
while ($seq = $seq_in->next_seq() ) {
    
    # blast against ep mito and print result...
    my $blast_report = $factory->blastn(
                                    -query   => $seq,
                                    -outfile => "_blast_report",
        );
    
    # clean up the temporary files created...
    $factory->cleanup; 
    
    # process the result file using Bio::SearchIO...
    my $blast_report_in = Bio::SearchIO->new(
                                            -format => 'blast' ,
                                            -file   => 'blast_report',
        );
    
    # create a result object... 
    my $result = $blast_report_in->next_result;
    
    # find out if there is a hit...
    my $hits = $result->num_hits;
    
    
    if ($hits == 0) {
        
        # if not a hit move onto next...
        next;                       
    } 
    else {
        
        # else get the query name and print it to output file...
        my $query_name = $result->query_name;
        print $outfile "$query_name\n";
    }
}



