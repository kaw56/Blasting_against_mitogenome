#!usr/bin/perl
# annotate_hits.pl
# version 1.5 for blasting against mitochondrial genome
use warnings;
use strict;

use Bio::SearchIO;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlastPlus;

#####################
# setting things up #
#####################

die 'usage: annotate_hits.pl <query sequence> <query sequence format> <database for comparison> <definition file>' 
    unless (@ARGV == 3);

my $sequence_filename   = shift;
#my $format              = shift;
my $database_filename   = shift;
my $annotation_filename = shift;
my $report_filename     = "$sequence_filename-RESULTS.txt";

# check whether the file to search against exists...
die "Check that the database/file to be searched against exists!\n" 
    unless (-e $database_filename);

# open output for printing results...
open (my $report, ">", $report_filename) or die "Can't open $report_filename for writing, $!";

# create seqIO object...
my $seq_in  = Bio::SeqIO->new(
                            -format => 'fastq',
                            -file   => $sequence_filename,
    );

# create a blast factory, set up ep mito as subject...
my $factory = Bio::Tools::Run::StandAloneBlastPlus->new(
                            -db_data => $database_filename,
                            -create  => 1,
    ); 

# create a hash of annotation positions...
open (my $annotation, '<', $annotation_filename) or die "Can't open $annotation_filename, $!";

my %annotation_for = ();

while (my $line = <$annotation>) {
    
    # split line on space into an array...
    my @gene_information = split(" ", $line);
    
    # take the gene name as the key for the hash...
    my $gene_name = shift @gene_information;
    
    # reference to the 2 position values...   
    my $positions_ref = \@gene_information;
    
    # put into the hash...     
    $annotation_for{$gene_name} = $positions_ref;
    
}

close ($annotation);

#########
# blast #
#########

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
                                            -file   => '_blast_report',
        );
    
    # create a result object... 
    my $blast_result = $blast_report_in->next_result;
    
    # find out if there is a hit...
    my $number_of_hits = $blast_result->num_hits;
    
    if ($number_of_hits == 0) {
        
        # if not a hit move onto next...
        next;         
    } 
    else {
        
        # else get the query name...
        my $read_identifier = $blast_result->query_name;
        
        # look at the hits (even if there is more than one)
        while (my $blast_hit = $blast_result->next_hit) {
            
            while (my $hsp = $blast_hit->next_hsp) {
                # get start point and end point in an array
                my @range = $hsp->range('hit');
                
                # compare to the hash of annotations
                foreach my $gene (keys %annotation_for) {
    
                    # use hash of arrays to get range  
                    # and ask if at least one of the query numbers is in the range...
                    if ( ($range[0] >= ${$annotation_for{$gene}}[0] && $range[0] <= ${$annotation_for{$gene}}[1]) 
                            || ($range[1] >= ${$annotation_for{$gene}}[0] && $range[1] <= ${$annotation_for{$gene}}[1]) ) {
                       
                       # print to the report
                        print $report "$read_identifier $gene ", join( " ", @range), "\n";
                    }
                    else {
                        next;
                    }
                }
            }
        }
    }   
}

# remove temporary blast report
unlink "_blast_report";


