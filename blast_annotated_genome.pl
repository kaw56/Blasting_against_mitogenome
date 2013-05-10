#!usr/bin/perl
# blast_annotated_genome.pl
# version 2 for blasting against mitochondrial genome
# this time using a database of mitochondrial annotations
use warnings;
use strict;

use Bio::SearchIO;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlastPlus;


##########
# set up #
##########

die 'usage: annotate_hits.pl <query> <query format> <database for comparison>' 
    unless (@ARGV == 3);

# filenames    
my $sequence_filename   = shift;
my $format              = shift;
my $database_filename   = shift;
my $report_filename     = $sequence_filename . "RESULTS_summary.txt";

die "Check that the file to be searched against is in this directory\n" 
    unless (-e $database_filename);

# open output for printing results...
open (my $summary, ">", $report_filename) or die "Can't open $report_filename for writing, $!";
        
# create seqIO object...
my $seq_in  = Bio::SeqIO->new(
                            -format => $format,
                            -file   => $sequence_filename,
    );

# create a blast factory, set up database as subject...
my $factory = Bio::Tools::Run::StandAloneBlastPlus->new(
                            -db_data => $database_filename,
                            -create  => 1
    ); 

#########
# blast # 
#########

# initialise reporting hash
my %information_for = ();

# run a blast on all the query sequences in that file
while (my $seq = $seq_in->next_seq() ) {
       
    # blast against ep mito and print result to temp file...
    my $blast_report = $factory->blastn(
                                    -query   => $seq,
                                    -outfile => "_blast_report",
        );
    
    # clean up the temporary files created by blast plus...
    $factory->cleanup; 
    
    # process the result file using Bio::SearchIO...
    my $blast_report_in = Bio::SearchIO->new(
                                            -format => 'blast' ,
                                            -file   => '_blast_report',
        );
    
    # create a result object... 
    my $blast_result = $blast_report_in->next_result();
    
    # find out if there is a hit...
    my $number_of_hits = $blast_result->num_hits();
    
    if ($number_of_hits == 0) {
        # if not a hit move onto next query sequence...
        next;         
    } 
    else {
        # get all the things to go into both reports in a hash...
        $information_for{'query_name'}       = $blast_result->query_name;
        $information_for{'algorithm'}        = $blast_result->algorithm; 
        $information_for{'database_entries'} = $blast_result->database_entries;
        $information_for{'hit_num'}          = $blast_result->num_hits();
        
        # open an individual report file for output
        my $query_id = $information_for{'query_name'} . "_REPORT.txt";
        open(my $ind_report_out, '>', $query_id) 
                or die "can't open $query_id for writing, $!"; 
        
        # print the file header... 
        my $information_for_ref = \%information_for; # passing the hash by ref
        print_header($ind_report_out, $information_for_ref); 
        
        while (my $blast_hit = $blast_result->next_hit() ) {
            # add hit specific things to hash
            $information_for{'hit_name'}    = $blast_hit->name;
            $information_for{'hit_length'}  = $blast_hit->length;
            $information_for{'hit_sig'}     = $blast_hit->significance;
            $information_for{'num_hsp'}     = $blast_hit->num_hsps;
            
            # print hit related information
            $information_for_ref = \%information_for;
            print_hit($ind_report_out, $information_for_ref);
            
            while(my $hsp = $blast_hit->next_hsp() ) {
                # add hsp stuff to the hash
                $information_for{'evalue'}          = $hsp->evalue;
                $information_for{'frac_identical'}  = $hsp->frac_identical;
                $information_for{'alignment'}       = $hsp->get_aln;
                $information_for{'start'}           = $hsp->start('hit');
                $information_for{'end'}             = $hsp->end('hit');
                $information_for{'query_string'}    = $hsp->query_string;
                $information_for{'hit_string'}      = $hsp->hit_string;
                $information_for{'homology_string'} = $hsp->homology_string;
                
                if ($information_for{'evalue'} != 0) {
                    next;
                } else {
                    # print to the individual report
                    $information_for_ref = \%information_for;
                    print_hsp($ind_report_out, $information_for_ref);
                }
                              
                # print to summary file
                print $summary "$information_for{'query_name'} $information_for{'hit_name'} $information_for{'evalue'}\n"
                
            }
        }    
        
          
        
      
    
    
    
    }  
}

# remove temporary blast report
unlink "_blast_report";






#######################################
#######################################
##Subroutines
##

#######################################
# print_header
#
# prints the result level information to individual output file 
# print_header( <output file>, <reference to hash>)

sub print_header {
    
    my $outfile  = shift @_;
    my $info_ref = shift @_;
    
    my $header = <<EOF;
$info_ref->{'algorithm'} report for $info_ref->{'query_name'}

Database entries:   $info_ref->{'database_entries'}
Number of hits:     $info_ref->{'hit_num'}

EOF
    
    print $outfile $header;
    
}


#######################################
# print_hit
#
# prints the hit level information to output file
# print_hit( <output file>, <reference to hash>)
sub print_hit {
    my $outfile  = shift @_;
    my $info_ref = shift @_;
    
    my $hit = <<EOF;
HITS

Hit name:           $info_ref->{'hit_name'}
Length:             $info_ref->{'hit_length'}
Hit significance:   $info_ref->{'hit_sig'}
Number of HSPs:     $info_ref->{'num_hsp'}

EOF

    print $outfile $hit;
}
    
#######################################
# print_hsp
#
# prints hsp level information to individual output file
# print_hsp( <output file>, <reference to hash>)

sub print_hsp {
    my $outfile  = shift @_;
    my $info_ref = shift @_;
    
    my $hsp = <<EOF;
HSP
    evalue:                 $info_ref->{'evalue'}
    Fraction identical:     $info_ref->{'frac_identical'}
    Start:                  $info_ref->{'start'}
    End                     $info_ref->{'end'}
    
Query:  $info_ref->{query_string}
        $info_ref->{homology_string}
Sbjct:  $info_ref->{hit_string}


EOF

    print $outfile $hsp;
}



