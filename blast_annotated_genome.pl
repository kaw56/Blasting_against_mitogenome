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
    if ($seq->length == 0) {
        next;
    } else {

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
        } else {
            # get all the things to go into both reports in a hash...
            $information_for{'query_name'}       = $blast_result->query_name;
            $information_for{'hit_num'}          = $blast_result->num_hits();


            while (my $blast_hit = $blast_result->next_hit() ) {
                # add hit specific things to hash
                $information_for{'hit_name'}    = $blast_hit->name;
                $information_for{'hit_sig'}     = $blast_hit->significance;
                $information_for{'num_hsp'}     = $blast_hit->num_hsps;

                # adjust the hit name
                $information_for{'hit_name'} =~ s/lcl\|(\w{4})/$1/;


                while(my $hsp = $blast_hit->next_hsp() ) {
                    # add hsp stuff to the hash
                    $information_for{'start'}           = $hsp->start('hit');
                    $information_for{'end'}             = $hsp->end('hit');
                    $information_for{'score'}           = $hsp->score;
                    $information_for{'length'}          = $hsp->length('total');

                    # print to summary file
                    print $summary "$information_for{'query_name'} $information_for{'hit_name'} $information_for{'score'} $information_for{'length'} $information_for{'start'} $information_for{'end'} \n";

                }
            }
        }
    }
}


$factory->cleanup;

# remove temporary blast report
unlink "_blast_report";


