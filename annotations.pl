#!usr/bin/perl
# annotations.pl
use warnings;
use strict;

die "usage annotations.pl <first integer> <second integer>\n" unless (@ARGV == 2);

open (my $infile, '<', 'mito_anno.txt') or die "Can't open mito_anno.txt, $!";

my %annotation_for  = ();
my @range           = @ARGV;

# read in annotation file line by line...
while (my $line = <$infile>) {
    
    # split line on space into an array...
    my @gene_information = split(" ", $line);
    
    # take the gene name as the key for the hash...
    my $gene_name = shift @gene_information;
    
    # reference to the 2 position values...   
    my $positions_ref = \@gene_information;
    
    # put into the hash...     
    $annotation_for{$gene_name} = $positions_ref;
    
}

# take 2 numbers and see if they are within the range

foreach my $gene (keys %annotation_for) {
    
    # use hash of arrays to get range  
    # and ask if at least one of the query numbers is in the range...
    if ( ($range[0] >= ${$annotation_for{$gene}}[0] && $range[0] <= ${$annotation_for{$gene}}[1]) 
            || ($range[1] >= ${$annotation_for{$gene}}[0] && $range[1] <= ${$annotation_for{$gene}}[1]) ) {
        print "$gene \n";
    }
    else {
        next;
    }
}


