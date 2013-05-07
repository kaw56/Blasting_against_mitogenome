#!/usr/bin/perl
# blast_against_ep_mito.pl
use strict;
use warnings;

# provide the query file on command line 
my $query = shift @ARGV;

# run a blast search using the specified file against the Ep mitochondrial genome
my $align = `blastn -subject /Users/oj/Documents/Katherine/Sequence/Fasta_files/partial_Ep_mitogenome -query $query`;

# Print results to a text file
my $filename = "$query aligned.txt";
open (OUT, '>', $filename);
print OUT $align;
close (OUT);
