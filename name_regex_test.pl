#!usr/bin/perl
# name_regex_test.pl
use warnings;
use strict;

my $gene_name = "lcl|nad1";

$gene_name =~ s/lcl\|(\w{4})/$1/;

print "$gene_name\n";
