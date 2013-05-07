#!usr/bin/perl
# compare.pl
use warnings;
use strict;

die "usage compare <first integer> <second integer>\n" unless (@ARGV == 2);

my @range  = @ARGV;
my @master = (1, 10);

if ( !($range[0] >= $master[0] && $range[0] <= $master[1]) && !($range[1] >= $master[0] && $range[1] <= $master[1]) ) {
    print "Boo\n";
    
} 
else { 
    print "yay\n";
}
