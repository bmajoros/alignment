#!/usr/bin/perl
use strict;

my $x=0;
my $y;
while(<STDIN>) {
    if(/likelihood=(\S+)/) {
	$y=$1;
	print "$x $y\n";
	++$x;
    }
    elsif(/a score=(\S+)/) {
	$y=$1;
	next unless $y<0;
	print "$x $y\n";
	++$x;
    }
    elsif(/REJECT/) {
	print "$x $y\n";
	++$x;
    }
}
