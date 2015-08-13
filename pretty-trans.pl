#!/usr/bin/perl
use strict;

while(<STDIN>) {
  chomp;
  my @fields=split/\s+/,$_;
  foreach my $field (@fields) {
    if($field eq "-inf") {$field="0.000"}
    else {
      $field=exp($field);
      if($field=~/^0\.(...)/) {
	my $x=$1;
	$field="0.$x";
      }
      elsif($field=~/(.*)\.(.*)e-(.*)/) {
	my ($x,$y,$z)=($1,$2,$3);
	$field="0." . "0"x($3-1) . $x . $y;
	if($field=~/^0\.(...)/) {
	  my $x=$1;
	  $field="0.$x";
	}
      }
      elsif($field eq "1") {$field="1.000"}
      else {die $field}
      if($field<0.001) {$field=0.001}
    }
    print "$field\t";
  }
  print "\n";
}

