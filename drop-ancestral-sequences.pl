#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.fasta> <out.fasta>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

my $reader=new FastaReader($infile);
my $writer=new FastaWriter;
open(OUT,">$outfile") || die "can't open file for writing: $outfile\n";
while(1) {
  my ($defline,$sequence)=$reader->nextSequence();
  last unless defined($defline) && length($defline)>0;
  $defline=~/>(\S+)/ || die "Can't parse defline: $defline\n";
  next if($1=~/^A\d+$/ || $1 eq "root");
  $writer->addToFasta($defline,$sequence,\*OUT);
}
close(OUT);


