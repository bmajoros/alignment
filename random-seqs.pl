#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <length> <num-seqs>\n" unless @ARGV==2;
my ($length,$n)=@ARGV;

for(my $i=0 ; $i<$n ; ++$i) {
  my $seq=generate($length);
  print ">$i\n$seq\n";
}

#=================================================

sub generate
  {
    my ($length)=@_;
    my $seq;
    for(my $i=0 ; $i<$length ; ++$i) {
      my $c;
      my $r=rand(1);
      if($r<.25) {$c='A'}
      elsif($r<0.50) {$c='C'}
      elsif($r<0.75) {$c='G'}
      else {$c='T'}
      $seq.=$c;
    }
    return $seq;
  }



