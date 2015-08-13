#!/usr/bin/perl
use strict;
use ProgramName;
use TempFilename;

my $MIN_VALUE=20;

my $name=ProgramName::get();
die "$name <infile> <outfile>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

my $tmp=TempFilename::generate();
$tmp.=".txt";

$tmp="temptemptemp.txt";

my (%array,$largestX,$largestY);
open(IN,$infile) || die "can't open $infile";
while(<IN>) {
  my @fields=split/\s+/,$_;
  next unless @fields>=4;
  my ($x,$y,$state,$value)=@fields;
  if(!defined($array{$x}->{$y}) || $value>$array{$x}->{$y}) {
    $array{$x}->{$y}=$value;
  }
  if($x>$largestX) {$largestX=$x}
  if($y>$largestY) {$largestY=$y}
}
close(IN);

open(OUT,">$tmp") || die "can't create file: $tmp\n";
print OUT "# ImageMagick pixel enumeration: $largestX,$largestY,255,RGB\n";
for(my $x=0 ; $x<=$largestX ; ++$x) {
  for(my $y=0 ; $y<$largestY ; ++$y) {
    my $p=$array{$x}->{$y};
    my $z=defined($p) ? int($p*255) : 0;
    if(defined($p) && $z<$MIN_VALUE) {$z=$MIN_VALUE}
    $z=255-$z;
    print OUT "$x,$y: ( $z, $z, $z)\n";
  }
}
close(OUT);

system("convert $tmp $outfile");

