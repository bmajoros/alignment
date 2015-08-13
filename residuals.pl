#!/usr/bin/perl
use strict;
use SummaryStats;

my @residuals;
for(my $i=0 ; $i<10000 ; ++$i) {
  open(IN,"test-mixture-model .2 10 $i |");
  my $listen=0;
  while(<IN>) {
    if(/residuals/) {$listen=1;next}
    if($listen) {
      chomp;
      my @r=split/\s+/,$_;
      foreach my $x (@r) {
	push @residuals,abs($x);
      }
    }
  }
  close(IN);
}


my ($mean,$stddev,$min,$max)=SummaryStats::summaryStats(\@residuals);
print "$mean +/ $stddev ($min - $max)\n";

