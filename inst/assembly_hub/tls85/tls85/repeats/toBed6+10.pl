#!/usr/bin/env perl

use strict;
use warnings;

my $sixteenSamples = 0;
my $fifteenSamples = 0;

my $argc = scalar(@ARGV);
if ($argc < 1) {
  printf STDERR "usage: toBed6+10.pl rmskClass/file.tab > ricCom1.class.tab\n";
  printf STDERR "e.g.: ./toBed6+10.pl rmskClass/rRNA.tab > ricCom1.rmsk.RNA.tab\n";
  exit 255;
}

while (my $file = shift) {

open (FH, "<$file") or die "can not read $file";
while (my $line = <FH>) {
  next if ($line =~ m/^\s+SW\s+perc.*/);
  next if ($line =~ m/^score\s+div.*/);
  next if ($line =~ m/^$/);
  chomp $line;
  my @a = split('\t', $line);
  shift @a;
  die "ERROR: not finding 16 fields in a line" if (scalar(@a) != 16);
  printf "%s\t%d\t%d\t%s\t0\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
    $a[4], $a[5], $a[6], $a[9], $a[8], $a[0],
    $a[1], $a[2], $a[3], $a[7], $a[10], $a[11], $a[12], $a[13], $a[14];
}
close (FH);
}
