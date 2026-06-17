#!/usr/bin/perl -w
use strict;
my $file = $ARGV[0]."_S1_L001_".$ARGV[1]."_001.fastq.gz";
my $out= $ARGV[0]."_rename_S1_L001_".$ARGV[1]."_001.fastq";
open(OUTPUT,">$out");
open(INPUT,"zcat $file|");
while(my $line = <INPUT>){
chomp $line;
my @info =split(/\s+/,$line);
my @info1 = split(/\./,$info[0]);
$info[0]=$info1[0].".".$info1[1];
my $new_line=join(" ", @info);
print OUTPUT "$new_line\n";
$line = <INPUT>;
print OUTPUT "$line";
}
`gzip $out`;
