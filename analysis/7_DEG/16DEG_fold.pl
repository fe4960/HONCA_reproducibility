#!/usr/bin/perl -w
use strict;

my @dem=("age"); #,"raceEur","raceAfrican");
my @cell=("all","Oligodendrocyte","Endothelial_cell","Astrocyte","Fibroblast");
for my $dem (@dem){
for my $ct (@cell){
my $subfile="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/ct_list_"."$ct";

my $dir = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/$ct"."_subclass_DEG/"; #$lines"; #"/Artery_Arteriole_DEG_age";
my $output = $dir.$ct."_DEG_p_fc_padj0.05";
my %hash=();
my %hash_up=();
my %hash_down=();
open(OUTPUT,">$output");

open(INPUTs,"$subfile");
while(my $lines = <INPUTs>){
chomp $lines;
my $file = $dir.$lines."/$lines"."_DEG_age";

open(INPUT,$file);
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
#if($info[5] < 0.05){
$hash{$info[0]}->{$lines}=1;
print OUTPUT "$info[0]\t$lines\t$info[1]\t$info[4]\t$info[5]\n";
#}
}
}
}
}

