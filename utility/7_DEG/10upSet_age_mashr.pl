#!/usr/bin/perl -w
use strict;

my @dem=("age"); #,"raceEur","raceAfrican");
my @cell=("all");
#my @cell = ("Oligodendrocyte");
#my @cell=("Endothelial_cell","Astrocyte","Fibroblast");
for my $dem (@dem){
for my $ct (@cell){
my $dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/".$ct."_subclass_DEG/";
my %hash=();
my %hash_up=();
my %hash_down=();
my $file = $dir.$dem."_DEG_res_cpm_".$dem."_lfsr_ct_ancFix";

open(INPUT,$file);
my $cell=<INPUT>;
chomp $cell;
my @cell = split(/\s+/,$cell);
my $n=0;
while(my $line = <INPUT>){
chomp $line;
$n++;
my @info = split(/\s+/,$line);
for(my $i=0; $i<=$#info; $i++){
if($info[$i] < 0.01){
$hash{$n}->{$cell[$i]}=1;
}
}
}
my %sum=();
for my $gene ( keys %hash){
my @cell1= sort keys %{$hash{$gene}};
my $cell_line = join("&",@cell1);
$sum{$cell_line}->{$gene}=1;
}


#my $dir1 = "/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/region_race_batch_ageNum_dream_atlas_sample";
my $output = "$dir/".$ct."_upset_list_$dem"."_all_ct";

my %count=();

open(OUTPUT,">$output");
my $line="";
for my $key ( keys %sum){
my $gene_num = keys %{$sum{$key}};
if($gene_num >=10){
$line .= "\'$key\'=$gene_num,";
my @key=split("&",$key);
if($#key==0){
$count{$key}->{only}=$gene_num;
$count{$key}->{all}+=$gene_num;
}else{
for my $key1 ( @key){
$count{$key1}->{all}+=$gene_num;
}
}
}
}


my $new_line = substr($line,0,length($line)-1);
print OUTPUT "$new_line";
}
}
