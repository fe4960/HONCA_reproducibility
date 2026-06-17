#!/usr/bin/perl -w
use strict;
#HCA_ON/data/7_DEG/major_majorclass_DEG/Oligodendrocyte_precursor_cell/Oligodendrocyte_precursor_cell_DEG_age_cpm1
#
my @dem=("age_cpm1"); #,"raceEur","raceAfrican");
#my @cell=("major"); #"all","Oligodendrocyte","Endothelial_cell","Astrocyte","Fibroblast");
#my @cell = ("Oligodendrocyte");
my @cell = ("Astrocyte");
for my $dem (@dem){
for my $ct (@cell){
	#my $subfile="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/ct_list_"."$ct";
my $subfile="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/ct_list_"."$ct";

#####my $dir = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/Oligodendrocyte_subclass_DEG/"; #$lines"; #"/Artery_Arteriole_DEG_age";
my $dir = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/Astrocyte_subclass_DEG/"; #$lines"; #"/Artery_Arteriole_DEG_age";

#my $dir = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/major_majorclass_DEG/"; #$lines"; #"/Artery_Arteriole_DEG_age";
my %hash=();
my %hash_up=();
my %hash_down=();


open(INPUTs,"$subfile");
while(my $lines = <INPUTs>){
chomp $lines;
my $file = $dir.$lines."/$lines"."_DEG_age_cpm1";

#my $dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/".$ct."_subclass_DEG/";

open(INPUT,$file);
<INPUT>;
#my $cell=<INPUT>;
#chomp $cell;
#my @cell = split(/\s+/,$cell);
#my $n=0;
while(my $line = <INPUT>){
chomp $line;
#$n++;
my @info = split(/\s+/,$line);
#for(my $i=0; $i<=$#info; $i++){
if($info[-1] < 0.05){

#if($info[5] < 0.05){
$hash{$info[0]}->{$lines}=1;
}
#}
}
}
my %sum=();
for my $gene ( keys %hash){
my @cell1= sort keys %{$hash{$gene}};
my $cell_line = join("&",@cell1);
$sum{$cell_line}->{$gene}=1;
}


#my $dir1 = "/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/region_race_batch_ageNum_dream_atlas_sample";
my $output = "$dir/".$ct."_upset_list_$dem"."_all_ct_nomashr";
my $output1 = "$dir/".$ct."_upset_list_$dem"."_all_ct_nomashr_tab";

my %count=();

open(OUTPUT,">$output");
open(OUTPUT1,">$output1");

my $line="";
for my $key ( keys %sum){
my $gene_num = keys %{$sum{$key}};
if($gene_num >=10){
$line .= "\'$key\'=$gene_num,";
print OUTPUT1 "$key\t$gene_num\n";
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
