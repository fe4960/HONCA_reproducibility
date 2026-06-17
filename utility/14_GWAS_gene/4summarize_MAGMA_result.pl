#!/usr/bin/perl -w
#my $control_file = "/storage/chenlab/Users/junwang/human_meta/GWAS_file_list_new_control_MAGMA";
my $control_file = "/dfs3b/ruic20_lab/junw42/eye_QTL/scripts/2_GWAS/MAGMA/File.list1";
my %control;
my $s =$ARGV[0];
#print $s;
#exit;
open(INPUT,$control_file);
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
$control{$info[0]} = $info[1];
#print "$info[0]\n";;
}
#my $output = $control_file."_summary_list_new"."_$s";
my $output = $control_file."_summary_list_new"."_MAGMA_$s";

open(OUTPUT,">$output");
#print OUTPUT "FDR\tCell_type\tTrait\n";

print OUTPUT "FDR\tP-val\tCell_type\tTrait\n";
#my @list = ("MAGMA_list_major"); #Astrocyte"); #_Oligodendrocyte"); #,"MAGMA_list1");
#my @list = ("MAGMA_list_celltype"); #Astrocyte"); #_Oligodendrocyte"); #,"MAGMA_list1");
my @list= ("MAGMA_list_$s");

#my @list= ("MAGMA_list_Astrocyte");
#my @list = ("MAGMA_list"); #,"MAGMA_list1");
for my $list (@list){
#my $dir = "/storage/chen/home/jw29/human_meta/scripts/GWAS/$list";
my $dir="eye_QTL/scripts/2_GWAS/MAGMA/MAGMA_Files/$list";	
open(INPUT,$dir);
while(my $line = <INPUT>){
chomp $line;

#my $file = "$line/MAGMA_result_atlas"; #ad
my $file = "$line"; #."$s";
print "$file\n";
open(INPUT1,$file)|| die ("$file");
<INPUT1>;
while(my $line1 = <INPUT1>){
chomp $line1;
my @info1 = split(/\t/,$line1);
#$info1[1] =~ s/_MungeSumstats.txt.35UP.10DOWN//g;

#$info1[1] =~ s/_MungeSumstats.txt_reform.35UP.10DOWN//g;
#$info1[1] =~ s/_MungeSumstats.txt1_reform.35UP.10DOWN//g;
#$info1[1] =~ s/_MungeSumstats_rmFRQ.txt_reform.35UP.10DOWN//g;
#$info1[1] =~ s/_MungeSumstats.txt_rmFQS_reform.35UP.10DOWN//g;
#if($info1[-3] eq "Linear"){
if($info1[-3] eq "retina_linear.Linear"){

		print OUTPUT  "$info1[-2]\t$info1[8]\t$info1[-1]\t$control{$info1[1]}\n";
	#		print OUTPUT  "$info1[-2]\t$info1[-1]\t$control{$info1[1]}\n";

#print OUTPUT "$info1[7]\t$info1[-1]\t$control{$info1[1]}\n";
}
}
}
}
