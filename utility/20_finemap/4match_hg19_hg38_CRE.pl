#!/user/bin/perl -w
use strict;
#my $file = "/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/g2p_cA_final_major_full_max_LCRE_uniq_clean";
#my $file = "/dfs3b/ruic20_lab/jinjingj/human_atlas_atac/finemap/data/g2p_cA_final_LCRE_uniq_clean";
#my $file = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/g2p_cA_final_major_full_max_LCRE_uniq_clean";
my $file = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/g2p_cA_final_major_full_max_LCRE_uniq_clean";
open(INPUT,$file);
<INPUT>;
my %hash;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\s+/,$line);
$hash{$info[1]}->{$info[0]}=1;
}

#my $hg38="/dfs3b/ruic20_lab/jinjingj/human_atlas_atac/finemap/data/g2p_cA_final_LCRE_uniq_clean_hg38_hg19_out";
my $hg38="/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/g2p_cA_final_major_full_max_LCRE_uniq_clean_hg38_hg19_out";
#my $hg38="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/g2p_cA_final_major_full_max_LCRE_uniq_clean_hg38_hg19_out";
my $output = $hg38."_CRE";
open(OUTPUT,">$output");
print OUTPUT "hg38\thg19\tgene\n";
open(INPUT1,$hg38);
<INPUT1>;
#hg38    hg19
#chr10:100000226-100000726       chr10:101759983-101760483
while(my $line = <INPUT1>){
chomp $line;
my @info = split(/\s+/,$line);
for my $gene ( keys %{$hash{$info[0]}}){
print OUTPUT "$info[0]\t$info[1]\t$gene\n";
}
}
