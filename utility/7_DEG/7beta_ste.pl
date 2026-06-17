#!/user/bin/perl -w
use strict;
#my @trait = ("raceAfrican", "raceAsian", "raceHispanic", "age", "gender", "ageGender");
my @trait=("age");
#my @cell = ("Rod", "Cone", "BC", "HC", "AC", "RGC", "MG"); #, "Astrocyte", "Microglia"); 
my $maincell=$ARGV[0];
my $subcell=$ARGV[1];
my @cell;
open(INPUT,$subcell);
while(my $line = <INPUT>){
	chomp $line;
	push(@cell,$line);
}
close(INPUT);
#my $main="/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/region_race_batch_ageNum_dream_atlas_sample";
my $main="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/".$maincell."_subclass_DEG/";
for my $tr ( @trait){
my %gene=();
for(my $i=0; $i<=$#cell; $i++){
my $dir=$main."/$cell[$i]"."/";

my $file = "$dir/$cell[$i]"."_DEG_".$tr;
print "$file";
#exit;
open(INPUT,$file);
<INPUT>;
while(my $line = <INPUT>){
chomp $line;
my @info = split(/\t/,$line);
#my $beta = 2**$info[1];
my $beta = $info[1];

my $ste= abs($beta/$info[3]);
$gene{$info[0]}->{$cell[$i]}->{"ste"}=$ste;
$gene{$info[0]}->{$cell[$i]}->{"beta"}=$beta;
$gene{$info[0]}->{$cell[$i]}->{"z"}=$info[-2];

}
}

my $outb = "$main/$tr"."_DEG_res_cpm_".$tr."_beta_ancFix";
my $outs = "$main/$tr"."_DEG_res_cpm_".$tr."_ste_ancFix";
my $outz = "$main/$tr"."_DEG_res_cpm_".$tr."_z_ancFix";


open(OUTPUTb,">$outb");
open(OUTPUTs,">$outs");
open(OUTPUTz,">$outz");

for(my $i=0; $i<=$#cell; $i++){
print OUTPUTb "\t$cell[$i]";
}
print OUTPUTb "\n";

for(my $i=0; $i<=$#cell; $i++){
print OUTPUTs "\t$cell[$i]";
}
print OUTPUTs "\n";

for(my $i=0; $i<=$#cell; $i++){
print OUTPUTz "\t$cell[$i]";
}
print OUTPUTz "\n";



for my $gene ( keys %gene){
print OUTPUTb "$gene";
print OUTPUTs "$gene";
print OUTPUTz "$gene";

for(my $i=0; $i<=$#cell; $i++){
if(!(defined $gene{$gene}->{$cell[$i]}->{"beta"})){
$gene{$gene}->{$cell[$i]}->{"beta"}=0;
#print "$cell[$i]\t$gene\t$tr\n";

}
if(!(defined $gene{$gene}->{$cell[$i]}->{"ste"})){
$gene{$gene}->{$cell[$i]}->{"ste"}=0;
}

if(!(defined $gene{$gene}->{$cell[$i]}->{"z"})){
$gene{$gene}->{$cell[$i]}->{"z"}=0;
}

print OUTPUTb "\t$gene{$gene}->{$cell[$i]}->{beta}";
print OUTPUTs "\t$gene{$gene}->{$cell[$i]}->{ste}";
print OUTPUTz "\t$gene{$gene}->{$cell[$i]}->{z}";

}
print OUTPUTb "\n";
print OUTPUTs "\n"; 
print OUTPUTz "\n"; 

}
}
