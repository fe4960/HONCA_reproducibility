#!/user/bin/perl -w
use strict;

my $seq = uc("gctgagctaactggaAttgactgctcaatgga");
my $seq1 = uc("gctgagctaactggaGttgactgctcaatgga");


#my $seq = uc("gCcttcgtccggccgcagagg");
#my $seq1 = uc("gCcttcgtccagccgcagagg");
#`mkdir sc_human_retina/scripts/ATAC_data/`;
my $dir = "HCA_ON/data/13_GWAS/motif/";
`mkdir $dir`;
my @nt = ("A","C","G","T");

my $output = $dir."ref_SOX9_rs7636574";
my $output1 = $dir."alt_SOX9_rs7636574";

open(OUTPUT,">$output");
open(OUTPUT1,">$output1");

for(my $i=1; $i<=length($seq); $i++){
print OUTPUT "\t$i";
print OUTPUT1 "\t$i";
}
my @info = split(//,$seq);
my @info1 = split(//,$seq1);
print OUTPUT "\n";
print OUTPUT1 "\n";
for(my $i=0; $i<=$#nt; $i++){
print OUTPUT "$nt[$i]\t";
print OUTPUT1 "$nt[$i]\t";
for(my $j=0; $j<=$#info; $j++){
if($info[$j] eq $nt[$i]){
print OUTPUT "0.7\t";
}else{
print OUTPUT "0\t";
}

if($info1[$j] eq $nt[$i]){
print OUTPUT1 "0.7\t";
}else{
print OUTPUT1 "0\t";
}


}
print OUTPUT "\n";
print OUTPUT1 "\n";
}

exit;


