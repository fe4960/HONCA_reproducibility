awk -F "," '{print $1","$2","$3}' /dfs3b/ruic20_lab/junw42/reference/cellphonedb/cellphonedb-data/data/gene_input.csv | sort -u > /dfs3b/ruic20_lab/junw42/reference/cellphonedb/cellphonedb-data/data/gene_input_uniq.csv

[junw42@login-i17:/dfs3b/ruic20_lab/junw42] $awk -F "," '{print $1}' /dfs3b/ruic20_lab/junw42/reference/cellphonedb/cellphonedb-data/data/gene_input_uniq.csv | sort | uniq -c | awk '{if($1>1){print}}'
      2 NRXN1
[junw42@login-i17:/dfs3b/ruic20_lab/junw42] $less /dfs3b/ruic20_lab/junw42/reference/cellphonedb/cellphonedb-data/data/gene_input_uniq.csv | grep NRXN1
NRXN1,P58400,NRXN1
NRXN1,Q9ULB1,NRXN1
