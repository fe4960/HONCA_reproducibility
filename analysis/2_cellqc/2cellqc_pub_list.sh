#!/bin/sh

list="/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/cellranger/cellranger_list"
main="/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/cellranger/"
output=${list}_full
echo "sample	cellranger" >> $output

while read line
do
echo "$line	${main}/${line}/${line}/outs" >> $output
done < $list
