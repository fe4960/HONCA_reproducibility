ll /dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE267301/*_4.fastq.gz | awk -F "/" '{print $NF}' | sed -e "s/_4.fastq.gz//g" > /dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE267301_list
