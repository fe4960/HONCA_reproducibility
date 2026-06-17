awk '{if(length($1)==7){print "0"$_}else{print $_}}' /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/5_refine_major/interval_to_time.txt_simple > /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/5_refine_major/interval_to_time.txt_simple_format

less /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/5_refine_major/Table_HCA_donorinfo.txt_format | cut -f 1,4,7 | sort -u > /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/5_refine_major/Table_HCA_donorinfo.txt_format_uniq
