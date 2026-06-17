#!/bin/sh
source ~/.condainit
conda activate escape
R < /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/18_escape/1_opc.R --no-save > /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/18_escape/1_opc.out 2> /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/18_escape/1_opc.err
