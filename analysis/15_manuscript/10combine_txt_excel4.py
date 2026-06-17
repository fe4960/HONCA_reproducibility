import pandas as pd
import glob
import os

# Folder containing text files
input_folder = "/dfs3b/ruic20_lab/junw42/reference/GWAS/finemap_sumstat/suise/"

# Output Excel file
output_excel = "/dfs3b/ruic20_lab/junw42/reference/GWAS/finemap_sumstat/suise/combined_files.xlsx"

# Find all txt files
txt_files = glob.glob(os.path.join(input_folder, "*pip_var_all_merged_uniform_anno_ON"))

# Create Excel writer
with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
    
    for file in txt_files:
        
        # Read text file
        # Change sep="\t" if tab-delimited
        df = pd.read_csv(file, sep="\t", header=0, index_col=0 )
        
        # Sheet name from file name
        sheet_name = os.path.splitext(os.path.basename(file))[0]
        
        # Excel sheet names must be <=31 characters
        sheet_name = sheet_name[:31]
        
        # Write to Excel sheet
        df.to_excel(writer, sheet_name=sheet_name, index=True)

print(f"Saved to {output_excel}")
