# load functions and packages
import sys
#sys.path.append('/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/from_Tingting/')
sys.path.append('/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/11_spatial/')

import get_transcript_cor

# Inputs
#xenium_bundle = "/dfs3b/ruic20_lab/rawdata/Xenium/20241018__Human_Eye_MDF_FFPE_lung_probe/FFPE/output-XETG00221__0029326__24-0550_OS_Peripheral__20241018__000131/"
xenium_bundle = "/dfs3b/ruic20_lab/rawdata/Xenium/20250620__203138__Human_70y_ASRetina_20250620/output-XETG00221__0068480__25-019211OS_PP_02__20250620__203217/"

#xenium_bundle = "/dfs3b/ruic20_lab/rawdata/Xenium/20250620__203138__Human_70y_ASRetina_20250620/output-XETG00221__0068480__25-019211OS_PP_01__20250620__203217/"
#polygon_file = "/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/24-0550_OS_Peripheral_retina.csv"
polygon_file = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/output-XETG00221__0068480__25-019211OS_PP_02__20250620__203217/Selection_1_coordinates.csv"

#polygon_file = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/output-XETG00221__0068480__25-019211OS_PP_01__20250620__203217/Selection_2_coordinates.csv"
#sudo_bulk_PATH = "/dfs3b/ruic20_lab/tingty7/projects/humen_eye_xenium/retina_sanes_H3_gene_total_counts.csv"
sudo_bulk_PATH="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/ONONH_all_rawcount_rmRet_rawgene_count_ON.csv"
#sudo_bulk_PATH="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/ONONH_all_rawcount_rmRet_rawgene_count.csv"

# get all the transcripts in the selected area(s)
transcripts_in_polygons, counts_df = get_transcript_cor.get_transcripts_in_polygon(xenium_bundle, polygon_file)

counts_df.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/output-XETG00221__0068480__25-019211OS_PP_02__20250620__203217/Selection_1_coordinates_ct.csv")
# get transcript density in the selected area(s)
get_transcript_cor.get_transcripts_density_in_polygon(transcripts_in_polygons, polygon_file)

# get transcripts level correlation
#get_transcript_cor.get_transcripts_cor_polygon(counts_df, sudo_bulk_PATH, savepath="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/output-XETG00221__0068480__25-019211OS_PP_01__20250620__203217/Selection_2_cor.pdf")
#get_transcript_cor.get_transcripts_cor_polygon(counts_df, sudo_bulk_PATH, savepath="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/output-XETG00221__0068480__25-019211OS_PP_01__20250620__203217/Selection_2_cor_ON.pdf")
#get_transcript_cor.get_transcripts_cor_polygon(counts_df, sudo_bulk_PATH, savepath="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/output-XETG00221__0068480__25-019211OS_PP_01__20250620__203217/Selection_1_cor_ON.pdf")
#get_transcript_cor.get_transcripts_cor_polygon(counts_df, sudo_bulk_PATH, savepath="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/output-XETG00221__0068480__25-019211OS_PP_01__20250620__203217/Selection_1_cor.pdf")

get_transcript_cor.get_transcripts_cor_polygon(counts_df, sudo_bulk_PATH, savepath="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/output-XETG00221__0068480__25-019211OS_PP_02__20250620__203217/Selection_1_cor_ON.pdf")

