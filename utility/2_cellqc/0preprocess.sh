#!/bin/sh

#dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/"

#mkdir ${dir}

#for list1 in file_list_Hann 
#do

ls /dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/ > /dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/cellranger_list

list="/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/cellranger_list"

#list="/storage/chenlab/Data/HCA/human_ancestry_ret/RNA/cellranger/${list1}_dir"
#list=/storage/chenlab/Data/HCA/human_ancestry_ret/RNA/cellranger/file_list_cellqc_all_rmHRCA_list
while read line
do
#arr=($line)
arr=/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/${line}
mkdir ${arr}/filtered_feature_bc_matrix
cd ${arr}/filtered_feature_bc_matrix
tar -xvf ${arr}/filtered_feature_bc_matrix.tar.gz

mkdir ${arr}/raw_feature_bc_matrix
cd ${arr}/raw_feature_bc_matrix
tar -xvf ${arr}/raw_feature_bc_matrix.tar.gz

mkdir ${arr}/analysis
cd ${arr}/analysis
tar -xvf ${arr}/analysis.tar.gz
mv clustering/gene_expression_graphclust     clustering/graphclust
mv clustering/gene_expression_kmeans_2_clusters  clustering/kmeans_2_clusters
mv clustering/gene_expression_kmeans_3_clusters  clustering/kmeans_3_clusters
mv clustering/gene_expression_kmeans_4_clusters  clustering/kmeans_4_clusters
mv clustering/gene_expression_kmeans_5_clusters  clustering/kmeans_5_clusters
mv clustering/gene_expression_kmeans_6_clusters  clustering/kmeans_6_clusters
mv clustering/gene_expression_kmeans_7_clusters  clustering/kmeans_7_clusters
mv clustering/gene_expression_kmeans_8_clusters  clustering/kmeans_8_clusters
mv clustering/gene_expression_kmeans_9_clusters  clustering/kmeans_9_clusters
mv clustering/gene_expression_kmeans_10_clusters  clustering/kmeans_10_clusters

mv diffexp/gene_expression_graphclust  diffexp/graphclust
mv diffexp/gene_expression_kmeans_2_clusters diffexp/kmeans_2_clusters
mv diffexp/gene_expression_kmeans_3_clusters diffexp/kmeans_3_clusters
mv diffexp/gene_expression_kmeans_4_clusters diffexp/kmeans_4_clusters
mv diffexp/gene_expression_kmeans_5_clusters diffexp/kmeans_5_clusters
mv diffexp/gene_expression_kmeans_6_clusters diffexp/kmeans_6_clusters
mv diffexp/gene_expression_kmeans_7_clusters diffexp/kmeans_7_clusters
mv diffexp/gene_expression_kmeans_8_clusters diffexp/kmeans_8_clusters
mv diffexp/gene_expression_kmeans_9_clusters diffexp/kmeans_9_clusters
mv diffexp/gene_expression_kmeans_10_clusters diffexp/kmeans_10_clusters
mv pca/gene_expression_10_components  pca/10_components
mv tsne/gene_expression_2_components tsne/2_components
mv umap/gene_expression_2_components umap/2_components
done < ${list} 

#done
