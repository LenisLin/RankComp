## RankCompV2 main 
source("./RankCompV2-functions.R")

## 1. Load test dataset
dataset_prefix="FOLFOX_GSE104645"
test_data=load_test_data(prefix=paste0("./dataset/",
                         dataset_prefix))

## split the data into two groups: sensitive and resistant
test_data=test_data_split(test_data)

## 2. calculate the rank of two genes and select the stable gene pairs in different groups
gene_rank_list=list()
#gene_rank_list[["sensitive"]]=calculate_gene_rank_v1(test_data[["sensitive"]])
#gene_rank_list[["resistant"]]=calculate_gene_rank_v1(test_data[["resistant"]])

gene_rank_list[["sensitive"]]=calculate_gene_rank_v2(test_data[["sensitive"]])
gene_rank_list[["resistant"]]=calculate_gene_rank_v2(test_data[["resistant"]])

## 3. select background gene pairs
background_genepairs=select_background_gene_paris(gene_rank_list)
rm(gene_rank_list)

## 4. for each genes calculate Fisher's exact test matrix
Fisher_test_matrix=form_Fisher_test_matrix(background_genepairs)

## 5. Fisher's exact test
Fisher_test_result=Fisher_exact_test(Fisher_test_matrix)

## 6. identify the up-regulated or down regulated
DEGs_result=identify_the_genes_regulated(Fisher_test_result,background_genepairs)

saveRDS(DEGs_result,file = paste0(dataset_prefix,"_DEGs.RDS"))
  