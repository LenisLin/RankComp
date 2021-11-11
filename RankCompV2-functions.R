## RankCompV2 functions

## 1. load test data
## the test data contain two parts: test_exp is the expression matrix and test_clinical. The rows of 
## test_exp are genes' symbol and columns are samples. Test_clinical's rows are samples and columns
## are the information of each sample
load_test_data=function(prefix){
  test_exp=readRDS(paste0(prefix,"_exp.rds"))
  test_clinical=readRDS(paste0(prefix,"_clinical.rds"))

  return_list=list()
  return_list[[1]]=test_exp
  return_list[[2]]=test_clinical
  names(return_list)=c("exp","clinical")
  
  return(return_list)
}

## split the data into two groups: sensitive and resistant
test_data_split=function(dataset_list){
  exp=dataset_list[["exp"]]
  clinical=dataset_list[["clinical"]]
  
  sensitive_index=grep(clinical$Drug_Response,pattern = "sensitive")
  resistant_index=grep(clinical$Drug_Response,pattern = "resistant")
  
  sensitive_exp=exp[,sensitive_index]
  resistant_exp=exp[,resistant_index]
  
  return_list=list()
  return_list[["sensitive"]]=sensitive_exp
  return_list[["resistant"]]=resistant_exp
  
  return(return_list)
}

## 2. calculate the rank of two genes
calculate_gene_rank_v1=function(exp){
  ## to save the working memory, using the dynamic planning to solve every matrix
  ## for the first sample
  result_exp=matrix(data = NA,nrow = nrow(exp),ncol = nrow(exp))
  rownames(result_exp)=rownames(exp)
  colnames(result_exp)=rownames(exp)
  i=1
  for (j in 1:nrow(exp)) {
    tem_vector=rep(exp[j,i],(nrow(exp)-j+1))
    result_exp[j,j:nrow(exp)]=tem_vector-exp[j:nrow(exp),i]
  }
  result_exp=ifelse(result_exp>0,1,0)
  
  for(i in 2:ncol(exp)){
    tem_matrix=matrix(data = NA,nrow = nrow(exp),ncol = nrow(exp))
    rownames(tem_matrix)=rownames(exp)
    colnames(tem_matrix)=rownames(exp)
    for (j in 1:nrow(exp)) {
      tem_vector=rep(exp[j,i],(nrow(exp)-j+1))
      tem_matrix[j,j:nrow(exp)]=tem_vector-exp[j:nrow(exp),i]
    }
    tem_matrix=ifelse(tem_matrix>0,1,0)
    
    ## find the same gene expression order gene pairs
    result_exp=find_same_REO_gene_pairs(result_exp,tem_matrix)
    cat("Iteration",i,"\n")
    cat("The result matrix is: \n")
    cat(result_exp[1:5,1:5],'\n')
  }
  
  ## reduce the result_exp as a 2-dimentional matrix
  result_exp[1:5,1:5]
  candidatic_gene_pairs=select_the_gene_pairs(result_exp)
  
  return(candidatic_gene_pairs)
}

## find the same gene expression order gene pairs
find_same_REO_gene_pairs=function(result_exp,tem_matrix){
  
  ## initiate a matrix to save the index
  mapping_matrix=matrix(data = NA,nrow = nrow(result_exp),ncol = ncol(result_exp))
  rownames(mapping_matrix)=rownames(result_exp)
  colnames(mapping_matrix)=rownames(result_exp)
  
  ## using the nxor to calculate the reverse order gene pairs
  mapping_matrix=xor(result_exp,tem_matrix)
  mapping_matrix[1:6,1:6]
  
  ## TURE in mapping_matrix represents the reverse order in a group which will be replaced -1
  result_exp_test=result_exp
  for(k in 1:ncol(result_exp)){
    result_exp_test[,k][mapping_matrix[,k]]=(-1)
  }
  return(result_exp_test)
}

## reduce the result_exp as a list
select_the_gene_pairs=function(result_exp){
  gene_pairs=list()
  for (i_ in 1:ncol(result_exp)){
    tem=result_exp[,i_]
    zero_index=which(tem==0)
    one_index=which(tem==1)
    
    tem=tem[c(zero_index,one_index)]
    gene_pairs[[i_]]=tem
  }
  names(gene_pairs)=colnames(result_exp)
  
  return(gene_pairs)
}

## 2. plus
calculate_gene_rank_v2=function(exp,significant_threshold=0.05,method="t_test"){
  ## t_test may be the speed-limited step
  if(method=="t_test"){
    result_list=list()
    for(i in 1:(nrow(exp)-1)){
      ## using a matrix to save the rank order
      tem_matrix=matrix(nrow = 1)
      k=2
      for(j in (i+1):nrow(exp)){
        gene_pair_tem=as.numeric(exp[i,])-as.numeric(exp[j,])
        gene_pair_rank_tem=ifelse(gene_pair_tem>0,1,0)
        
        p_value=t.test(x=gene_pair_tem,
                           y=rep(0,ncol(exp)))$p.value

        if(p_value>significant_threshold){
          next;
        }
        
        else{
          tem_matrix=cbind(tem_matrix,
                           matrix(data = ifelse(length(which(gene_pair_rank_tem==1))>length(which(gene_pair_rank_tem==0))
                                                ,1,0),nrow = 1,ncol = 1))
          colnames(tem_matrix)[k]=rownames(exp)[j]
          k=k+1
        }
      }
      
      tem_matrix=tem_matrix[,-1]
      cat('Gene ',i,' ',ncol(tem_matrix),'\n')
      result_list[[i]]=tem_matrix
    }
    names(result_list)=rownames(exp)[1:(nrow(exp)-1)]
    return(result_list)
  }
}

## 3. select background gene pairs
select_background_gene_paris=function(gene_rank_list){
  for(i in 1:length(gene_rank_list[[1]])){
    tem_s=gene_rank_list[["sensitive"]][[i]]
    tem_r=gene_rank_list[["resistant"]][[i]]
    
    stable_pairs=intersect(names(tem_s),names(tem_r))
    
    gene_rank_list[["sensitive"]][[i]]=gene_rank_list[["sensitive"]][[i]][stable_pairs]
    gene_rank_list[["resistant"]][[i]]=gene_rank_list[["resistant"]][[i]][stable_pairs]
  }
  return(gene_rank_list)
}

## 4. for each genes calculate Fisher's exact test matrix
form_Fisher_test_matrix=function(background_genepairs){

  ## form a list to save all the genes fisher_matrix
  resu_list=list()
  ## each gene
  for(i in 1:length(background_genepairs[[1]])){
    sensitive_tem=background_genepairs[[1]][[i]]
    resistant_tem=background_genepairs[[2]][[i]]
    
    G1_higher_Gi_index=which(sensitive_tem==1)
    G1_minor_Gj_index=which(sensitive_tem==0)
    
    ## in the control group(sensitive)
    a=length(G1_higher_Gi_index)
    b=length(G1_minor_Gj_index)
    
    ## in the disease group(resistant)
    x=length(which(resistant_tem[G1_higher_Gi_index]==0))
    y=length(which(resistant_tem[G1_minor_Gj_index]==1))
    
    mat_=matrix(data = c(a,a-x+y,b,b-y+x),nrow = 2,ncol = 2)
    rownames(mat_)=c("sensitive","resistant")
    colnames(mat_)=c("G1>Gi","G1<Gj")
    
    resu_list[[i]]=mat_
  }
  names(resu_list)=names(background_genepairs[[1]])
  return(resu_list)
  
}

## 5. Fisher's exact test
Fisher_exact_test=function(Fisher_test_matrix,FDR_threshold=0.1){
  result_matrix=matrix(ncol = 2,nrow = length(Fisher_test_matrix))
  colnames(result_matrix)=c("gene","p_value")
  for(i in 1:length(Fisher_test_matrix)){
    tem_mat=Fisher_test_matrix[[i]]
    
    test_resu=fisher.test(tem_mat)
    result_matrix[i,]=c(names(Fisher_test_matrix)[i],test_resu$p.value)
  }
  result_matrix=result_matrix[order(result_matrix[,"p_value"],decreasing = F),]
  FDR=p.adjust(result_matrix[,"p_value"],method="fdr",n=nrow(result_matrix))
  result_matrix=cbind(result_matrix,"FDR"=FDR)
  
  result_matrix=result_matrix[which(result_matrix[,"FDR"]<=FDR_threshold),]
  result_matrix[,"p_value"]=round(as.numeric(result_matrix[,"p_value"]),digits = 8)
  result_matrix[,"FDR"]=round(as.numeric(result_matrix[,"FDR"]),digits = 4)
  
  return(as.data.frame(result_matrix))
}
  
## 6. identify the up-regulated or down regulated
identify_the_genes_regulated=function(Fisher_test_result,background_genepairs){
  resistant_background=background_genepairs[["resistant"]]
  sensitive_background=background_genepairs[["sensitive"]]
  
  DEGs=as.character(Fisher_test_result$gene) 
  type_vector=c()
  for(i in 1:length(DEGs)){
    r_tem=resistant_background[[DEGs[i]]]
    s_tem=sensitive_background[[DEGs[i]]]
    
    if(length(which(s_tem==1))>length(which(r_tem==1))){type="up_regulated"}
    else if(length(which(s_tem==0))>length(which(r_tem==0))){type="down_regulated"}
    else{type="none"}
    type_vector=c(type_vector,type)
    
  }
  Fisher_test_result$regulate_type=type_vector
  return(Fisher_test_result)
}
  
