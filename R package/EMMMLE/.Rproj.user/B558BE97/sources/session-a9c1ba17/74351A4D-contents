##cluster_init
cluster_init<-function(expression_profile,
                        celltypes_num,
                        celltypes_ref = NULL,
                        ls_est = "TSS",
                        gene_GRN,
                        HVG_model_num,
                        zero_GRN = NULL,
                        preprocess_Control = list(HVG_num = 2000,npc = 50,run_umap = FALSE,label_umap = NULL, cluster_method = "SNN"),
                        core_num = 1){
  #Note: user must provide either "celltypes_num" or "celltypes_ref" in "Cluster_Control".
  #Note: user must gene_GRN

  ##Step0: Hyper-parameter definition
  ##-----------------------------------------------------------------------------

  ##preprocess_Control
  #######################################
  HVG_num<-ifelse(is.null(preprocess_Control$HVG_num),2000,preprocess_Control$HVG_num)
  HVG_num <- min(nrow(expression_profile),HVG_num)
  npc<-ifelse(is.null(preprocess_Control$npc),50,preprocess_Control$npc)
  npc<-min(npc,nrow(expression_profile))
  run_umap<-ifelse(is.null(preprocess_Control$run_umap),FALSE,preprocess_Control$run_umap)
  if(is.null(preprocess_Control$label_umap)){
    label_umap<-NULL
  }else{
    label_umap<-preprocess_Control$label_umap
  }
  verbose<-FALSE
  cluster_method<-ifelse(is.null(preprocess_Control$cluster_method),"SNN",preprocess_Control$cluster_method)

  if(cluster_method == "SNN"){
    resolution<-0.8
  }else{
    resolution<-NULL
  }

  ##
  if(is.null(zero_GRN)){
    zero_GRN<-matrix(0,nrow = 1,ncol = 2)
    zero_GRN_use = FALSE
  }else{
    zero_GRN_use = TRUE
  }

  #######################################

  ##-----------------------------------------------------------------------------


  #Step1 : Seurat standard workflow
  ##-----------------------------------------------------------------------------
  if(is.null(rownames(expression_profile))){
    rownames(expression_profile)<-paste("Gene",1:nrow(expression_profile),sep = "")
  }
  if(is.null(colnames(expression_profile))){
    colnames(expression_profile)<-paste("Cell",1:ncol(expression_profile),sep = "")
  }

  EmPLNet_intseurat<-CreateSeuratObject(counts = expression_profile,
                                          project = "EmPLNet_intseurat")
  EmPLNet_intseurat<-NormalizeData(EmPLNet_intseurat, verbose = verbose)

  EmPLNet_intseurat<-FindVariableFeatures(EmPLNet_intseurat, selection.method = "vst",
                                            nfeatures = max(HVG_num,HVG_model_num), verbose = verbose)

  HVG_list<-VariableFeatures(EmPLNet_intseurat)
  if(ls_est == "TSS"){
    S_depth<-as.vector(Matrix::colSums(expression_profile))/1000
  }
  if(ls_est == "GMPR"){
    S_depth<-as.vector(GMPR(as.matrix(Matrix::t(EmPLNet_intseurat@assays$RNA@counts[which(rownames(EmPLNet_intseurat@assays$RNA@counts) %in% HVG_list[1:HVG_num]),]))))
  }
  if(ls_est == "none"){
    S_depth<-rep(1,length(S_depth))
  }

  EmPLNet_intseurat <- ScaleData(EmPLNet_intseurat, verbose = verbose)
  suppressWarnings(EmPLNet_intseurat<- RunPCA(EmPLNet_intseurat, npcs = npc, verbose = verbose))
  pca_embedding<-NULL
  if(cluster_method == "Kmeans"){
    pca_embedding<-matrix(as.vector(EmPLNet_intseurat@reductions$pca@cell.embeddings),nrow = dim(EmPLNet_intseurat@reductions$pca@cell.embeddings)[1],ncol = dim(EmPLNet_intseurat@reductions$pca@cell.embeddings)[2],byrow = FALSE)
  }
  ##
  p_umap<-NULL
  umap_embedding<-NULL
  if(run_umap){
    suppressWarnings(EmPLNet_intseurat<- RunUMAP(EmPLNet_intseurat, reduction = "pca", dims = 1:min(npc,dim(EmPLNet_intseurat@reductions$pca@cell.embeddings)[2]),verbose = verbose))

    if(!is.null(label_umap)){
      EmPLNet_intseurat$seurat_clusters<-label_umap
      p_umap <- DimPlot(EmPLNet_intseurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
                        repel = TRUE) + NoLegend()
    }
    umap_embedding<-matrix(as.vector(EmPLNet_intseurat@reductions$umap@cell.embeddings),nrow = dim(EmPLNet_intseurat@reductions$umap@cell.embeddings)[1],ncol = dim(EmPLNet_intseurat@reductions$umap@cell.embeddings)[2],byrow = FALSE)
  }

  ##-----------------------------------------------------------------------------


  #Step2: Initial the U and pi
  ##-----------------------------------------------------------------------------
  resolution_0<-NULL
  if(is.null(celltypes_ref)){
    if(cluster_method == "SNN"){
      #1. use seurat SNN graph algorithm to initial celltypes_label
      EmPLNet_intseurat<-FindNeighbors(EmPLNet_intseurat, dims = 1:npc,verbose = verbose)
      EmPLNet_intseurat<-FindClusters(EmPLNet_intseurat, resolution = resolution,verbose = verbose)
      resolution_0<-resolution
      if(!is.null(celltypes_num)){
        while (length(unique(EmPLNet_intseurat$seurat_clusters)) != celltypes_num) {
          if(length(unique(EmPLNet_intseurat$seurat_clusters)) > celltypes_num){
            resolution_0<-resolution_0 * 0.9
            EmPLNet_intseurat<-FindClusters(EmPLNet_intseurat, resolution = resolution_0,verbose = verbose)
          }else{
            resolution_0<-resolution_0 * 1.1
            EmPLNet_intseurat<-FindClusters(EmPLNet_intseurat, resolution = resolution_0,verbose = verbose)
          }
        }
      }
      celltypes_label<-as.numeric(EmPLNet_intseurat$seurat_clusters)
      pi_vec<-as.vector(table(celltypes_label))/length(celltypes_label)

    }
    if(cluster_method == "Kmeans"){
      suppressWarnings(celltypes_label<-kmeans(pca_embedding,centers = celltypes_num,nstart = celltypes_num)$cluster)
      EmPLNet_intseurat$seurat_clusters<-celltypes_label
      pi_vec<-as.vector(table(celltypes_label))/length(celltypes_label)
    }
  }else{
    celltypes_label<-as.numeric(as.factor(celltypes_ref))
    pi_vec<-as.vector(table(celltypes_label))/length(celltypes_label)

  }
  celltypes_num<-length(unique(celltypes_label))
  ##-----------------------------------------------------------------------------


  ##Step3: Determine the feature for GRN of interest
  ##-----------------------------------------------------------------------------
  HVG_model<-HVG_list[0:HVG_model_num]
  if(length(HVG_model) == 0){
    gene_model<-gene_GRN
  }else{
    gene_model<-c(gene_GRN,setdiff(HVG_model,gene_GRN))
  }
  gene_model_index<-c()
  for(j in 1:length(gene_model)){
    gene_model_index<-c(gene_model_index,which(rownames(EmPLNet_intseurat@assays$RNA@counts) == gene_model[j]))
  }
  ##-----------------------------------------------------------------------------


  ##Step4: Initial the mean, M and S of each group
  ##-----------------------------------------------------------------------------
  expression_count_choose<-as(Matrix::t(EmPLNet_intseurat@assays$RNA@counts[gene_model_index,]),"sparseMatrix")
  rownames(expression_count_choose)<-colnames(EmPLNet_intseurat@assays$RNA@counts)
  rm(EmPLNet_intseurat);gc()
  colnames(expression_count_choose)<-gene_model
  expression_normal_choose<-as.matrix(log((expression_count_choose + 1) / as.vector(S_depth)))

  mu_mat<-matrix(NA,nrow = celltypes_num,ncol = ncol(expression_normal_choose))
  for(g in 1:celltypes_num){
    mu_mat[g,]<- colMeans(expression_normal_choose[which(celltypes_label == g),])
  }
  ##
  sigma_mat_list<-list()

  for(g in 1:celltypes_num){
    observed_data_group<-expression_normal_choose[which(celltypes_label == g),]
    num_cell_group<-dim(observed_data_group)[1]
    error_mat_group<-observed_data_group-matrix(rep(mu_mat[g,],num_cell_group),num_cell_group,dim(observed_data_group)[2],byrow = TRUE)
    sigma_mat_list[[g]]<-(t(error_mat_group) %*% error_mat_group)/num_cell_group
  }
  ##
  rm(expression_normal_choose);gc()

  gene_GRN_index_use<-c(1:length(gene_GRN))
  EmPLNet_list = list(celltypes_label = celltypes_label,
                        mu_mat = mu_mat,
                        sigma_mat_list=sigma_mat_list,
                        pi_vec = pi_vec,
                        obs_mat = expression_count_choose,
                        S_depth = S_depth,
                        HVG_list = HVG_list,
                        gene_GRN = gene_GRN, gene_GRN_index_use = gene_GRN_index_use,
                        p_umap = p_umap,
                        zero_GRN = zero_GRN, zero_GRN_use = zero_GRN_use)
  return(EmPLNet_list)
}
