# Function to install or update needed packages

installation_options <- function(install_copyKat = FALSE, install_devtools = FALSE, update_copyKat = FALSE){
  if (!is.logical(install_copyKat) || length(install_copyKat) != 1) {
    stop("install_copyKat is not a boolean.")
  }
  if (!is.logical(install_devtools) || length(install_devtools) != 1) {
    stop("Argument 2 is not a boolean.")
  }
  if (!is.logical(update_copyKat) || length(update_copyKat) != 1) {
    stop("update_copyKat is not a boolean.")
  }
  if(install_devtools){
    install.packages("devtools")
  }
  
  if(install_copyKat){
    library(devtools)
    install_github("navinlabcode/copykat", force = TRUE)
  }
  
  if(update_copyKat){
    remove.packages("copykat")
    detach("package:copykat")
  }
}

# Function to run copykat
run_one_copykat <- function(test_number,
                            different_folder = TRUE,
                            sample_slices_vector = c('MEND160'),
                            path = "/Documents/DATA/", 
                            copykat_distance = "euclidean", # "euclidean" # "pearson" # "spearman"
                            reference_option = "None", # "None" # "Normal" # "GG4" # "MEND156" # "allBenignRdata" # "stroma"
                            reference_path = "",
                            rdata_file_name = "",
                            adata_option = "epi", # "se" # "everything"
                            ngene_perchr = 5,
                            copykat_ks = 0.1,
                            copykat_ws = 25,
                            copykat_lowdr = 0.05,
                            copykat_cores = 50,
                            subClusters = 6){
  
  if(different_folder){
    new_folder  <- paste("copyKatTest",test_number,sep="")
    
    if (!dir.create(new_folder, showWarnings = TRUE)) {
      stop("Failed to create directory")
    }
    
    setwd(new_folder) # setting new wd
  }
  
  for(sample_slice in sample_slices_vector){
    
    print(sample_slice)
    
    
    if(sample_slice == 'MEND156'){
      slice_name = '^H2_5_'
    }
    if(sample_slice == 'MEND158'){
      slice_name = '^H2_1_'
    }
    if(sample_slice == 'MEND159'){
      slice_name = '^H1_5_'
    }
    if(sample_slice == 'MEND160'){
      slice_name = '^H1_4_'
    }
    if(sample_slice == 'MEND161'){
      slice_name = '^H1_2_'
    }
    
    if(reference_option != "None"){
      if(reference_option == "copykat"){
        benign_cells <- read.csv("/Documents/CopyKatReferenceCells.csv")
        reference_benign_cells <- benign_cells[benign_cells$copykat.pred == "diploid",]
        reference_benign_cells <- reference_benign_cells$cell.names
      }else if(reference_option == "GG4"){
        reference_file = '_ReferenceGG4.csv'
        benign_cells <- read.csv(paste(path,sample_slice,reference_file,sep=""))
        reference_benign_cells <- benign_cells[benign_cells$IsReference == "True",]
        reference_benign_cells <- sub(slice_name, "", reference_benign_cells$V1)
      }else if(reference_option == "stroma"){
        groupCells = read.csv(paste(path,'groupFiles.txt',sep=""), sep = '\t', header = FALSE)
        reference_benign_cells = groupCells[groupCells$V2 %in% c('Stroma'),]
        reference_benign_cells <- reference_benign_cells[grepl(slice_name, reference_benign_cells$V1), ]
        reference_benign_cells <- sub(slice_name, "", reference_benign_cells$V1) 
      }else if(reference_option == "MEND156"){
        reference_file <- '_Reference.csv'
        benign_cells <- read.csv(paste(path,reference_option,reference_file,sep=""))
        slice_name2 <- '^H2_5_'
        reference_benign_cells <- benign_cells[benign_cells$IsReference == "True",]
        reference_benign_cells <- sub(slice_name2, "", reference_benign_cells$V1)        
      }else if(reference_option == "allBenignRdata"){
        load(paste0(path,rdata_file_name,".RData"))
        reference_benign_cells <- se@meta.data[se@meta.data$Label %in% c("Benign"),"Barcode"]
      }else{
        reference_file = '_Reference.csv'
        benign_cells <- read.csv(paste(path,sample_slice,reference_file,sep=""))
        reference_benign_cells <- benign_cells[benign_cells$IsReference == "True",]
        reference_benign_cells <- sub(slice_name, "", reference_benign_cells$V1)
      }
    } else{
      reference_benign_cells = ""
    }
    
    if(adata_option == "epi"){
      groupCells = read.csv(paste(path,'groupFiles.txt',sep=""), sep = '\t', header = FALSE)
      epi_cells = groupCells[groupCells$V2 %in% c('Benign', 'Benign*', 'GG1', 'GG4', 'GG2',
                                                  'GG4 Cribriform', 'PIN', 'Transition_State'),]
      epi_cells <- epi_cells[grepl(slice_name, epi_cells$V1), ]
      epi_cells <- sub(slice_name, "", epi_cells$V1)
      
      
      adata <- read.csv(paste(path,sample_slice, '.csv', sep = ''), header = TRUE, row.names = 'X')
      adata <- adata[epi_cells,]
      
      # transpose
      t_adata <- transpose(adata)
      
      # get row and colnames in order
      colnames(t_adata) <- rownames(adata)
      rownames(t_adata) <- colnames(adata)
    }else if (adata_option == "se"){
      raw_counts <- se@assays$Spatial@counts
      # If you want to convert it to a regular matrix
      t_adata <- as.matrix(raw_counts)
      
    }else if(adata_option == "everything"){
      adata <- read.csv(paste(path,sample_slice, '.csv', sep = ''), header = TRUE, row.names = 'X')
    # transpose
    t_adata <- transpose(adata)
    
    # get row and colnames in order
    colnames(t_adata) <- rownames(adata)
    rownames(t_adata) <- colnames(adata)
    }
    
    for(beta in c(1)){
      print(beta)
      try({
        results <- copykat(rawmat = t_adata, 
                           norm.cell.names = reference_benign_cells,
                           distance = copykat_distance, # "euclidean" # "pearson" # "spearman"
                           ngene.chr = ngene_perchr,
                           KS.cut = copykat_ks,
                           win.size = copykat_ws,
                           LOW.DR= copykat_lowdr,
                           n.cores = 50) # Omica has 72
        hc <- results$hclustering
        
        
        CNVsData <- data.frame("class" = NaN, 
                               "subclone" = cutree(hc, k = subClusters))
        
        allnames <- rownames(CNVsData)
        for(i in 1:length(allnames)){
          allnames[i] <- paste(substr(allnames[i],1, nchar(allnames[i])-2),"-1",sep="")  
        }
        rownames(CNVsData) <- allnames
        
        preds <- results$prediction
        rownames(preds) = preds$cell.names
        
        CNVsData$class <- preds[allnames,"copykat.pred"]
        CNVsData$class <- ifelse(CNVsData$class== "diploid", "Normal", "Tumor")
        
        matrix_results <- t(results$CNAmat[,4:length(results$CNAmat)])
        rownames(matrix_results) <- allnames
        
        CNVsData$means <- rowMeans(matrix_results)
        CNVsData$mins <- apply(matrix_results, 1, min)
        CNVsData$maxs <- apply(matrix_results, 1, max)
        CNVsData$positmeans <- apply(matrix_results, 1, function(x) mean(x[x >= 0], na.rm = TRUE))
        CNVsData$negativemeans <- apply(matrix_results, 1, function(x) mean(x[x < 0], na.rm = TRUE))
        
        BigCNV <- cbind(CNVsData,matrix_results)
        
        CNVbyClassGroup2 <- BigCNV %>%
          group_by(subclone) %>%
          summarise(across(7:(ncol(BigCNV)-1), mean))
        
        distances_centroids <- rowMeans(CNVbyClassGroup2^2)
        orderClusters <- order(distances_centroids)
        df_order <- data.frame(orderClusters,distances_centroids[orderClusters])
        
        un <- unique(CNVsData[CNVsData$class=="Normal","subclone"])
        CNVsData$isBest = FALSE
        for(val in orderClusters){
          if(val %in% un){
            CNVsData[CNVsData$subclone == val,"isBest"] = TRUE
            break
          }
        }
        
        write.csv(CNVsData, paste(sample_slice, '_', beta, '_CNV_CopyKatWithRef',test_number,'.csv', sep = ''))
        write.csv(df_order, paste(sample_slice, '_', beta, '_CNV_OrderClusters',test_number,'.csv', sep = ''))
      })
    }
  }
  if(different_folder){
    setwd("..") # going back to parent wd
  }
}

# Installing packages
installation_options(TRUE,TRUE,FALSE)

# Loading packages
library(stats)
library(copykat)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# Running copykat
run_one_copykat(test_number = 85,
                different_folder = FALSE,
                sample_slices_vector = c('MEND160'),
                ngene_perchr = 5,
                copykat_distance = "pearson",
                reference_option = "allBenignRdata",
                rdata_file_name = "seurat_st_H1_4_annotated2",
                adata_option = "se",
                copykat_lowdr = 0.08,
                subClusters = 6)

