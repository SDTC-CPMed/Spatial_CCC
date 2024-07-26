devtools::install('SCEVAN-main')

library(SCEVAN)
library(data.table)

for(sample_slice in c('MEND160')){
#for(sample_slice in c('TENX40', TENX41, TENX40)){
  
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
  
  
  epi_cells = read.csv('groupFiles.txt', sep = '\t', header = FALSE)
  epi_cells = epi_cells[epi_cells$V2 %in% c('Benign', 'Benign*', 'GG1', 'GG4', 'GG2',
                                              'GG4 Cribriform', 'PIN', 'Transition_State'),]
  epi_cells <- epi_cells[grepl(slice_name, epi_cells$V1), ]
  epi_cells <- sub(slice_name, "", epi_cells$V1)
  

  
  reference_benign_cells = read.csv('Reference_selection_group_from_benign_spots.csv')
  reference_benign_cells = reference_benign_cells[reference_benign_cells$group == 'Reference',]
  reference_benign_cells <- reference_benign_cells[grepl(slice_name, reference_benign_cells$X), ]
  reference_benign_cells <- sub(slice_name, "", reference_benign_cells$X)
  
  stroma_cells = read.csv('groupFiles.txt', sep = '\t', header = FALSE)
  stroma_cells = stroma_cells[stroma_cells$V2 %in% c('Stroma'),]
  stroma_cells <- stroma_cells[grepl(slice_name, stroma_cells$V1), ]
  stroma_cells <- sub(slice_name, "", stroma_cells$V1)
  
  
  
  
  adata <- read.csv(paste('hest/input/', sample_slice, '.csv', sep = ''), header = TRUE, row.names = 'X')
  
  # In case you want to keep only epithelial cells: 
  #adata = adata[epi_cells,]  
  
  # transpose
  t_adata <- transpose(adata)
  
  # get row and colnames in order
  colnames(t_adata) <- rownames(adata)
  rownames(t_adata) <- colnames(adata)
  
  for(beta in c(1)){
    print(beta)
    try({
      results <- pipelineCNA(t_adata, par_cores = 50, beta_vega = beta, norm_cell = stroma_cells,
                             sample = paste(sample_slice, beta, sep = '_'), perc_genes = 1, SUBCLONES = TRUE,
                             SCEVANsignatures = FALSE, ClonalCN = FALSE, plotTree = FALSE)
      write.csv(results, paste('hest/output/', sample_slice, '_', beta, '_CNV_test_withReference.csv', sep = ''))
    })
  }
}

