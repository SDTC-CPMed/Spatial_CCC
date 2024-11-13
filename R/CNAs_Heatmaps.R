library(parallelDist)
library(viridis)
library(copykat)

test_numbers = c('INT25' = "22",
                 'INT26' = "23",
                 'INT27' = "21",
                 'INT28' = "24",
                 'TENX40' = "25",
                 'TENX46' = "26",
                 'MEND151' = "19",
                 'MEND156' = "17",
                 'MEND158' = "18",
                 'MEND160' = "16",
                 'MEND161' = "20")

data_path = 'C:/Users/dandia/OneDrive - Karolinska Institutet/Documents/Github/Spatial_CCC/DATA/'
samples = c('INT25','INT26','INT27','INT28','TENX40','TENX46','MEND151','MEND156','MEND158','MEND160','MEND161')
for(sample in samples){
  test_number = test_numbers[sample]

  input_path = paste0('C:/Users/dandia/OneDrive - Karolinska Institutet/Documents/Github/Spatial_CCC/DATA/CopyKatDataPlots/',sample, '/',test_number,'/')
  CNVMatrix = read.csv(paste0(input_path,sample,"_1_CNVsMatrix",test_number,".csv"))
  rownames(CNVMatrix) = CNVMatrix$X
  CNVMatrix$X = NULL
  
  PT_Data = read.csv(paste0(data_path,sample,"_HA_PT.csv"))
  rownames(PT_Data) = PT_Data$X
  PT_Data$X = NULL
  PT_Data$Label_orig2 = NULL
  
  chr <- unname(as.matrix(CNVMatrix)[1, ]) %% 2+1
  rbPal1 <- colorRampPalette(c('black','grey'))
  CHR <- rbPal1(2)[as.numeric(chr)]
  chr1 <- cbind(CHR,CHR)
  
  CNVMatrix = CNVMatrix[-c(1,2,3),]
  CNVMatrix$Index <- rownames(CNVMatrix)
  PT_Data$Index <- rownames(PT_Data)
  
  
  # Perform a left join, keeping all rows from CNVMatrix
  result <- merge(CNVMatrix, PT_Data, by = "Index", all.x = TRUE)
  rownames(result) = result$Index
  result$Index = NULL
  result$Label_orig[is.na(result$Label_orig)] = ""
  
  col_order <- names(result)
  col_order[match(c("Label_orig", "dpt_pseudotime"), col_order)] <- c("dpt_pseudotime", "Label_orig")
  result <- result[, col_order]
  
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 999)
  
  col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))
  labels_orig = result[["Label_orig"]]
  labels_orig[labels_orig == ""] = NA
  label_colors_5 <- c('Benign' = "blue", 
                      'GG1' = "peachpuff", 
                      'GG2' = "#FF6666", 
                      'GG3' = "red", 
                      'GG4 Cribriform' = "purple")
  
  colors_5 <- label_colors_5[labels_orig]
  
  viridis_colors <- viridis(100)  # 100 colors from the viridis palette
  num_colors <- viridis_colors[cut(result[["dpt_pseudotime"]], breaks = seq(0, 1, length.out = 101), labels = FALSE)]
  
  colors <- rbind(colors_5,num_colors) #colors_2,
  rownames(colors) = c("Histological Annotations", "Pseudotime")
  
  jpeg(paste0(data_path,"heatmap_",sample,"_copykat.jpeg"), height=2500, width=4000, res=100)
  heatmap.3(as.matrix(result[,1:(dim(result)[2] - 2)]),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = "pearson")), hclustfun = function(x) hclust(x, method="ward.D"),
            ColSideColors=chr1,RowSideColors=colors,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, main=paste0("Copy Number Profiles (",sample,")"), cex.main=4, margins=c(15,15))
  dev.off()
}
