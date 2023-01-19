library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(grid)
library(stringr)
library(patchwork)
library(Polychrome)
library(Matrix)
library(optparse)
# library(limma)
# Broad's R-4.1 doesn't have limma

# Default values
cluster_res <- c(0.2, 0.4, 0.6)
# 1:30 is often good enough.
# dimensions <- c(1:30)

##### Data Loading Functions ----
# GetNewCellNames function takes in a vector of chracters and a number.
# Replaces the barcode postfix (-1) with the sample number (-<sample-number>)
GetNewCellNames <- function(old_cell_names, sample_no){
  for(i in c(1:length(old_cell_names))){
    old_cell_names[i] <- gsub("-.*", 
                              paste0("-", sample_no), 
                              old_cell_names[i])
  }
  return (old_cell_names)
}


LoadData <- function(sample_names, file_names, directory){
  sample_list <- list()
  # Read in 10x data
  sample_list <- lapply(file.path(directory, paste0(file_names, "/filtered_feature_bc_matrix")), FUN = Read10X)
  # Rename entries in list
  names(sample_list) <- sample_names
  # Create seurat object
  for(i in c(1:length(sample_names))){
    sample_list[[i]] <- CreateSeuratObject(sample_list[[i]], sample_names[i])}
  # Rename cells in list
  sample_list <- lapply(sample_list, FUN = function(X){ 
    old_cell_names <- rownames(X@meta.data)
    sample_num <- str_extract(levels(X$orig.ident), "\\d.*")
    new_cell_names <- GetNewCellNames(old_cell_names, sample_num)
    X <- RenameCells(X, new.names = new_cell_names)
  })
  return(sample_list)
}
##### Data Loading Functions ####

##### QC Functions ----

countCell_vector <- function(seurat_obj, colName){
  # extract all unique samples in the seurat object
  samples <- Idents(seurat_obj)
  samples <- unique(samples)
  # create dataframe with unique samples as row names
  countDf <- data.frame(matrix(ncol = 2, nrow = 0))
  for (sample in samples){
    count <- c(sample, dim(subset(seurat_obj@meta.data, Idents(seurat_obj) == sample))[1])
    countDf <- rbind(countDf, count)
  }
  # rename
  colnames(countDf) <- c("sample", paste0(colName, "Freq"))
  countDf[,2] <- as.numeric(countDf[,2])
  return (countDf)
}

countCell_plot <- function(cellCountDF, x_axis, y_axis, title){
  plot <- ggplot(cellCountDF, aes(x = cellCountDF[,x_axis], y = cellCountDF[,y_axis])) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    geom_label(aes(label = cellCountDF[,y_axis]), y = 100, size=3) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_hline(yintercept = mean(cellCountDF[,y_axis]), color = "red")+
    ylab("Freq")+
    xlab("Sample")+
    ggtitle(title)
  return(plot)
}

feature_facetPlot <- function(sample.metadata, feature, cutoff.df){
  plot <- ggplot(sample.metadata, aes(x=sample.metadata[,feature], fill = orig.ident)) + 
    geom_density(alpha = 0.2, show.legend = FALSE) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    xlab(paste(feature,"(log scaled)")) +
    facet_wrap(vars(orig.ident),scales="free_x") +
    geom_vline(data = cutoff.df,
               aes(xintercept = X3SD_above), 
               colour = "limegreen") +
    geom_vline(data = cutoff.df,
               aes(xintercept = X2SD_below), 
               colour = "limegreen")
  return(plot)
}

CreateSDTable <- function(seurat_obj, metrics_list){
  sdTable <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(sdTable) <- c("mean", "SD", "2SD_below",
                         #"1SD_below", "1SD_above","2SD_above",
                         "3SD_above"
                         )        
  for(metric in metrics_list){
    meanVal <- mean(as.numeric(seurat_obj@meta.data[,metric]))
    sdVal <- sd(as.numeric(seurat_obj@meta.data[,metric]))
    sdTable[metric,] <- c(meanVal, sdVal,meanVal - 2*sdVal, 
                          #meanVal - 1*sdVal, meanVal + 1*sdVal, meanVal + 2*sdVal,
                          meanVal + 3*sdVal
                          )
  }
  return(sdTable)
}

# Extract SD values into one vector
ExtractSDs_intoDF <- function(SD_list, specific_metric, specific_SD){
  specificSD_allSamples <- c()
  for(sample_name in names(SD_list)){
    specificSD_allSamples <- c(specificSD_allSamples, SD_list[[sample_name]][specific_metric, specific_SD])
  }
  return(specificSD_allSamples)
}

# Filter for cells with a list of requirements
ApplyFilter <- function(sample_list, cutoff_list){
  # Create an empty list with same length as list of tumor samples
  sample_list.filtered <- vector("list", length = length(sample_list))
  # Update entries in new list
  for(i in c(1:length(sample_list))){
    sample_list.filtered[[i]] <- subset(sample_list[[i]], cells = cutoff_list[[i]])
    
  }
  names(sample_list.filtered) <- names(sample_list)
  return(sample_list.filtered)
}

## Select cells with > or < X nFeature/nCount
Select_SD <- function(sample_list, SD_list, rowName = "nCount_RNA", colName = "3SD_above", limit_dir = "lower"){
  # Create an empty list with same length as list of tumor samples
  pass_SD_cutoff <- vector("list", length = length(sample_list))
  # Update entries in new list
  for(i in c(1:length(sample_list))){
    if(rowName == "nCount_RNA"){
      if(limit_dir == "upper"){
        pass_SD_cutoff[[i]] <- WhichCells(sample_list[[i]], expression = nCount_RNA < SD_list[[i]][rowName,colName])}
      else if(limit_dir == "lower"){
        pass_SD_cutoff[[i]] <- WhichCells(sample_list[[i]], expression = nCount_RNA > SD_list[[i]][rowName,colName])}
    }
    else if(rowName == "nFeature_RNA"){
      if(limit_dir == "upper"){
        pass_SD_cutoff[[i]] <- WhichCells(sample_list[[i]], expression = nFeature_RNA < SD_list[[i]][rowName,colName])}
      else if(limit_dir == "lower"){
        pass_SD_cutoff[[i]] <- WhichCells(sample_list[[i]], expression = nFeature_RNA > SD_list[[i]][rowName,colName])}
    }
  }
  names(pass_SD_cutoff) <- names(sample_list)
  return(pass_SD_cutoff)
}


Update_lower_SD <- function(SD_table){
  count_features <- SD_table%>%
    filter(rownames(.)%in%c("nFeature_RNA","nCount_RNA"))%>%
    mutate(`2SD_below` = case_when(
      `2SD_below` > 0 ~ 100,
      is.na(`2SD_below`) ~ 100,
      TRUE ~ as.numeric(`2SD_below`)))
  mito <- SD_table%>%
    filter(rownames(.)==("percent.mt"))%>%
    mutate(`2SD_below` = 0)
  SD_table <- rbind(count_features, mito)
  return(SD_table)
}

##### QC Functions #####

##### Clustering Functions ----
PlotCluster <- function(sample, integrated = FALSE, res = "ident", title = NA){
  if(res == "ident"){
    num_cols <- length(unique(Idents(sample)))
    col_palette = as.vector(createPalette(num_cols,  c("#ff0000", "#00ff00", "#0000ff")))
    
    clusterPlot <- DimPlot(sample, cols = col_palette, 
                           group.by = "ident", 
                           label.size = 3, 
                           label = TRUE, repel = TRUE) +
      ggtitle(if(is.na(title)){paste0(deparse(substitute(sample)))} else {title})
  } else{
    res_num <- paste0(
      if(integrated == TRUE)
      {"integrated_snn_res."} 
      else {"SCT_snn_res."
      }, res)
    num_cols <- length(unique(sample[[res_num]][,1]))
    col_palette = as.vector(createPalette(num_cols,  c("#ff0000", "#00ff00", "#0000ff")))
    
    clusterPlot <- DimPlot(sample, cols = col_palette, 
                           group.by = res_num, 
                           label.size = 3,
                           label = TRUE, repel = TRUE) + 
      ggtitle(if(is.na(title)){paste0("res: ", res, ", ", deparse(substitute(sample)))} else {title})
  }
  print(clusterPlot)
}

CellsByCluster <- function(sample, integrated = FALSE, resolution){
  if (integrated == TRUE){
    cells_per_sample <- as.data.frame.matrix(t(table(sample@meta.data[,"orig.ident"], 
                                                     sample@meta.data[,paste0("integrated_snn_res.", resolution)])))
  } else {
    cells_per_sample <- as.data.frame.matrix(t(table(sample@meta.data[,"orig.ident"], 
                                                     sample@meta.data[,paste0("SCT_snn_res.", resolution)])))
  }
  return (cells_per_sample)
}

PlotCellsByCluster <- function(cells_per_sample, resolution){
  cells_per_sample <- rownames_to_column(cells_per_sample, "cluster")
  cells_per_sample <- pivot_longer(cells_per_sample, cols = -cluster)
  colnames(cells_per_sample) <- c("cluster", "sample", "count")
  
  clusterPlot <- ggplot(cells_per_sample, aes(x = cluster, y = count, fill = sample)) + 
    geom_bar(stat = "identity", position="fill") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme_classic() +
    ggtitle(paste0("Proportion of cells in each sample and cluster at res ", resolution))
  
  print(clusterPlot)
}

##### Clustering Functions #####

####### Big Load Data and QC Function ----
Loading_and_QC <- function(output_dir,sample_names,file_names,data_directory, species){
  sample_list <- LoadData(sample_names, file_names, data_directory)
  
  if (species == "mouse"){
    mt_pattern = "^mt-"
    hb_pattern = "^Hb(a|b)-"
  } else if (species == "human"){
    mt_pattern = "^MT-"
    hb_pattern = "^HB(a|b)-"
  } else{
    stop("Unknown species. Specify mouse or human only.\n")
  }

  for (i in c(1:length(sample_list))){
    # check for % of counts in our cells that map to the mitochondrial genome
    sample_list[[i]][["percent.mt"]] <- PercentageFeatureSet(sample_list[[i]], pattern = mt_pattern)
    # check for % of counts in our cells that have red blood cell genes
    sample_list[[i]][["percent.hb"]] <- PercentageFeatureSet(sample_list[[i]], pattern = hb_pattern)
  }
  # concatenates gene exp matricies
  sample_merged <- merge(x = sample_list[[1]],
                      y = sample_list[c(2:length(sample_list))])

  ###### Start QC
  ## count number of cells in each replicate
  cellCount <- countCell_vector(sample_merged, "")
  
  # Plot number of cells in each replicate
  print(countCell_plot(cellCount,"sample","Freq","Number of cells in samples before filtering"))
  
  ## Violin plots to visualize differences in count/feature/percent mt
  count_vln <- VlnPlot(sample_merged, features = "nCount_RNA", 
                       pt.size = FALSE, log = TRUE, ) + 
    NoLegend() + 
    theme(axis.text = element_text(size = 5), 
          plot.title = element_text(size = 7), 
          axis.title.x = element_text(size = 5))

  feature_vln <- VlnPlot(sample_merged, features = "nFeature_RNA", 
                         pt.size = FALSE, log = TRUE) + 
    NoLegend() + 
    theme(axis.text = element_text(size = 5), 
          plot.title = element_text(size = 7), 
          axis.title.x = element_text(size = 5))
  
  mt_vln <- VlnPlot(sample_merged, features = "percent.mt", 
                    pt.size = FALSE) + 
    NoLegend() + 
    theme(axis.text = element_text(size = 5), 
          plot.title = element_text(size = 7), 
          axis.title.x = element_text(size = 5))
  
  ## printing violin plots in two pages
  print(count_vln + feature_vln)
  
  print(mt_vln)

  ### Creating density plots for nCount/Feature/MT
  testMetrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
  
  SD_list <- lapply(sample_list, CreateSDTable, testMetrics)
  
  count_vline <- data.frame(orig.ident = unique(sample_merged@meta.data$orig.ident),  # Create data for lines
                            X2SD_below = ExtractSDs_intoDF(SD_list, "nCount_RNA", "2SD_below"), 
                            X3SD_above = ExtractSDs_intoDF(SD_list, "nCount_RNA", "3SD_above"))
  # change all 0's into NA, since can't log negative numbers
  count_vline[(count_vline) < 0] <- NA
  
  print(feature_facetPlot(sample_merged@meta.data,"nCount_RNA",count_vline))
  
  
  feature_vline <- data.frame(orig.ident = unique(sample_merged@meta.data$orig.ident),  # Create data for lines
                              X2SD_below = ExtractSDs_intoDF(SD_list, "nFeature_RNA", "2SD_below"),
                              X3SD_above = ExtractSDs_intoDF(SD_list, "nFeature_RNA", "3SD_above"))
  # change all 0's into NA, since can't log negative numbers
  feature_vline[(feature_vline) < 0] <- NA
  
  print(feature_facetPlot(sample_merged@meta.data,"nFeature_RNA",feature_vline))
  
  ##### Filtering
  # Select cells with < 25 %mit
  pass_mt_cutoff <- lapply(sample_list, WhichCells, expression = percent.mt < 25)
  pass_hb_cutoff <- lapply(sample_list, WhichCells, expression = percent.hb < 10)
  
  # Filter cells with < 25 %mit and < 10% hb
  sample_list.mt_filtered <- ApplyFilter(sample_list, pass_mt_cutoff)
  sample_list.mt_filtered <- ApplyFilter(sample_list.mt_filtered, pass_hb_cutoff)
  
  # concatenates gene exp matricies
  sample_merged.mt_filtered <- merge(x = sample_list.mt_filtered[[1]],
                                     y = sample_list.mt_filtered[c(2:length(sample_list.mt_filtered))])
  
  # Filter for cells that pass lower SD cutoffs
  # If the lower SD (for both count and feature) are negative, set lower SD to 100
  SD_list <- lapply(SD_list, Update_lower_SD)
  
  pass_LowernCount_cutoff <- Select_SD(sample_list.mt_filtered, SD_list, "nCount_RNA", "2SD_below", "lower")
  pass_LowernFeature_cutoff <- Select_SD(sample_list.mt_filtered, SD_list, "nFeature_RNA", "2SD_below", "lower")
  
  sample_list.lowerSD_filtered <- ApplyFilter(sample_list.mt_filtered, pass_LowernCount_cutoff)
  sample_list.lowerSD_filtered <- ApplyFilter(sample_list.lowerSD_filtered, pass_LowernFeature_cutoff)
  
  # Filter for cells that pass upper SD cutoffs
  pass_UppernCount_cutoff <- Select_SD(sample_list.lowerSD_filtered, SD_list, "nCount_RNA", "3SD_above", "upper")
  pass_UppernFeature_cutoff <- Select_SD(sample_list.lowerSD_filtered, SD_list, "nFeature_RNA", "3SD_above", "upper")
  
  sample_list.fully_filtered <- ApplyFilter(sample_list.lowerSD_filtered, pass_UppernCount_cutoff)
  sample_list.fully_filtered <- ApplyFilter(sample_list.fully_filtered, pass_UppernFeature_cutoff)
  
  # concatenates gene exp matricies
  sample_merged.fully_filtered <- merge(x = sample_list.fully_filtered[[1]],
                                        y = sample_list.fully_filtered[c(2:length(sample_list.fully_filtered))])
  
  # count number of cells in each replicate left
  cellCount.fully_filtered <- countCell_vector(sample_merged.fully_filtered, "")
  
  # Plot number of cells left in each replicate
  print(countCell_plot(cellCount.fully_filtered,"sample","Freq","Number of cells in samples after filtering"))
  
  # # Plot nCounts
  # feature_facetPlot(sample_merged.fully_filtered@meta.data,"nCount_RNA",count_vline)
  # 
  # # Plot nFeatures
  # feature_facetPlot(sample_merged.fully_filtered@meta.data,"nFeature_RNA",feature_vline)
  
  return(sample_list.fully_filtered)
}
####### Load Data and QC ######

####### Big SCTransform and RPCA Integration Function ----
Normalize_and_Integrate <- function(sample_list.fully_filtered){
  sample_list.fully_filtered <- lapply(X = sample_list.fully_filtered, FUN = SCTransform)
  features <- SelectIntegrationFeatures(object.list = sample_list.fully_filtered, nfeatures = 3000)
  sample_list.fully_filtered <- PrepSCTIntegration(object.list = sample_list.fully_filtered, 
                                                   anchor.features = features)
  sample_list.fully_filtered <- lapply(X = sample_list.fully_filtered, FUN = RunPCA, features = features)
  anchors <- FindIntegrationAnchors(object.list = sample_list.fully_filtered, normalization.method = "SCT",
                                    anchor.features = features, reduction = "rpca")
  integrated_obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  integrated_obj <- RunPCA(integrated_obj, verbose = FALSE)
  integrated_obj <- RunUMAP(integrated_obj, reduction = "pca", dims = 1:30)
  integrated_obj <- FindNeighbors(integrated_obj, dims = 1:30)
  
  print(DimPlot(integrated_obj, reduction = "umap", order = TRUE, group.by = "orig.ident") + ggtitle("Integrated samples")) 
  print(DimPlot(object = integrated_obj, reduction = "umap", split.by = "orig.ident", 
          order = TRUE, group.by = "orig.ident") + 
    NoLegend() + 
    facet_wrap(vars(orig.ident))
    + ggtitle("Integrated samples"))
  
  return(integrated_obj)
}
####### Big SCTransform and RPCA Integration Function #####

####### Big SCTransform and Merging Function ----
Normalize_and_Merge <- function(sample_list.fully_filtered){
  sample_list.fully_filtered <- lapply(X = sample_list.fully_filtered, FUN = SCTransform)
  merged_obj <- merge(sample_list.fully_filtered[[1]], y = sample_list.fully_filtered[c(2:length(sample_list.fully_filtered))])
  VariableFeatures(merged_obj[["SCT"]]) <- rownames(merged_obj[["SCT"]]@scale.data)
  merged_obj <- RunPCA(merged_obj, verbose = FALSE)
  merged_obj <- RunUMAP(merged_obj, reduction = "pca", dims = 1:30)
  merged_obj <- FindNeighbors(merged_obj, dims = 1:30)
  
  print(DimPlot(merged_obj, reduction = "umap", order = TRUE, group.by = "orig.ident") + 
          ggtitle("Merged samples"))
  print(DimPlot(object = merged_obj, reduction = "umap", split.by = "orig.ident", 
                order = TRUE, group.by = "orig.ident") + 
          NoLegend() + 
          facet_wrap(vars(orig.ident)) +
          ggtitle("Merged samples"))
  
  return(merged_obj)
}
####### Big SCTransform and Merging Function #####

####### Big Cluster and DEG extraction Function ----
Cluster_and_Annotate <- function(full_obj, gene_vec, output_dir, integration_bool){
  for(gene in gene_vec){
    DefaultAssay(full_obj) <- "RNA"
    print(FeaturePlot(full_obj, gene))
  }
  for(i in cluster_res){
    if(integration_bool == TRUE){DefaultAssay(full_obj) <- "integrated"} else{DefaultAssay(full_obj) <- "SCT"}
    full_obj<-FindClusters(object = full_obj, verbose = FALSE, resolution = i)
    PlotCluster(full_obj, res = i, integrated = integration_bool)
    if(integration_bool == TRUE){
      Idents(full_obj) <- paste0("integrated_snn_res.", i)
    } else {
      Idents(full_obj) <- paste0("SCT_snn_res.", i)
      # PrepSCTFindMarkers function is not available in Seurat 4.0.3, which is what the Broad's R-4.1 DotKit has loaded.
      # full_obj <- PrepSCTFindMarkers(full_obj)
    }
    cells_per_sample <- CellsByCluster(sample = full_obj, integrated = integration_bool, resolution = i)
    # write.csv(cells_per_sample, file.path(output_dir, paste0("res_", i, "_cells_per_sample.csv")),
    #           row.names = FALSE)
    PlotCellsByCluster(cells_per_sample = cells_per_sample, resolution = i)
  }
  return(full_obj)
}
####### Big Cluster and DEG extraction Function ####

# Find_All_Markers <- function(full_obj, ){
#   DefaultAssay(full_obj) <- "RNA"
#   marker_list <- FindAllMarkers(full_obj)
#   write.csv(marker_list, file.path(output_dir, paste0("res_", i, "_cluster_markers.csv")),
#             row.names = FALSE)
# }

### Handle input args ----
option_list <- list(
  make_option("--sampleNames",action="store"),
  make_option("--fileNames",action="store"),
  make_option("--outputDir",action="store"),
  make_option("--dataDir",action="store"),
  make_option("--pdfName",action="store"),
  make_option("--geneVec",action="store"),
  make_option("--integrate",action="store",default="false"),
  make_option("--species", action= "store")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Must give a value for every single argument, except integrate
if(is.null(opt$sampleNames) | is.null(opt$fileNames) | is.null(opt$geneVec) | is.null(opt$outputDir) |
   is.null(opt$dataDir) | is.null(opt$pdfName) | is.null(opt$species)) { 
  stop("A mandatory argument not found\n")
}

# Reassign variables
sampleNames = opt$sampleNames
fileNames = opt$fileNames
geneVec = opt$geneVec
output_dir = opt$outputDir
data_dir = opt$dataDir
pdf_name = opt$pdfName
integration_bool = as.logical(opt$integrate)
species = opt$species

# Parse strings into character vectors
sample_names <- str_split(sampleNames, ",")[[1]]
file_names <- str_split(fileNames, ",")[[1]]
gene_vec <- str_split(geneVec, ",")[[1]]

# Number of sample names and files given must be the same
if (length(sample_names) != length(file_names)){ 
  stop("Mismatch between the number of sample names and files given\n")
}

# Print out
print(paste0('Sample names: ', toString(sample_names)))
print(paste0('File names: ', toString(file_names)))
print(paste0('Output directory: ', output_dir))
print(paste0('Data directory: ', data_dir))
print(paste0('Name of PDF: ', pdf_name))
print(paste0('Genes to plot: ', toString(gene_vec)))
print(paste0('Integrate? ', integration_bool))
print(paste0("Species: ", species))

### Run commands ---
pdf(paste0(output_dir,pdf_name), width = 6, height = 6)

# Testing
# sample_names = c("spleen1", "spleen2", "spleen3")
# file_names = c("Pool115-1_2-Spleen1CNT", "Pool115-1_3-Spleen2CNT", "Pool115-1_4-Spleen3CNT")
# data_dir = "/Users/chloetu/Desktop/Anna_exp/cDNA/"
# output_dir = "/Users/chloetu/Desktop/test_folder"
# species = "mouse"
# gene_vec = c("Atg7", "Cd8a")
# pdf_name = "test.pdf"
# integration_bool = TRUE

filtered_samples.list <- Loading_and_QC(output_dir,sample_names,file_names, data_dir, species)
print("Loading and QC done")

if(integration_bool == TRUE){
  full_obj <- Normalize_and_Integrate(filtered_samples.list)
  print("Normalization and integration done")
} else{
  full_obj <- Normalize_and_Merge(filtered_samples.list)
  print("Normalization and merging done")
}

rm(filtered_samples.list)
gc()

full_obj <- Cluster_and_Annotate(full_obj, gene_vec, output_dir, integration_bool)
print("Clustering done")

dev.off()

saveRDS(full_obj, file.path(output_dir, "seurat_obj_after_analysis.RDS"))
print("Object saved")




