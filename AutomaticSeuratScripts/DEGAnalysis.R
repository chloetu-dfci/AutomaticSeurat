library(Seurat)
library(optparse)
# library(limma)
# Broad's R-4.1 doesn't have limma

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

Analyse_Clusters <- function(full_obj, resolution, output_dir, integration_bool){
  Idents(full_obj) <- "orig.ident"
  cells_per_sample <- CellsByCluster(sample = full_obj, integrated = integration_bool, resolution = resolution)
  write.csv(cells_per_sample, file.path(output_dir, paste0("res_", resolution, "_cells_per_sample.csv")))
  DefaultAssay(full_obj) <- "RNA"
  if(integration_bool == TRUE){
    Idents(full_obj) <- paste0("integrated_snn_res.", resolution)
  } else {
    Idents(full_obj) <- paste0("SCT_snn_res.", resolution)
  }
  marker_list <- FindAllMarkers(full_obj)
  write.csv(marker_list, file.path(output_dir, paste0("res_", resolution, "_cluster_markers.csv")))
}

### Handle input args ----
option_list <- list(
  make_option("--seuratObject",action="store"),
  make_option("--outputDir",action="store"),
  make_option("--integrate",action="store",default=FALSE),
  make_option("--resolution",action="store")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Must give a value for every single argument, except integrate
if(is.null(opt$resolution) | is.null(opt$outputDir)) { 
  stop("A mandatory argument not found\n")
}


# Reassign variables
seuratObject = opt$seuratObject
output_dir = opt$outputDir
integration_bool = opt$integrate
resolution = opt$resolution

# Testing:
# seuratObject = "seurat_obj_after_analysis.RDS"
# output_dir = "/Users/chloetu/Desktop/test_folder"
# resolution = 0.2

### Run commands ---
full_obj <- readRDS(file.path(output_dir, seuratObject))
Analyse_Clusters(full_obj, resolution, output_dir, integration_bool)

print("DEG analysis done")



