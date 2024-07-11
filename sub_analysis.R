suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(magrittr)))

## ===== define a function of sub analysis ===== ##
sub_analysis <- function(seurat_obj, population, output_dir){
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # subset the seurat object and prep for normalization
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_sub <- subset(seurat_obj, subset = annotation %in% population)
  seurat_sub_split <- SplitObject(seurat_sub, split.by = "orig.ident")
  
  # normalization and merge
  seurat_sub_split <- lapply(X = seurat_sub_split, FUN = function(x){
    x <- SCTransform(x,
                     vst.flavor = "v2",
                     variable.features.n = 3000,
                     return.only.var.genes = FALSE)
    return(x)
  })
  seu_sub_m <- merge(seurat_sub_split[[1]], y = seurat_sub_split[[2]])
  
  # prep for DE analysis and visualization
  seu_sub_m %<>% PrepSCTFindMarkers(assay = "SCT")
  
  # set variable features
  VariableFeatures(seu_sub_m) <- c(VariableFeatures(seurat_sub_split[[1]]),
                                   VariableFeatures(seurat_sub_split[[2]]))
  
  # run PCA
  seu_sub_m %<>% RunPCA()
  
  # save the seurat object
  saveRDS(seu_sub_m, file = file.path(output_dir, paste0(paste(population, collapse = "_"), "_seu_sub_m", ".rds")))
  
  # print message
  cat("Normalization, merge and runPCA completed: ", paste(population, collapse = ", "), ". Resulting seurat object saved to: ", output_dir, paste0(population, "_seu_sub_m", ".rds"))
  
  # test different PC and resolution
  i_range <- c(25, 30, 35)
  j_range <- seq(from = 0.2, to = .8, by = 0.2)
  plots_list <- list()
  for (i in i_range){
    for (j in j_range){
      set.seed(10086)
      
      seu_sub_m.c <- seu_sub_m %>%
        FindNeighbors(dims = 1:i) %>%
        FindClusters(resolution = j) %>%
        RunUMAP(dims = 1:i)
      
      plot_obj <- DimPlot(seu_sub_m.c, reduction = "umap", label = TRUE, repel = TRUE) +
        labs(caption = paste0("PC=", as.character(i), " res=", as.character(j)))
      plots_list[[paste0("PC", i, "res", j)]] <- plot_obj
    }
  }
  png(filename = file.path(output_dir, paste0(population, "_seu_sub_m", "pc_res_umaps.png")), width = 600*length(i_range), height = 500*length(j_range))
  wrap_plots(plots_list, ncol = length(i_range))
  dev.off()
  
  # print message
  cat("UMAP plots of combination of PCs and resolutions saved to: ", output_dir, paste0(paste(population, collapse = "_"), "_seu_sub_m", "pc_res_umaps.png"))
}

## ===== define options for the script ===== ##
description_text <- "This script is to run basic subcluster analysis on a Seurat object with specified cell populations. 
Note that this script is tailored to two samples, no integration (merge), Seurat object RData variable named merged.c.

Usage: Rscript sub_analysis.R -i /path/to/input.RData -p monocyte,macrophage -o /path/to/output_dir/"

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "Path to the input RData file"),
  make_option(c("-p", "--population"), type = "character", default = NULL, help = "Cell populations to subset separated by comma"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL, help = "Directory for the output files")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = description_text)
opt <- parse_args(opt_parser)

## ===== check if input file and populations are provided ===== ##
if (is.null(opt$input) || is.null(opt$population)) {
  print_help(opt_parser)
  stop("Input Seurat object RData file and cell populations are required.", call. = FALSE)
}

## ===== load the input file and run the sub analysis ===== ##
load(opt$input)
seurat_obj <- merged.c

# parse the populations argument into a vector
population <- unlist(strsplit(opt$population, ","))

# call the function and run the analysis
sub_analysis(seurat_obj, population, opt$output_dir)


#note that this pipeline is tailored to two samples, no integration (merge), seurat variable name merged.c
