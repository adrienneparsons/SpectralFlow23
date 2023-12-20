# Harmony integration of flow cytometry samples and LISI analysis
# Date: 12/20/2023

# Prerequisite packages
library(Seurat)
library(harmony)
library(flowCore)
library(Matrix)
library(ggplot2)
library(ggrastr)
library(dplyr)
library(tidyr)
library(forcats)

# Harmony integration of BM flow samples
# Load in the samples
samp1 <- read.csv("/Users/addie/Dropbox (Partners HealthCare)/Human-PBMC_Methods Paper/Data/Harmony integration of BM samples/20230522_BMsample_fromOMIQ_pre-scaling/20230522_BMsample_fromOMIQ_pre-scaling_livecells/020222 VAN GALEN PANEL unmixed-ALL_Unmixed.csv")
samp2 <- read.csv("/Users/addie/Dropbox (Partners HealthCare)/Human-PBMC_Methods Paper/Data/Harmony integration of BM samples/20230522_BMsample_fromOMIQ_pre-scaling/20230522_BMsample_fromOMIQ_pre-scaling_livecells/032822 VAN GALEN PANEL-ALL-Unmixed.csv")
samp3 <- read.csv("/Users/addie/Dropbox (Partners HealthCare)/Human-PBMC_Methods Paper/Data/Harmony integration of BM samples/20230522_BMsample_fromOMIQ_pre-scaling/20230522_BMsample_fromOMIQ_pre-scaling_livecells/07142022 BM ONLY.csv")

# Prepare for Harmony by merging the data and making a Seurat object
seu.ls <- vector(mode = "list", length = 3)

for(sample in c("samp1", "samp2", "samp3")){
  print(sample)
  df <- get(sample)
  
  # Remove the parameters that are not needed for generating the UMAP
  df[,c(1:7, 38, 41)] <- NULL

  # Name each cell
  rownames(df) <- paste0(sample, "_", rownames(df))
  
  # Create Seurat object
  seu <- CreateSeuratObject(counts = t(df), 
                            project = sample)
  
  # Add it to the list
  seu.ls[[which(sample == c("samp1", "samp2", "samp3"))]] <- seu
}

# Merge
seu.all <- merge(x = seu.ls[[1]], y = seu.ls[[2]])
seu.all <- merge(x = seu.all, y = seu.ls[[3]])  

# Make the data slot the same as the uploaded MFI data
seu.all@assays$RNA@data <- seu.all@assays$RNA@counts

# The values from scaling in OMIQ are already normalized, and variable features
# Are not necessary to compute (there are only 33 features to compare)
# Generate dimensionality reduction and view the PCs that contribute most
seu.all <- ScaleData(seu.all)
features <- rownames(seu.all)
seu.all <- RunPCA(seu.all, features = features, approx = FALSE)
ElbowPlot(seu.all)

# Run Harmony to batch correct by donor
seu.all.2 <- RunHarmony(object = seu.all, 'orig.ident')

# seu.all.harmony is running the UMAP with the harmony reduction
seu.all.harmony <- RunUMAP(seu.all.2, reduction = "harmony", dims = 1:20)

# To compare UMAP coordinates without integration, you can generate a UMAP
# without the Harmony integration
# Based on the above elbow plot, find neighbors based on first 20 PCs, then generate a UMAP
seu.all <- FindNeighbors(seu.all, dims = 1:20, nn.method = "annoy")
seu.all.umap <- RunUMAP(seu.all, dims = 1:20)

# Visualize the data with and without Harmony integration
DimPlot(seu.all.umap, reduction = "umap", group.by = "orig.ident", shuffle = T, raster = T)+ theme(aspect.ratio = 1)
DimPlot(seu.all.harmony, reduction = "umap", group.by = "orig.ident", shuffle = TRUE, raster = T)+ theme(aspect.ratio = 1)

# Save the Harmony coordinates from the merged seurat object
UMAP_coords <- as.data.frame(seu.all.harmony[["umap"]]@cell.embeddings)

# Assign the coordinates to new variables for each sample
samp1_coords <- UMAP_coords[1:nrow(samp1),]
samp2_coords <- UMAP_coords[(nrow(samp1)+1):(nrow(samp1)+nrow(samp2)),]
samp3_coords <- UMAP_coords[(nrow(samp1)+nrow(samp2)+1):nrow(UMAP_coords),]

# Add the UMAP coordinates to the data
samp1 <- cbind(samp1, samp1_coords)
samp2 <- cbind(samp2, samp2_coords)
samp3 <- cbind(samp3, samp3_coords)


# Make the sample data matrices
samp1a <- as.matrix(samp1)
samp2a <- as.matrix(samp2)
samp3a <- as.matrix(samp3)
rownames(samp1a) <- NULL
rownames(samp2a) <- NULL
rownames(samp3a) <- NULL

# Re-format the columnn names so OMIQ will read them
colnames <- c("Time", "SSC-H", "SSC-A",	"FSC-H", "FSC-A", "SSC-B-H", "SSC-B-A",
                     "BV421-A",
                     "cFluor V450-A",	"BV480-A",	"BV510-A",
                     "cFluor V547-A",	"BV570-A",	"BV605-A",	
                     "BV650-A",	"BV711-A",	"BV750-A",	"BV785-A",
                     "BB515-A",	"cFluor B548-A","NovaFL Blue 610-70S-A",	"cFluor B677-A",
                     "BB700-A",	"BB755-A",	"cFluor BYG575-A",
                     "cFluor YG584-A",	"cFluor BYG610-A", "cFluor YG610-A",
                     "PE-Fire 640-A",	"cFluor BYG667-A",	"cFluor BYG710-A", "PerCP eFl 710-A",
                     "cFluor BYG750-A",	"cFluor BYG781-A",	"cFluor R659-A",	"cFluor R685-A",
                     "cFluor R720-A",	"ViaDye Red-A", "cFluor R780-A",	"cFluor R840-A", "AF-A",
                     "UMAP_1", "UMAP_2")

colnames(samp1a) <- colnames
colnames(samp2a) <- colnames
colnames(samp3a) <- colnames

# Write the new .fcs files for the three BM samples
fframe_s1 <- flowFrame(exprs = samp1a)
fframe_s2 <- flowFrame(exprs = samp2a)
fframe_s3 <- flowFrame(exprs = samp3a)

# Save
setwd("/Users/addie/Dropbox (Partners HealthCare)/Human-PBMC_Methods Paper/Data/Harmony integration of BM samples")
write.FCS(fframe_s1,"20230522_BMSample1.fcs", what="numeric", delimiter = "\\")
write.FCS(fframe_s2,"20230522_BMSample2.fcs", what="numeric", delimiter = "\\")
write.FCS(fframe_s3,"20230522_BMSample3.fcs", what="numeric", delimiter = "\\")

######################################################################
#LISI proportion plotting, pre-Harmony

#For bone marrow, calculate the proportion of pre-Harmony LISI scores 
# that fall above n-1 (2)

# First, add the labels to the results
res_BM_preharm$Sample <- all_sample_preharm_labels$Label

# Calculate the proportions
props_list_BM <- list()
for(donor in unique(res_BM_preharm$Sample)){
  df <- res_BM_preharm[res_BM_preharm$Sample == donor,] 
  props <- sum(df$Label >= length(unique(res_BM_preharm$Sample))-1) / nrow(df)
  props_list_BM <- append(props_list_BM, props)
}


#Generate a pre-Harmony UMAP for each PBMC UMAP, to be used later
# Argument passed to function is the directory where pre-Harmony UMAP coordinates
# are found
preharm_umap <- function(dir){
setwd(dir)
  
# List the files and remove the extra one
files <- list.files()
files <- files[-1]

# Read in the CSV files and add the resulting data frames to a list
my_data <- list()
my_data <- lapply(files, read.csv)

# Format the data frames so the sample ID is annotated
for(df in 1:length(my_data)){
  data <- my_data[[df]]
  data$Sample <- paste0("samp", df)
  my_data[[df]] <- data
}

# Make a final concatenated data frame of all 3 samples
final_df <- rbind(my_data[[1]], my_data[[2]])
final_df <- rbind(final_df, my_data[[3]])

# Pull just the UMAP coordinates from the data
preharm_UMAP_coords <- final_df[, c("umap_1", "umap_2", "Sample")]

# To make the UMAPs more visually meaningful and to prevent one cell type covering another,
#shuffle the order of the rows
shuf_index <- sample(1:nrow(preharm_UMAP_coords), replace = F)
preharm <- preharm_UMAP_coords[shuf_index,]

# generate a UMAP, rasterize points, and save
p <- ggplot(preharm, aes(x = umap_1, y = umap_2, color = Sample))+
  geom_point(size = 0.5, alpha = 0.75)+
  theme(aspect.ratio = 1)

rasterize(p, layers = 'Point')

ggsave("preharmony_umap_by_sample.pdf", p, device = "pdf")}

# Generate the pre-Harmony UMAPs for each PBMC UMAP
MND_lymph_umap <- preharm_umap("/Users/addie/Downloads/MND_Lymphocytes-UMAP_coordinates")
MND_allcells_umap <- preharm_umap("/Users/addie/Downloads/MND_AllCells-UMAP_coordinates")
TB_lymphs_umap <- preharm_umap("/Users/addie/Downloads/TB_Lymphocyte-UMAP_coordinates 2")
MND_Monos_umap <- preharm_umap("/Users/addie/Downloads/MND_Monocytes-UMAP_coordinates")

# Function for calculating proportion of high LISI cells per UMAP
# Argument passed to function is the directory where the data is
LISI_props <- function(dir){
  setwd(dir)
  
  # Get the files with the cell-specific data
  files <- list.files(dir)
  files <- files[-1]
  
  # Read eaach data file in and append it to a list
  my_data <- list()
  my_data <- lapply(files, read.csv)
  
  # Make a list of the number of rows in each data frame (nummber of cells)
  n.list <- list()
  for(df in 1:length(my_data)){
    n.list <- append(n.list, nrow(my_data[[df]]))
  }
  
  # Make the list a numeric vector
  n.list <- unlist(n.list)
  
  # Concatenate the data into a single data frame
  final_df <- rbind(my_data[[1]], my_data[[2]])
  final_df <- rbind(final_df, my_data[[3]])
  
  # Get the UMAP coordinates from the data
  preharm_UMAP_coords <- final_df[, c("umap_1", "umap_2")]
  
  # Get the labels needed for LISI analysis
  preharm_labels <- data.frame("Label" = c(rep("samp1", n.list[1]), rep("samp2", n.list[2]), 
                                           rep("samp3", n.list[3])))
  
  # Run LISI analysis
  res_preharm <- compute_lisi(preharm_UMAP_coords, preharm_labels,"Label")
  
  # Add the labels to the LISI data
  res_preharm$Sample <- preharm_labels$Label
  
  # Calculate the proportion of high-LISI cells per sample
  props_list <- list()
  for(donor in unique(res_preharm$Sample)){
    df <- res_preharm[res_preharm$Sample == donor,] 
    props <- sum(df$Label >= length(unique(res_preharm$Sample))-1) / nrow(df)
    props_list <- append(props_list, props)
  }
  
  return(props_list)}

# Run the LISI proportion calculation function for each PBMC UMAP
MND_lymph <- LISI_props("/Users/addie/Downloads/MND_Lymphocytes-UMAP_coordinates")
MND_allcells <- LISI_props("/Users/addie/Downloads/MND_AllCells-UMAP_coordinates")
TB_lymphs <- LISI_props("/Users/addie/Downloads/TB_Lymphocyte-UMAP_coordinates 2")
MND_Monos <- LISI_props("/Users/addie/Downloads/MND_Monocytes-UMAP_coordinates")

# New function for calculating high LISI proportions excluding donor with
# a lymphoma

# argument passed to function is the directory where the data is
nolymphoma <- function(dir){
  setwd(dir)
  
  # Get the files names for all of the samples except lymphoma donor
  files <- list.files(dir)
  files <- files[-1]
  files <- files[-grep("3146", files)]
  
  # Read in the other files and append the data frame to a list
  my_data <- list()
  my_data <- lapply(files, read.csv)
  
  # Get the number of rows per data frame in the list
  n.list <- list()
  for(df in 1:length(my_data)){
    n.list <- append(n.list, nrow(my_data[[df]]))
  }
  
  # Make the list of nrows a numeric vector
  n.list <- unlist(n.list)
  
  # Concatenate the data
  final_df <- rbind(my_data[[1]], my_data[[2]])
  
  # Get the UMAP coordinates from the data
  preharm_UMAP_coords <- final_df[, c("umap_1", "umap_2")]
  
  # Generate the labels for LISI
  preharm_labels <- data.frame("Label" = c(rep("samp1", n.list[1]), rep("samp2", n.list[2])))
  
  # Compute LISI
  res_preharm <- compute_lisi(preharm_UMAP_coords, preharm_labels,"Label")
  
  # Add the labels to the LISI results
  res_preharm$Sample <- preharm_labels$Label
  
  #Calculate the proportion of high LISI cells for both samples
  # NOTE: the calculation is different because calculating n-1 for two samples
  # would mean every cell had a high proportion LISI
  props_list <- list()
  for(donor in unique(res_preharm$Sample)){
    df <- res_preharm[res_preharm$Sample == donor,] 
    props <- sum(df$Label > 1+ 0.3333*length(unique(res_preharm$Sample))) / nrow(df)
    props_list <- append(props_list, props)
  }
  
  return(props_list)
}

# Run the function for all PBMC UMAPs
MND_lymph_nolymphoma <- nolymphoma("/Users/addie/Downloads/MND_Lymphocytes-UMAP_coordinates")
MND_allcells_nolymphoma <- nolymphoma("/Users/addie/Downloads/MND_AllCells-UMAP_coordinates")
TB_lymphs_nolymphoma <- nolymphoma("/Users/addie/Downloads/TB_Lymphocyte-UMAP_coordinates 2")
MND_Monos_nolymphoma <- nolymphoma("/Users/addie/Downloads/MND_Monocytes-UMAP_coordinates")

# Make a data frame with all of the LISI data
final_LISIdata <- data.frame("No_Lymphoma_TB_lymphocytes" = c(NA, unlist(TB_lymphs_nolymphoma)),
                             "No_Lymphoma_MND_All_cells" = c(NA, unlist(MND_allcells_nolymphoma)),
                             "No_Lymphoma_MND_Just_lymphocytes" = c(NA, unlist(MND_lymph_nolymphoma)),
                             "No_Lymphoma_MND_Just_monocytes" = c(NA, unlist(MND_Monos_nolymphoma)),
                             "TB_Lymphocytes" = unlist(TB_lymphs),
                             "MND_All_cells" = unlist(MND_allcells),
                             "MND_Just_lymphocytes" = unlist(MND_lymph),
                             "MND_Just_monocytes" = unlist(MND_Monos),
                             "Bone_Marrow_Cells" = unlist(props_list_BM))

# Add the means to the data frame for plotting
final_LISIdata[4,5:9] <- colMeans(final_LISIdata[,5:9])
final_LISIdata[4, 1:4] <- colMeans(final_LISIdata[2:3, 1:4])
rownames(final_LISIdata) <- c("3146", "3156", "958", "Average")

# Format the data for plotting
final_LISIdata2 <- pivot_longer(final_LISIdata, colnames(final_LISIdata))
final_LISIdata2$Label <- c(rep("PBMC_3146", 8), "BMC_1", rep("PBMC_3156", 8), 
                           "BMC_2", rep("PBMC_958", 8), "BMC_3", rep("Average", 9))

# Add an annotation for whether the data is BMC, PBMC, or average
final_LISIdata2$Shape <- "PBMCs"
final_LISIdata2$Shape[final_LISIdata2$Label == "Average"] <- "Average"
final_LISIdata2$Shape[c(9, 18, 27)] <- "BMC"

# Plot the data
plot <- ggplot(na.omit(final_LISIdata2), aes(x = fct_relevel(name,
                                                             c("Bone_Marrow_Cells",
                                                               "TB_Lymphocytes",
                                                               "MND_All_cells",
                                                               "MND_Just_lymphocytes",
                                                               "MND_Just_monocytes",
                                                               "No_Lymphoma_TB_lymphocytes",
                                                               "No_Lymphoma_MND_All_cells",
                                                               "No_Lymphoma_MND_Just_lymphocytes",
                                                               "No_Lymphoma_MND_Just_monocytes")), 
                                             y = value, color = fct_relevel(Label,
                                                                            "PBMC_3146",
                                                                            "BMC_1",
                                                                            "PBMC_3156",
                                                                            "BMC_2",
                                                                            "PBMC_958",
                                                                            "BMC_3",
                                                                            "Average"),
                                             shape = Shape))+
  geom_point(size = 4)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(values = c("green3", "green3", "orange", "orange", "deepskyblue2",
                                "deepskyblue2", "red"))+
  scale_shape_manual(values = c(16, 17, 15))+
  theme(aspect.ratio = 1)+
  xlab("Sample")+
  ylab("Proportion high LISI cells")+
  theme(legend.title = element_blank())+
  theme(axis.text = element_text(size = 16))+
  theme(axis.title = element_text(size = 16))+
  theme(legend.text=element_text(size=16))

# Save the plot
setwd("/Users/addie/desktop")
ggsave("highLISI_preharmony.pdf", plot, device = "pdf", height = 8, width = 8)
