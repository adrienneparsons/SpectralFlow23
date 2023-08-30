# Adrienne Parsons
# May 22, 2023


# Prerequisites
library(Seurat)
library(harmony)
library(flowCore)
library(Matrix)
library(ggplot2)
library(ggrastr)

# Harmony integration of BM flow samples

samp1 <- read.csv("~/Data/Harmony integration of BM samples/20230522_BMsample_fromOMIQ_pre-scaling/20230522_BMsample_fromOMIQ_pre-scaling_livecells/020222 VAN GALEN PANEL unmixed-ALL_Unmixed.csv")
samp2 <- read.csv("~/Data/Harmony integration of BM samples/20230522_BMsample_fromOMIQ_pre-scaling/20230522_BMsample_fromOMIQ_pre-scaling_livecells/032822 VAN GALEN PANEL-ALL-Unmixed.csv")
samp3 <- read.csv("~/Data/Harmony integration of BM samples/20230522_BMsample_fromOMIQ_pre-scaling/20230522_BMsample_fromOMIQ_pre-scaling_livecells/07142022 BM ONLY.csv")

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
  seu.ls <- append(seu.ls, seu)
}

# Merge
seu.all <- merge(x = seu.ls[[4]], y = seu.ls[[5]])
seu.all <- merge(x = seu.all, y = seu.ls[[6]])  

# Make the data slot the same as the uploaded MFI data
seu.all@assays$RNA@data <- seu.all@assays$RNA@counts

# The values from scaling in OMIQ are already normalized, and variable features
# Are not necessary to compute (there are only 33 features to compare)
# Generate dimensionality reeduction and view the PCs that contribute most
seu.all <- ScaleData(seu.all)
features <- rownames(seu.all)
seu.all <- RunPCA(seu.all, features = features, approx = FALSE)
ElbowPlot(seu.all)

# Run Harmony to batch correct by donor (currently ignored because of runtimes)

# seu.all.2 is a test of running Harmony based on PCA, not for generating the UMAP
seu.all.2 <- RunHarmony(object = seu.all, reduction = "pca", 'orig.ident')

# seu.all.harmony is running the UMAP with the harmony reduction
seu.all.harmony <- RunUMAP(seu.all.2, reduction = "harmony", dims = 1:20)

# Visualize the data with Harmony integration
DimPlot(seu.all.harmony, reduction = "umap", group.by = "orig.ident", shuffle = TRUE, raster = T)+ theme(aspect.ratio = 1)

# Save the Harmony coordinates from the merged seurat object
UMAP_coords <- as.data.frame(seu.all.harmony[["umap"]]@cell.embeddings)

# Assign the coordinates to new variables for each sample
samp1_coords <- UMAP_coords[1:60534,]
samp2_coords <- UMAP_coords[60535:144044,]
samp3_coords <- UMAP_coords[144045:298197,]

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


setwd("~/Data/Harmony integration of BM samples")
write.FCS(fframe_s1,"20230522_BMSample1.fcs", what="numeric", delimiter = "\\")
write.FCS(fframe_s2,"20230522_BMSample2.fcs", what="numeric", delimiter = "\\")
write.FCS(fframe_s3,"20230522_BMSample3.fcs", what="numeric", delimiter = "\\")

##################################################################################