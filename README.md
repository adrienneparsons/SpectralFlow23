# SpectralFlow23
Harmony integration accompanying Spasic et al. 2023

This repository contains a script and files used to implement Harmony integration of bone marrow samples collected for Spasic et al. 2023. The script is designed to generate new UMAP coordinates for each cell following Harmony integration, and then it adds those coordinates to the data and exports a new .fcs file. 

## Set up for script
High-dimensional spectral flow cytometry data was collected on a Cytek Aurora flow cytometer. FCS files were exported from SpectroFlo and imported into the cytometry analysis platform OMIQ.ai ([app.omiq.ai](https://app.omiq.ai/)). In OMIQ, quality filters were applied as described in the accompanying manuscript (Spasic et al, under review). Cell types were identified in OMIQ using conventional two-dimensional gating. Then, raw fluorescence data was exported for all parameters for every high quality, singlet, live, CD45+ cell.

## Running the script
These csv files are loaded into R as described in "Harmony integration of BM samples.R". Parameters that are not desired to be included in the UMAP coordinate generation (i.e. Time, FSC/SSC, Viability, and autofluorescence) are removed, and each csv is turned into an indiviual Seurat object with the name of the sample being the Project name. The Seurat objects are then merged into a single object.

From there, the merged Seurat object has its data scaled and dimensionality reduction in performed. Harmony is then run on the Seurat object to regress out the sample-specific batch effects. This results in the following UMAP (may vary slighly depending on software version):

![alt text](/integrated_umap.png "Integrated UMAP")

The UMAP coordinates for the entire Seurat object are extracted and subsetted to be the coordinates of the cells in each of the three samples, and then the UMAP coordinates are appended to the raw fluorescence values.

OMIQ does not recognize the fluorophore names as they are exported, so a final step before writing a new .fcs file from the UMAP coordinate-appended data frames is to rename the columns such that OMIQ will recognize the names upon re-uploading. Finally, new .fcs files with UMAP coordinates are written.

These new .fcs files can be uploaded into OMIQ, and using the "Figure" function, scatterplots can be drawn with the saved UMAP coordinates.
