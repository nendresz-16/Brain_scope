
setwd("P1_Brain_scope")

#Load packages 
library(readxl)
library(readr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)

meta <- read_tsv("/project/rrg-shreejoy/nendresz/Brain_scope/PEC2_sample_metadata.txt")
meta <- as.data.frame(meta)
meta <- meta %>%
  dplyr::rename(
    Diagnosis = Disorder,
    Sex       = Biological_Sex,
    Age       = Age_death
  )

Dev  <- readRDS("/scratch/nendresz/P1_Brain_scope/Files/LT_DevBrain_Supertypes.rds")
Velm <- readRDS("/scratch/nendresz/P1_Brain_scope/Files/LT_Velmeshev_Supertypes.rds")
UCLA <- readRDS("/scratch/nendresz/P1_Brain_scope/Files/LT_UCLA_Supertypes.rds")

#Dev 


# make sure rownames of donor metadata are donor IDs
rownames(meta) <- meta$Individual_ID

# for each cell, pull the donor metadata based on its "individualID"
donor_info <- meta[Dev$individualID, ]

# add these donor-level columns to the Seurat metadata
Dev <- AddMetaData(Dev, donor_info)

saveRDS(Dev, "Files/LT_DevBrain_Final.rds" )

#Velm


# make sure rownames of donor metadata are donor IDs
rownames(meta) <- meta$Individual_ID

# for each cell, pull the donor metadata based on its "individualID"
donor_info <- meta[Velm$individualID, ]

# add these donor-level columns to the Seurat metadata
Velm <- AddMetaData(Velm, donor_info)

saveRDS(Velm, "Files/LT_Velmeshev_Final.rds" )



#UCLA


# make sure rownames of donor metadata are donor IDs
rownames(meta) <- meta$Individual_ID

# for each cell, pull the donor metadata based on its "individualID"
donor_info <- meta[UCLA$individualID, ]

# add these donor-level columns to the Seurat metadata
UCLA <- AddMetaData(UCLA, donor_info)

saveRDS(UCLA, "Files/LT_UCLA_Final.rds" )
