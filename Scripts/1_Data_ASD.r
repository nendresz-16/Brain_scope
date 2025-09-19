library(readr)
library(dplyr)
library(purrr)

# Extract tar.gz into the current working directory
#untar("/project/rrg-shreejoy/nendresz/Brain_scope/Cell_metadata.tar.gz", 
      #exdir = "/project/rrg-shreejoy/nendresz/Brain_scope/")

#unzip("/project/rrg-shreejoy/nendresz/Brain_scope/snrna_expr_matrices.zip",
      #exdir = "/project/rrg-shreejoy/nendresz/Brain_scope/snrna_expr_matrices")


# Now check what files got extracted
#list.files("/project/rrg-shreejoy/nendresz/Brain_scope/")


#Read sample metadata
meta <- read_tsv("/project/rrg-shreejoy/nendresz/Brain_scope/PEC2_sample_metadata.txt")

# Filter for ASD cases
meta_asd <- meta %>%
  filter(Disorder == "ASD")

# Preview
meta_asd %>% count(Cohort)

#Cohorts w ASD = DevBrain, UCLA-ASD, Velmeshev_et_al

/project/rrg-shreejoy/nendresz/Brain_scope/snrna_expr_matrices/snrna_expr_matrices/DevBrain

devbrain_files <- list.files(
  "/project/rrg-shreejoy/nendresz/Brain_scope/snrna_expr_matrices/snrna_expr_matrices/DevBrain",
  pattern = "annotated_matrix.txt.gz$",
  full.names = TRUE
)

/project/rrg-shreejoy/nendresz/Brain_scope/snrna_expr_matrices/snrna_expr_matrices/UCLA-ASD
/project/rrg-shreejoy/nendresz/Brain_scope/snrna_expr_matrices/snrna_expr_matrices/Velmeshev_et_al

#cell metadata
Dev <- read_tsv("/project/rrg-shreejoy/nendresz/Brain_scope/Cell_metadata/DevBrain_cell_metadata.tsv")
UCLA <- read_tsv("/project/rrg-shreejoy/nendresz/Brain_scope/Cell_metadata/UCLA-ASD_cell_metadata_mismatches_removed.tsv")
Krig <- read_tsv("/project/rrg-shreejoy/nendresz/Brain_scope/Cell_metadata/Kriegstein_cell_metadata.tsv")



HSB189 <- read_tsv("/project/rrg-shreejoy/nendresz/Brain_scope/snrna_expr_matrices/snrna_expr_matrices/DevBrain/HSB189_DFC-annotated_matrix.txt.gz")

HSB189_meta <- Dev%>%filter(individualID == "HSB189_DFC")


HSB189_meta_df <- as.data.frame(HSB189_meta)   # convert tibble → data.frame
rownames(HSB189_meta_df) <- HSB189_meta_df$barcodekey



HSB189_mat <- as.data.frame(HSB189)

# Set rownames = genes
rownames(HSB189_mat) <- HSB189_mat$featurekey
HSB189_mat$featurekey <- NULL

# Replace column names with barcodekeys
colnames(HSB189_mat) <- rownames(HSB189_meta_df)




process_donor <- function(file_path, meta) {
  # Get donor ID from filename
  donor_id <- sub("-annotated_matrix.txt.gz", "", basename(file_path))
  
  # Expression matrix
  expr <- read_tsv(file_path)
  
  # Metadata for that donor
  meta_donor <- meta %>% filter(individualID == donor_id)
  meta_df <- as.data.frame(meta_donor)
  rownames(meta_df) <- meta_df$barcodekey
  
  # Clean expression
  expr_df <- as.data.frame(expr)
  rownames(expr_df) <- expr_df$featurekey
  expr_df$featurekey <- NULL
  colnames(expr_df) <- rownames(meta_df)
  
  list(expr = expr_df, meta = meta_df, donor_id = donor_id)
}


devbrain_files <- list.files(
  "/project/rrg-shreejoy/nendresz/Brain_scope/snrna_expr_matrices/snrna_expr_matrices/DevBrain",
  pattern = "annotated_matrix.txt.gz$",
  full.names = TRUE
)

devbrain_list <- lapply(devbrain_files, process_donor, meta = Dev)
names(devbrain_list) <- sapply(devbrain_list, function(x) x$donor_id)


# Extract expression matrices
expr_list <- lapply(devbrain_list, function(x) x$expr)

# Check gene order is consistent
all(sapply(expr_list, function(x) identical(rownames(x), rownames(expr_list[[1]]))))
# should be TRUE

# Combine (cbind keeps genes in rows, cells in columns)
expr_list2 <- lapply(expr_list, function(x) {
  m <- as.matrix(x)
  rownames(m) <- rownames(x)   # gene names
  m
})

expr_all <- do.call(cbind, expr_list2)
dim(expr_all)


meta_all <- do.call(rbind, lapply(devbrain_list, function(x) x$meta))

#Save 
saveRDS(expr_all,"/project/rrg-shreejoy/nendresz/Brain_scope/DevBrain_counts.rds")



###### UCLA #######


process_donor <- function(file_path, meta) {
  # Get donor ID from filename (before "-annotated_matrix.txt.gz")
  donor_id <- sub("-annotated_matrix.txt.gz", "", basename(file_path))
  
  # Expression matrix
  expr <- read_tsv(file_path)
  
  # Metadata for that donor
  meta_donor <- meta %>% filter(individualID == donor_id)
  meta_df <- as.data.frame(meta_donor)
  rownames(meta_df) <- meta_df$barcodekey
  
  # Clean expression
  expr_df <- as.data.frame(expr)
  rownames(expr_df) <- expr_df$featurekey
  expr_df$featurekey <- NULL
  colnames(expr_df) <- rownames(meta_df)
  
  list(expr = expr_df, meta = meta_df, donor_id = donor_id)
}

ucla_files <- list.files(
  "/project/rrg-shreejoy/nendresz/Brain_scope/snrna_expr_matrices/snrna_expr_matrices/UCLA-ASD",
  pattern = "annotated_matrix.txt.gz$",
  full.names = TRUE
)

ucla_list <- lapply(ucla_files, process_donor, meta = UCLA)
names(ucla_list) <- sapply(ucla_list, function(x) x$donor_id)


# Expression
expr_list <- lapply(ucla_list, function(x) x$expr)
expr_list2 <- lapply(expr_list, function(x) {
  m <- as.matrix(x)
  rownames(m) <- rownames(x)
  m
})

stopifnot(all(sapply(expr_list2, function(m) identical(rownames(m), rownames(expr_list2[[1]])))))

expr_all <- do.call(cbind, expr_list2)

# Metadata
meta_all <- do.call(rbind, lapply(ucla_list, function(x) x$meta))


#Save 
saveRDS(expr_all,"/project/rrg-shreejoy/nendresz/Brain_scope/UCLA_counts.rds")




######### Velmshev #############

process_donor <- function(file_path, meta) {
  # Get donor ID from filename (before "-annotated_matrix.txt.gz")
  donor_id <- sub("-annotated_matrix.txt.gz", "", basename(file_path))
  
  # Expression matrix
  expr <- read_tsv(file_path)
  
  # Metadata for that donor
  meta_donor <- meta %>% filter(individualID == donor_id)
  meta_df <- as.data.frame(meta_donor)
  rownames(meta_df) <- meta_df$barcodekey
  
  # Clean expression
  expr_df <- as.data.frame(expr)
  rownames(expr_df) <- expr_df$featurekey
  expr_df$featurekey <- NULL
  colnames(expr_df) <- rownames(meta_df)
  
  list(expr = expr_df, meta = meta_df, donor_id = donor_id)
}

# List all Velmeshev donor files
vel_files <- list.files(
  "/project/rrg-shreejoy/nendresz/Brain_scope/snrna_expr_matrices/snrna_expr_matrices/Velmeshev_et_al",
  pattern = "annotated_matrix.txt.gz$",
  full.names = TRUE
)

# Run process_donor for each file
vel_list <- lapply(vel_files, process_donor, meta = Krig)
names(vel_list) <- sapply(vel_list, function(x) x$donor_id)

# Expression
expr_list <- lapply(vel_list, function(x) x$expr)
expr_list2 <- lapply(expr_list, function(x) {
  m <- as.matrix(x)
  rownames(m) <- rownames(x)
  m
})


# Merge into one giant matrix
expr_all <- do.call(cbind, expr_list2)

# Metadata
meta_all <- do.call(rbind, lapply(vel_list, function(x) x$meta))

expr_all[1:5, 1:5]

saveRDS(expr_all,"/project/rrg-shreejoy/nendresz/Brain_scope/Velmeshev_counts.rds")

