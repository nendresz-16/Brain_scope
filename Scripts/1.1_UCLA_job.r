
###### UCLA #######
UCLA <- read_tsv("/project/rrg-shreejoy/nendresz/Brain_scope/Cell_metadata/UCLA-ASD_cell_metadata_mismatches_removed.tsv")

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


