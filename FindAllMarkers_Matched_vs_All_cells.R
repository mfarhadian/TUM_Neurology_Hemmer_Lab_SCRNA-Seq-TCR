


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(openxlsx)
  library(future)
  library(forcats)
  library(stringr)
})

# ----------------------- Paths ----------------------------------------------
rds_path   <- "/data_storage/Final_all_combo_final.rds"
xlsx_path  <- "/data_storage/TOP3 Motifs_Expression profiles and Annotations.xlsx"
xlsx_sheet <- "All_paires_133_Motifes"

paired_xlsx <- "/data_storage/TCR_Paired_chains_Final_Combined.xlsx"

outdir <- "/data_storage/1_Matched_vs_All_Cells/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

xlsx_out        <- file.path(outdir, "DE_PAIREDONLY_MatchedMotif_vs_AllCells.xlsx")
tsv_out         <- file.path(outdir, "DE_PAIREDONLY_MatchedMotif_vs_AllCells.tsv")
csv_out         <- file.path(outdir, "DE_PAIREDONLY_MatchedMotif_vs_AllCells.csv")
png_out         <- file.path(outdir, "Volcano_PAIREDONLY_MatchedMotif_vs_AllCells.png")
summary_txt     <- file.path(outdir, "SUMMARY_PAIREDONLY_MatchedMotif_vs_AllCells.txt")
miss_paired_csv <- file.path(outdir, "PairedTCR_not_found_in_all_combo.csv")

# ----------------------- (Optional) user definitions -------------------------
# Keeping these here in case you want to later restrict ident.1 to CD4 only again.
cd4_subtypes <- c(
  "CD4+ Central Memory",
  "CD4+ Effector Memory",
  "Tcm/ Naive CD4+ Tcells",
  "CD4+ Naive T cells",
  "Regulatory T cells"
)

# ----------------------- DE parameters ---------------------------------------
logfc_thresh        <- 0.00
min_pct             <- 0.05
alpha               <- 0.05
label_topN          <- 20
min_abs_fc_for_sig  <- 0.25

# ----------------------- FUTURE settings -------------------------------------
plan(sequential)
options(future.globals.maxSize = 50 * 1024^3)

# ----------------------- Helpers ---------------------------------------------
norm_barcode <- function(x){
  x <- str_trim(as.character(x))
  x <- str_replace(x, "-+$", "-1")
  x <- ifelse(str_detect(x, "-\\d+$"), x, paste0(x, "-1"))
  x
}

norm_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))

pick_col <- function(df, candidates, what = "column") {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    stop(
      "Could not find required ", what, ". Tried: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns include: ",
      paste(head(colnames(df), 60), collapse = ", "),
      if (ncol(df) > 60) " ... (truncated)" else ""
    )
  }
  hit[1]
}

safe_neglog10 <- function(p) {
  out <- rep(NA_real_, length(p))
  nz  <- which(p > 0 & is.finite(p))
  z   <- which(p == 0 | !is.finite(p))
  out[nz] <- -log10(p[nz])
  if (length(z) > 0) {
    pmin <- suppressWarnings(min(p[nz], na.rm = TRUE))
    base_eps <- ifelse(is.finite(pmin) && pmin > 0, pmin/10, 1e-300)
    ranks <- rank(seq_along(z), ties.method = "first")
    eps   <- base_eps * (ranks / (length(z) + 1))
    out[z] <- -log10(eps)
  }
  out
}

# ----------------------- Load Seurat -----------------------------------------
message("Reading Seurat object...")
all_combo <- readRDS(rds_path)
stopifnot(is(all_combo, "Seurat"))

meta0 <- all_combo@meta.data
meta0$.cell <- rownames(meta0)

# ----------------------- Load motif Excel ------------------------------------
message("Reading motif Excel...")
targets_xlsx <- openxlsx::read.xlsx(xlsx_path, sheet = xlsx_sheet)
stopifnot(ncol(targets_xlsx) >= 1)
colnames(targets_xlsx)[1] <- "Final_merged"
targets_xlsx <- targets_xlsx %>% distinct(Final_merged, .keep_all = TRUE)

target_keys <- as.character(targets_xlsx$Final_merged)

# Optional Motif_Number mapping (if present)
has_motif_number <- "Motif_Number" %in% colnames(targets_xlsx)
motif_map <- NULL
if (has_motif_number) {
  motif_map <- setNames(as.character(targets_xlsx$Motif_Number), target_keys)
}

# ----------------------- Load paired TCR Excel --------------------------------
message("Reading Paired_TCR Excel...")
Paired_TCR <- openxlsx::read.xlsx(paired_xlsx)

paired_batch_col <- pick_col(
  Paired_TCR,
  candidates = c(
    "Batch_Pool_sample", "Batch_Pool_info",
    "Contig_Sample", "Airr_Sample",
    "Contig_Sample_name_TRA_1", "Contig_Sample_name_TRB_1"
  ),
  what = "Paired_TCR batch/sample column"
)

paired_barcode_col <- pick_col(
  Paired_TCR,
  candidates = c("cell_join", "Airr_barcode", "Contig_barcode", "barcode", "Airr_Barcode", "Contig_Barcode"),
  what = "Paired_TCR barcode column"
)

paired_targets <- Paired_TCR %>%
  transmute(
    Batch_Pool_sample = as.character(.data[[paired_batch_col]]),
    cell_join         = norm_barcode(.data[[paired_barcode_col]])
  ) %>%
  distinct()

# ----------------------- Match Paired_TCR -> all_combo ------------------------
needed_meta_cols <- c("Batch_Pool_sample", "cell_join")
if (!all(needed_meta_cols %in% colnames(meta0))) {
  stop(
    "Seurat metadata must contain columns: ",
    paste(needed_meta_cols, collapse = ", "),
    "\nAvailable columns include: ",
    paste(head(colnames(meta0), 80), collapse = ", "),
    if (ncol(meta0) > 80) " ... (truncated)" else ""
  )
}

meta_key_paired <- meta0 %>%
  transmute(
    .cell,
    Batch_Pool_sample = as.character(.data$Batch_Pool_sample),
    cell_join         = norm_barcode(.data$cell_join)
  )

paired_hits <- meta_key_paired %>%
  inner_join(paired_targets, by = c("Batch_Pool_sample","cell_join"))

paired_cells_in_allcombo <- unique(paired_hits$.cell)

paired_missed <- paired_targets %>% anti_join(meta_key_paired, by = c("Batch_Pool_sample","cell_join"))
if (nrow(paired_missed) > 0) {
  write.csv(paired_missed, miss_paired_csv, row.names = FALSE)
}

paired_cells_in_allcombo <- intersect(paired_cells_in_allcombo, colnames(all_combo))

message(sprintf("Paired targets (unique batch+barcode): %d", nrow(paired_targets)))
message(sprintf("Paired cells matched into all_combo:   %d", length(paired_cells_in_allcombo)))
if (nrow(paired_missed) > 0) message(sprintf("Paired targets NOT found:               %d (saved: %s)", nrow(paired_missed), miss_paired_csv))

if (length(paired_cells_in_allcombo) < 2) {
  stop("Too few paired cells matched in all_combo. Check Paired_TCR columns and meta keys.")
}

# ----------------------- SUBSET Seurat to PAIRED cells ONLY -------------------
all_combo_paired <- subset(all_combo, cells = paired_cells_in_allcombo)
meta <- all_combo_paired@meta.data
meta$.cell <- rownames(meta)

# ----------------------- Find Final_merged key column in PAIRED metadata -------
meta_cols_norm <- norm_name(colnames(meta))
key_idx <- which(meta_cols_norm == "finalmerged")
if (length(key_idx) == 0) stop("Could not find a 'Final_merged' equivalent in PAIRED Seurat metadata.")
meta_key_col <- colnames(meta)[key_idx[1]]
message("Using metadata key column for motifs: ", meta_key_col)

# ----------------------- Motif matching within PAIRED cells -------------------
cell_keys <- as.character(meta[[meta_key_col]])
is_match  <- cell_keys %in% target_keys

meta$Motif_Number <- NA_character_
if (has_motif_number) {
  meta$Motif_Number <- unname(motif_map[cell_keys])
}
meta$Motif_Match <- ifelse(is_match, "matched", "unmatched")

all_combo_paired@meta.data <- meta

n_paired_total          <- ncol(all_combo_paired)
n_paired_motif_matched  <- sum(meta$Motif_Match == "matched", na.rm = TRUE)
n_paired_motif_unmatched<- sum(meta$Motif_Match == "unmatched", na.rm = TRUE)

# ----------------------- Define groups: matched motifs vs ALL other cells -----
# ident.1: ALL motif-matched paired cells (no CD4 restriction)
is_matched_motif <- meta$Motif_Match == "matched"

# ident.2: ALL other paired cells (this is "all cells" within PAIRED scope)
is_bg_all_other <- !is_matched_motif

group <- rep(NA_character_, nrow(meta))
group[is_matched_motif] <- "matched_motif"
group[is_bg_all_other]  <- "bg_all_other"

meta$DE_group_custom <- group
all_combo_paired@meta.data <- meta

# Keep only the two groups (in practice this keeps ALL paired cells)
keep_cells <- rownames(meta)[!is.na(meta$DE_group_custom)]
all_combo_de <- subset(all_combo_paired, cells = keep_cells)

Idents(all_combo_de) <- factor(
  all_combo_de@meta.data$DE_group_custom,
  levels = c("bg_all_other","matched_motif")
)

n_matched <- sum(Idents(all_combo_de) == "matched_motif")
n_bg      <- sum(Idents(all_combo_de) == "bg_all_other")

message("Group counts (PAIRED-only): matched_motif vs bg_all_other")
print(table(all_combo_paired@meta.data$DE_group_custom, useNA = "ifany"))
message(sprintf("Post-filter (PAIRED-only) -> matched_motif: %d | bg_all_other: %d | total used: %d",
                n_matched, n_bg, ncol(all_combo_de)))

if (n_matched < 1) stop("No matched motif cells (PAIRED-only) - nothing to compare.")
if (n_bg < 1)      stop("No background cells (PAIRED-only) - nothing to compare.")


# ----------------------- DE (Wilcoxon) ---------------------------------------
message("Running DE (Wilcoxon) on PAIRED-only cells: matched_motif vs bg_all_other ...")

suppressPackageStartupMessages(library(presto))

de <- FindMarkers(
  all_combo_de,
  ident.1 = "matched_motif",
  ident.2 = "bg_all_other",
  test.use = "wilcox",
  logfc.threshold = logfc_thresh,
  min.pct = min_pct,
  only.pos = FALSE,
  verbose = FALSE
) %>%
  rownames_to_column("gene") %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  mutate(
    neglog10p   = safe_neglog10(p_val),
    neglog10fdr = safe_neglog10(p_val_adj),
    regulation  = case_when(
      (p_val_adj < alpha) & (avg_log2FC >  min_abs_fc_for_sig) ~ "up",
      (p_val_adj < alpha) & (avg_log2FC < -min_abs_fc_for_sig) ~ "down",
      TRUE                                                     ~ "ns"
    )
  )

message("DE regulation counts:")
print(table(de$regulation))

# ----------------------- Exports (DE table) ----------------------------------
data.table::fwrite(de, tsv_out, sep = "\t")
data.table::fwrite(de, csv_out)

# ----------------------- Diagnostics / tables --------------------------------
diag_tbl <- tibble(
  metric = c(
    "analysis_scope","assay_used","alpha","min_pct","logfc_thresh","min_abs_fc_for_sig",
    "n_total_all_combo","n_paired_matched_in_all_combo","n_paired_total_after_subset",
    "n_paired_motif_matched","n_paired_motif_unmatched",
    "n_matched_motif_used","n_bg_all_other_used",
    "rds_source","motif_xlsx_source","motif_xlsx_sheet","paired_xlsx_source",
    "paired_batch_col","paired_barcode_col","meta_key_col_for_motifs",
    "ident1_definition","ident2_definition"
  ),
  value = c(
    "PAIRED-TCR ONLY", assay_used, alpha, min_pct, logfc_thresh, min_abs_fc_for_sig,
    ncol(all_combo), length(paired_cells_in_allcombo), n_paired_total,
    n_paired_motif_matched, n_paired_motif_unmatched,
    n_matched, n_bg,
    rds_path, xlsx_path, xlsx_sheet, paired_xlsx,
    paired_batch_col, paired_barcode_col, meta_key_col,
    "PAIRED & Motif_Match==matched (NO subtype filter)",
    "PAIRED & Motif_Match!=matched (ALL remaining paired cells)"
  )
)

matched_tbl <- all_combo_de@meta.data %>%
  mutate(.cell = rownames(all_combo_de@meta.data)) %>%
  filter(DE_group_custom == "matched_motif") %>%
  select(.cell,
         any_of(c("Annotation_cluster_based","Source.x","NN_Catgory","Status","Disease","Group")),
         Motif_Match, Motif_Number, !!sym(meta_key_col))

bg_tbl <- all_combo_de@meta.data %>%
  mutate(.cell = rownames(all_combo_de@meta.data)) %>%
  filter(DE_group_custom == "bg_all_other") %>%
  select(.cell,
         any_of(c("Annotation_cluster_based","Source.x","NN_Catgory","Status","Disease","Group")),
         Motif_Match, Motif_Number, !!sym(meta_key_col))

motif_counts_matched <- matched_tbl %>%
  filter(!is.na(Motif_Number), Motif_Number != "") %>%
  count(Motif_Number, name = "n_cells") %>%
  arrange(desc(n_cells))

motif_by_annot_paired <- all_combo_paired@meta.data %>%
  mutate(Annotation_cluster_based = as.character(Annotation_cluster_based)) %>%
  count(Annotation_cluster_based, Motif_Match, name = "n_cells") %>%
  arrange(desc(n_cells))

annot_all <- meta0 %>%
  count(Annotation_cluster_based, name = "n_all_cells")

annot_paired <- all_combo_paired@meta.data %>%
  count(Annotation_cluster_based, name = "n_paired_cells")

annot_compare <- full_join(annot_all, annot_paired, by = "Annotation_cluster_based") %>%
  mutate(
    n_all_cells    = ifelse(is.na(n_all_cells), 0L, n_all_cells),
    n_paired_cells = ifelse(is.na(n_paired_cells), 0L, n_paired_cells),
    frac_paired_within_annot = ifelse(n_all_cells > 0, n_paired_cells / n_all_cells, NA_real_)
  ) %>%
  arrange(desc(n_all_cells))

# ----------------------- Excel workbook --------------------------------------
wb <- createWorkbook()
addWorksheet(wb, "DE_full");                 writeData(wb, "DE_full", de)
addWorksheet(wb, "Diagnostics");             writeData(wb, "Diagnostics", diag_tbl)
addWorksheet(wb, "MatchedMotif_cells");      writeData(wb, "MatchedMotif_cells", matched_tbl)
addWorksheet(wb, "BG_AllOther_cells");       writeData(wb, "BG_AllOther_cells", bg_tbl)
addWorksheet(wb, "Motif_counts_in_matched"); writeData(wb, "Motif_counts_in_matched", motif_counts_matched)
addWorksheet(wb, "Motif_by_annotation");     writeData(wb, "Motif_by_annotation", motif_by_annot_paired)
addWorksheet(wb, "Annot_all_vs_paired");     writeData(wb, "Annot_all_vs_paired", annot_compare)
saveWorkbook(wb, xlsx_out, overwrite = TRUE)
