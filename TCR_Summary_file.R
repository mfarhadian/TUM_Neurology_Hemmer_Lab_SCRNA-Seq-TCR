############################################################
# TCR summary by CDR3alpha_CDR3beta + VDJ (both chains)
# + CSF/PBMC counts + unique patient counts
# + NEW: Patient.IDs_TCR.study (unique MS/CTR IDs per combo)

############################################################

library(data.table)
library(openxlsx)

############################################################
# 1) Load data
############################################################

DT <- read.xlsx("~/TCR_Paired_chains.xlsx")

# Key: CDR3alpha_CDR3beta
DT[, Merge_file := paste0(Contig_cdr3_TRA_1, "_", Contig_cdr3_TRB_1)]

# Clean Source
DT[, Source := toupper(trimws(Source))]

############################################################
# 2) Counts by Source (CSF/PBMC)
############################################################
DT_src <- DT[Source %in% c("CSF", "PBMC")]

source_counts <- DT_src[, .N, by = .(Merge_file, Source)]
source_counts <- dcast(source_counts, Merge_file ~ Source, value.var = "N", fill = 0)

if (!"CSF"  %in% names(source_counts)) source_counts[, CSF := 0L]
if (!"PBMC" %in% names(source_counts)) source_counts[, PBMC := 0L]

setnames(source_counts,
         old = c("CSF", "PBMC"),
         new = c("CSF_Count_Number", "PBMC_Count_Number"))

############################################################
# 3) Unique patient counts by Source (CSF/PBMC)
############################################################
patient_counts_by_source <- DT_src[, .(Unique_Patient_Count = uniqueN(NN_IshId_HLA)),
                                   by = .(Merge_file, Source)]
patient_counts_by_source <- dcast(patient_counts_by_source,
                                  Merge_file ~ Source,
                                  value.var = "Unique_Patient_Count",
                                  fill = 0)

if (!"CSF"  %in% names(patient_counts_by_source)) patient_counts_by_source[, CSF := 0L]
if (!"PBMC" %in% names(patient_counts_by_source)) patient_counts_by_source[, PBMC := 0L]

setnames(patient_counts_by_source,
         old = c("CSF", "PBMC"),
         new = c("Number of merged in Patient with CSF",
                 "Number of merged in Patient with PBMC"))

############################################################
# 4) Unique patients per diagnosis group (MS/CTRL)
############################################################
tmp_diag <- unique(DT[, .(Merge_file, NN_Diagnosis.group, NN_IshId_HLA)])
diag_counts <- tmp_diag[, .(Patient_Count = .N), by = .(Merge_file, NN_Diagnosis.group)]

phenotype_patient_counts <- dcast(diag_counts,
                                  Merge_file ~ NN_Diagnosis.group,
                                  value.var = "Patient_Count",
                                  fill = 0)

############################################################
# 4b) Unique MS/CTR
############################################################

patient_ids_study_by_combo <- DT[
  !is.na(Patient.IDs_TCR.study) & Patient.IDs_TCR.study != "",
  .(Patient.IDs_TCR.study = paste(sort(unique(Patient.IDs_TCR.study)), collapse = ",")),
  by = Merge_file
]

############################################################
# 5) MAIN SUMMARY: CDR3 + VDJ for BOTH CHAINS
############################################################
TCR_summary <- DT[, .(
  Count        = .N,
  Num_Patients = uniqueN(NN_IshId_HLA),
  
  # CDR3
  CDR3_ALPHA = Contig_cdr3_TRA_1[1],
  CDR3_BETA  = Contig_cdr3_TRB_1[1],
  
  # Contig gene calls
  TRAV = paste(unique(Contig_v_gene_TRA_1), collapse = ","),
  TRAJ = paste(unique(Contig_j_gene_TRA_1), collapse = ","),
  TRBV = paste(unique(Contig_v_gene_TRB_1), collapse = ","),
  TRBJ = paste(unique(Contig_j_gene_TRB_1), collapse = ","),
  
  # AIRR gene calls
  Airr_TRAV = paste(unique(Airr_v_call_TRA_1), collapse = ","),
  Airr_TRAJ = paste(unique(Airr_j_call_TRA_1), collapse = ","),
  Airr_TRBV = paste(unique(Airr_v_call_TRB_1), collapse = ","),
  Airr_TRBJ = paste(unique(Airr_j_call_TRB_1), collapse = ","),
  
  # Metadata lists
  NN_Diagnosis.group = paste(unique(NN_Diagnosis.group), collapse = ","),
  Airr_Sample        = paste(unique(Airr_Sample), collapse = ","),
  Source             = paste(unique(Source), collapse = ","),
  NN_IshId_HLA        = paste(unique(NN_IshId_HLA), collapse = ",")
), by = Merge_file]

############################################################
# 6) JOIN 
############################################################
TCR_summary_final <- merge(TCR_summary, phenotype_patient_counts, by = "Merge_file", all.x = TRUE)
TCR_summary_final <- merge(TCR_summary_final, source_counts, by = "Merge_file", all.x = TRUE)
TCR_summary_final <- merge(TCR_summary_final, patient_counts_by_source, by = "Merge_file", all.x = TRUE)


# add Patient.IDs_TCR.study summary per combo
TCR_summary_final <- merge(TCR_summary_final, patient_ids_study_by_combo, by = "Merge_file", all.x = TRUE)



