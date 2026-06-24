library(dplyr)
library(readr)

# ==========================================
# 1. SETUP
# ==========================================
patient_ids <- sprintf("WGS%03d", 1:22) # 22 Patients (Dropped 23)

# Subdirectories
wgs_dir <- "WGS"
vdwgs_cbs_dir <- "vdWGS-cbs"
vdwgs_hmm_dir <- "vdWGS-hmm"
qdna_dir <- "qdnaseq"

# Baseline noise filter (keeps the dataset from ballooning with neutral background blocks)
min_log2_variance <- 0.15

# Master list to hold all data
master_list <- list()

cat("Commencing Unified CNV Extraction...\n")

# ==========================================
# 2. EXTRACTION LOOP
# ==========================================
for (p in patient_ids) {
  
  # Define expected file paths
  wgs_file <- file.path(wgs_dir, paste0(p, "_hmm-germline.cns"))
  cbs_file <- file.path(vdwgs_cbs_dir, paste0(p, "_cbs.cns"))
  hmm_file <- file.path(vdwgs_hmm_dir, paste0(p, "_hmm-germline.cns"))
  qdna_file <- file.path(qdna_dir, paste0(p, "_segments.bed"))
  
  # ---------------------------------------------------------
  # A. Helper Function for CNVkit format (.cns)
  # ---------------------------------------------------------
  load_cnvkit <- function(filepath, source_name) {
    if (file.exists(filepath)) {
      df <- read.delim(filepath, stringsAsFactors = FALSE) %>%
        # Standardize and filter
        filter(!chromosome %in% c("chrX", "chrY", "X", "Y")) %>%
        filter(abs(log2) >= min_log2_variance) %>%
        mutate(
          patient = p,
          source = source_name,
          cnv_length = end - start,
          cnv_type = ifelse(log2 > 0, "Duplication", "Deletion")
        ) %>%
        # Keep only essential columns to merge cleanly
        select(patient, source, chromosome, start, end, cnv_length, cnv_type, log2, probes, depth, weight)
      return(df)
    }
    return(NULL)
  }
  
  # ---------------------------------------------------------
  # B. Helper Function for QDNAseq format (.bed)
  # ---------------------------------------------------------
  load_qdnaseq <- function(filepath, source_name) {
    if (file.exists(filepath)) {
      # QDNAseq BED has no header: V1=chr, V2=start, V3=end, V4=type, V5=log2
      df <- read.table(filepath, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                       col.names = c("chromosome", "start", "end", "cnv_type", "log2")) %>%
        # QDNAseq often lacks the "chr" prefix, add it for consistency
        mutate(chromosome = ifelse(grepl("^chr", chromosome), chromosome, paste0("chr", chromosome))) %>%
        filter(!chromosome %in% c("chrX", "chrY", "X", "Y")) %>%
        filter(abs(log2) >= min_log2_variance) %>%
        mutate(
          patient = p,
          source = source_name,
          cnv_length = end - start,
          # QDNAseq lacks probes, depth, and weight. Fill with NA so bind_rows() works perfectly.
          probes = NA_real_,
          depth = NA_real_,
          weight = NA_real_
        ) %>%
        select(patient, source, chromosome, start, end, cnv_length, cnv_type, log2, probes, depth, weight)
      return(df)
    }
    return(NULL)
  }
  
  # ---------------------------------------------------------
  # C. Extract & Append
  # ---------------------------------------------------------
  wgs_df  <- load_cnvkit(wgs_file, "WGS-HMM")
  cbs_df  <- load_cnvkit(cbs_file, "vdWGS-CBS")
  hmm_df  <- load_cnvkit(hmm_file, "vdWGS-HMM")
  qdna_df <- load_qdnaseq(qdna_file, "vdWGS-QDNAseq")
  
  # Combine all available dataframes for this patient into a single block
  patient_combined <- bind_rows(wgs_df, cbs_df, hmm_df, qdna_df)
  
  if (nrow(patient_combined) > 0) {
    master_list[[p]] <- patient_combined
  }
}

# ==========================================
# 3. COMPILE AND EXPORT
# ==========================================
cat("Compiling Master Dataframe...\n")
unified_master <- bind_rows(master_list)

cat(sprintf("Total CNVs extracted across all pipelines: %d\n", nrow(unified_master)))
table(unified_master$source) # Quick sanity check on how many calls came from each tool

# Save to a single clean CSV
write.csv(unified_master, "Unified_Cohort_CNVs.csv", row.names = FALSE)
cat("Successfully saved to 'Unified_Cohort_CNVs.csv'\n")
