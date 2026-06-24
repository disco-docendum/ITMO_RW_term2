# =========================================================================
# BIFURCATED GRID SEARCH: INDEPENDENT OPTIMIZATION FOR DEL AND DUP
# =========================================================================

library(dplyr)
library(GenomicRanges)

cat("Loading Target Regions and Filtering Database...\n")

# Load targeted regions BED file
targets_df <- read.table("targets_hg19_padded.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, 1:3]
colnames(targets_df) <- c("chromosome", "start", "end")
targets_df$chromosome <- ifelse(grepl("^chr", targets_df$chromosome), targets_df$chromosome, paste0("chr", targets_df$chromosome))
targets_gr <- GenomicRanges::reduce(makeGRangesFromDataFrame(targets_df))

# Load the unified CNV database
df_raw <- read.csv("Unified_Cohort_CNVs.csv", stringsAsFactors = FALSE) %>% 
  filter(patient != "WGS023") %>%
  filter(!chromosome %in% c("chrX", "chrY", "X", "Y", "chrM", "MT", "M"))

# Retain ONLY segments overlapping targeted exons
df_raw_gr <- makeGRangesFromDataFrame(df_raw, keep.extra.columns = TRUE)
hits <- findOverlaps(df_raw_gr, targets_gr, ignore.strand = TRUE)
df <- df_raw[unique(queryHits(hits)), ]

# Biological constraints
wgs_min_log2 <- 0.25
vdwgs_max_dup_log2 <- 1.25
max_cohort_freq <- 14 

# ID Hack & Weight Generation
df <- df %>% 
  mutate(
    seq_id = paste(patient, cnv_type, chromosome, sep="___"),
    Weight_Per_Probe = ifelse(probes > 0, weight / probes, 0)
  )

# FUSED ID HACK: UNIFIED SEQINFO
all_seq_ids <- unique(df$seq_id)
seq_info_univ <- Seqinfo(seqnames = all_seq_ids)

safe_make_gr <- function(df_sub) {
  if (nrow(df_sub) == 0) return(GRanges(seqinfo = seq_info_univ))
  reserved_keywords <- c("width", "strand", "seqnames", "ranges", "element")
  df_sub <- df_sub[, !(colnames(df_sub) %in% reserved_keywords), drop = FALSE]
  df_sub$seq_id <- factor(df_sub$seq_id, levels = all_seq_ids)
  gr <- makeGRangesFromDataFrame(df_sub, seqnames.field = "seq_id", keep.extra.columns = TRUE)
  seqinfo(gr) <- seq_info_univ
  return(gr)
}

# =========================================================================
# PRE-PROCESSING (stitching and artifact masking)
# =========================================================================
cat("Applying Pre-Processing Algorithms...\n")

stitch_hmm_fragments <- function(df_hmm, max_gap = 50000) {
  if(nrow(df_hmm) == 0) return(df_hmm)
  df_hmm <- df_hmm %>% filter(probes >= 3, (weight / probes) >= 0.50)
  gr <- safe_make_gr(df_hmm)
  merged_gr <- GenomicRanges::reduce(gr, min.gapwidth = max_gap)
  ov <- findOverlaps(merged_gr, gr)
  
  agg_df <- data.frame(m_idx = queryHits(ov), probes = df_hmm$probes[subjectHits(ov)], log2 = df_hmm$log2[subjectHits(ov)], len = width(gr)[subjectHits(ov)], weight = df_hmm$weight[subjectHits(ov)]) %>%
    group_by(m_idx) %>% summarise(total_probes = sum(probes), weighted_log2 = sum(log2 * len)/sum(len), total_weight = sum(weight), .groups = "drop")
  
  merged_df <- as.data.frame(merged_gr, row.names = NULL) %>%
    mutate(seq_id = as.character(seqnames), cnv_length = width, probes = agg_df$total_probes, log2 = agg_df$weighted_log2, weight = agg_df$total_weight, source = "vdWGS-HMM", Weight_Per_Probe = ifelse(probes>0, weight/probes, 0))
  
  parts <- strsplit(merged_df$seq_id, "___")
  merged_df$patient <- sapply(parts, `[`, 1); merged_df$cnv_type <- sapply(parts, `[`, 2); merged_df$chromosome <- sapply(parts, `[`, 3)
  return(merged_df)
}

mask_artifacts <- function(df_in, max_allowed=14) {
  if(nrow(df_in) == 0) return(df_in)
  reserved_keywords <- c("width", "strand", "seqnames", "ranges", "element")
  df_clean_cols <- df_in[, !(colnames(df_in) %in% reserved_keywords), drop = FALSE]
  
  gr <- makeGRangesFromDataFrame(df_clean_cols, seqnames.field = "chromosome", keep.extra.columns = TRUE)
  ov <- findOverlaps(gr, gr, ignore.strand = TRUE)
  patient_counts <- data.frame(qHit = queryHits(ov), patient = df_in$patient[subjectHits(ov)]) %>% group_by(qHit) %>% summarise(unique_patients = n_distinct(patient), .groups = "drop")
  df_in$cohort_frequency <- patient_counts$unique_patients
  return(df_in %>% filter(cohort_frequency <= max_allowed))
}

wgs_df  <- df %>% filter(source == "WGS-HMM", abs(log2) >= wgs_min_log2)
cbs_raw <- df %>% filter(source == "vdWGS-CBS")
hmm_raw <- stitch_hmm_fragments(df %>% filter(source == "vdWGS-HMM"))

# Clean all datasets including truth
wgs_df  <- mask_artifacts(wgs_df)
cbs_raw <- mask_artifacts(cbs_raw)
hmm_raw <- mask_artifacts(hmm_raw)

wgs_gr <- safe_make_gr(wgs_df)

# =========================================================================
# DECOUPLED PARAMETER GRIDS
# =========================================================================

del_grid <- expand.grid(
  rel_len = c(15000, 39000),          
  rel_prb = c(4, 5),                  
  rel_log2 = c(-0.25, -0.30),        
  rel_weight = c(0.85, 0.95),      
  high_len = c(39000),
  high_prb = c(6, 8),             
  high_log2 = c(-0.40, -0.50),       
  high_weight = c(1.00, 1.15)            
)

# Duplication grid: pushed to higher weights to combat the massive FP spike
dup_grid <- expand.grid(
  rel_len = c(15000, 39000, 50000),          
  rel_prb = c(4, 5, 6),                  
  rel_log2 = c(0.20, 0.25, 0.30),        
  rel_weight = c(0.95, 1.05, 1.15),      
  high_len = c(39000, 50000),
  high_prb = c(6, 8),             
  high_log2 = c(0.35, 0.45),       
  high_weight = c(1.10, 1.25)            
)

# =========================================================================
# ISOLATED EVALUATION ENGINE
# =========================================================================

eval_metrics <- function(query_gr, truth_gr) {
  if(length(query_gr) == 0) return(c(C_TP=0, C_FP=0, T_Found=0, T_Miss=length(truth_gr)))
  
  ov <- findOverlaps(query_gr, truth_gr)
  if(length(ov) == 0) return(c(C_TP=0, C_FP=length(query_gr), T_Found=0, T_Miss=length(truth_gr)))
  
  ov_int <- pintersect(query_gr[queryHits(ov)], truth_gr[subjectHits(ov)])
  ov_df <- data.frame(qHit = queryHits(ov), ov_len = width(ov_int)) %>% group_by(qHit) %>% summarise(total_ov = sum(ov_len), .groups="drop")
  
  is_tp <- rep(FALSE, length(query_gr))
  is_tp[ov_df$qHit] <- (ov_df$total_ov / width(query_gr)[ov_df$qHit]) >= 0.10
  
  C_TP <- sum(is_tp)
  C_FP <- length(query_gr) - C_TP
  
  hit_truths <- unique(subjectHits(ov)[queryHits(ov) %in% which(is_tp)])
  T_Found <- length(hit_truths)
  T_Miss <- length(truth_gr) - T_Found
  
  return(c(C_TP=C_TP, C_FP=C_FP, T_Found=T_Found, T_Miss=T_Miss))
}

calc_f1 <- function(res) {
  tp <- res["C_TP"]; fp <- res["C_FP"]; t_found <- res["T_Found"]; t_miss <- res["T_Miss"]
  prec <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
  rec  <- ifelse((t_found + t_miss) == 0, 0, t_found / (t_found + t_miss))
  f1 <- ifelse((prec + rec) == 0, 0, 2 * ((prec * rec) / (prec + rec)))
  return(f1)
}

run_grid_search <- function(grid_params, target_type) {
  cat(sprintf("\nRunning Isolated Optimization for: %s (Iterations: %d)\n", target_type, nrow(grid_params)))
  
  cbs_sub <- cbs_raw %>% filter(cnv_type == target_type)
  hmm_sub <- hmm_raw %>% filter(cnv_type == target_type)
  wgs_sub <- wgs_gr[wgs_gr$cnv_type == target_type]
  
  results_list <- list()
  
  for(i in 1:nrow(grid_params)) {
    p <- grid_params[i, ]
    if (i %% 50 == 0) cat(sprintf("  Processed %d / %d...\n", i, nrow(grid_params)))
    
    if (target_type == "Deletion") {
      c_rel <- cbs_sub %>% filter(cnv_length >= p$rel_len, probes >= p$rel_prb, Weight_Per_Probe >= p$rel_weight, log2 <= p$rel_log2)
      h_rel <- hmm_sub %>% filter(cnv_length >= p$rel_len, probes >= p$rel_prb, Weight_Per_Probe >= p$rel_weight, log2 <= p$rel_log2)
      c_high <- cbs_sub %>% filter(cnv_length >= p$high_len, probes >= p$high_prb, Weight_Per_Probe >= p$high_weight, log2 <= p$high_log2)
      h_high <- hmm_sub %>% filter(cnv_length >= p$high_len, probes >= p$high_prb, Weight_Per_Probe >= p$high_weight, log2 <= p$high_log2)
    } else {
      c_rel <- cbs_sub %>% filter(cnv_length >= p$rel_len, probes >= p$rel_prb, Weight_Per_Probe >= p$rel_weight, log2 >= p$rel_log2 & log2 <= vdwgs_max_dup_log2)
      h_rel <- hmm_sub %>% filter(cnv_length >= p$rel_len, probes >= p$rel_prb, Weight_Per_Probe >= p$rel_weight, log2 >= p$rel_log2 & log2 <= vdwgs_max_dup_log2)
      c_high <- cbs_sub %>% filter(cnv_length >= p$high_len, probes >= p$high_prb, Weight_Per_Probe >= p$high_weight, log2 >= p$high_log2 & log2 <= vdwgs_max_dup_log2)
      h_high <- hmm_sub %>% filter(cnv_length >= p$high_len, probes >= p$high_prb, Weight_Per_Probe >= p$high_weight, log2 >= p$high_log2 & log2 <= vdwgs_max_dup_log2)
    }
    
    cbs_rel_gr <- safe_make_gr(c_rel); hmm_rel_gr <- safe_make_gr(h_rel)
    cbs_high_gr <- safe_make_gr(c_high); hmm_high_gr <- safe_make_gr(h_high)
    
    cons_gr <- GenomicRanges::intersect(cbs_rel_gr, hmm_rel_gr, ignore.strand = TRUE)
    ensemble_gr <- GenomicRanges::union(cons_gr, cbs_high_gr)
    ensemble_gr <- GenomicRanges::union(ensemble_gr, hmm_high_gr)
    
    res_ens <- eval_metrics(ensemble_gr, wgs_sub)
    
    results_list[[i]] <- data.frame(
      Type = target_type,
      Rel_Len = p$rel_len, Rel_Prb = p$rel_prb, Rel_Log2 = p$rel_log2, Rel_Weight = p$rel_weight,
      High_Len = p$high_len, High_Prb = p$high_prb, High_Log2 = p$high_log2, High_Weight = p$high_weight,
      Ens_TP = res_ens["C_TP"], Ens_FP = res_ens["C_FP"], Ens_FN = res_ens["T_Miss"],
      Ens_F1 = calc_f1(res_ens)
    )
  }
  return(bind_rows(results_list))
}

# =========================================================================
# EXECUTE AND EXTRACT
# =========================================================================

del_results <- run_grid_search(del_grid, "Deletion")
dup_results <- run_grid_search(dup_grid, "Duplication")

write.csv(del_results, "Grid_Search_Optimal_Deletions.csv", row.names = FALSE)
write.csv(dup_results, "Grid_Search_Optimal_Duplications.csv", row.names = FALSE)

best_del <- del_results %>% arrange(desc(Ens_F1)) %>% head(1)
best_dup <- dup_results %>% arrange(desc(Ens_F1)) %>% head(1)

cat("\n========================================================================\n")
cat(" BEST PARAMETERS FOUND\n")
cat("========================================================================\n")
cat(sprintf("DELETIONS   | F1: %.3f | Rel Weight: %s | High Weight: %s\n", best_del$Ens_F1, best_del$Rel_Weight, best_del$High_Weight))
cat(sprintf("DUPLICATIONS| F1: %.3f | Rel Weight: %s | High Weight: %s\n", best_dup$Ens_F1, best_dup$Rel_Weight, best_dup$High_Weight))

# =========================================================================
# ASSEMBLE FINAL OPTIMIZED HYBRID ENSEMBLE
# =========================================================================
cat("\nAssembling Final Labeled Data Matrix based on bifurcated parameters...\n")

# Filter deletions with best del params
c_del_rel <- cbs_raw %>% filter(cnv_type == "Deletion", cnv_length >= best_del$Rel_Len, probes >= best_del$Rel_Prb, Weight_Per_Probe >= best_del$Rel_Weight, log2 <= best_del$Rel_Log2)
h_del_rel <- hmm_raw %>% filter(cnv_type == "Deletion", cnv_length >= best_del$Rel_Len, probes >= best_del$Rel_Prb, Weight_Per_Probe >= best_del$Rel_Weight, log2 <= best_del$Rel_Log2)
c_del_high <- cbs_raw %>% filter(cnv_type == "Deletion", cnv_length >= best_del$High_Len, probes >= best_del$High_Prb, Weight_Per_Probe >= best_del$High_Weight, log2 <= best_del$High_Log2)
h_del_high <- hmm_raw %>% filter(cnv_type == "Deletion", cnv_length >= best_del$High_Len, probes >= best_del$High_Prb, Weight_Per_Probe >= best_del$High_Weight, log2 <= best_del$High_Log2)

# Filter duplications with best dup params
c_dup_rel <- cbs_raw %>% filter(cnv_type == "Duplication", cnv_length >= best_dup$Rel_Len, probes >= best_dup$Rel_Prb, Weight_Per_Probe >= best_dup$Rel_Weight, log2 >= best_dup$Rel_Log2 & log2 <= vdwgs_max_dup_log2)
h_dup_rel <- hmm_raw %>% filter(cnv_type == "Duplication", cnv_length >= best_dup$Rel_Len, probes >= best_dup$Rel_Prb, Weight_Per_Probe >= best_dup$Rel_Weight, log2 >= best_dup$Rel_Log2 & log2 <= vdwgs_max_dup_log2)
c_dup_high <- cbs_raw %>% filter(cnv_type == "Duplication", cnv_length >= best_dup$High_Len, probes >= best_dup$High_Prb, Weight_Per_Probe >= best_dup$High_Weight, log2 >= best_dup$High_Log2 & log2 <= vdwgs_max_dup_log2)
h_dup_high <- hmm_raw %>% filter(cnv_type == "Duplication", cnv_length >= best_dup$High_Len, probes >= best_dup$High_Prb, Weight_Per_Probe >= best_dup$High_Weight, log2 >= best_dup$High_Log2 & log2 <= vdwgs_max_dup_log2)

# Unify back into the geometric engine
cbs_rel_gr <- safe_make_gr(bind_rows(c_del_rel, c_dup_rel))
hmm_rel_gr <- safe_make_gr(bind_rows(h_del_rel, h_dup_rel))
cbs_high_gr <- safe_make_gr(bind_rows(c_del_high, c_dup_high))
hmm_high_gr <- safe_make_gr(bind_rows(h_del_high, h_dup_high))

cons_gr <- GenomicRanges::intersect(cbs_rel_gr, hmm_rel_gr, ignore.strand = TRUE)
ensemble_gr <- GenomicRanges::union(cons_gr, cbs_high_gr)
ensemble_gr <- GenomicRanges::union(ensemble_gr, hmm_high_gr)

# Final labeling for export
raw_ensemble_df <- as.data.frame(ensemble_gr, row.names = NULL) %>% mutate(seq_id = as.character(seqnames))
ov <- findOverlaps(ensemble_gr, wgs_gr)
is_tp <- rep(FALSE, length(ensemble_gr))

if(length(ov) > 0) {
  ov_int <- pintersect(ensemble_gr[queryHits(ov)], wgs_gr[subjectHits(ov)])
  ov_df <- data.frame(qHit = queryHits(ov), ov_len = width(ov_int)) %>% group_by(qHit) %>% summarise(total_ov = sum(ov_len))
  is_tp[ov_df$qHit] <- (ov_df$total_ov / width(ensemble_gr)[ov_df$qHit]) >= 0.10
}

raw_ensemble_df$Classification <- ifelse(is_tp, "True Positive", "False Positive")
raw_ensemble_df$patient <- sapply(strsplit(raw_ensemble_df$seq_id, "___"), `[`, 1)
raw_ensemble_df$cnv_type <- sapply(strsplit(raw_ensemble_df$seq_id, "___"), `[`, 2)
raw_ensemble_df$chromosome <- sapply(strsplit(raw_ensemble_df$seq_id, "___"), `[`, 3)

final_export <- raw_ensemble_df %>% select(patient, chromosome, start, end, width, cnv_type, Classification) %>% rename(cnv_length = width)
write.csv(final_export, "Labeled_vdWGS_Calls_For_Analysis.csv", row.names = FALSE)

cat("[INFO] Bifurcated Analysis Complete. Matrices and Labeled data exported.\n")
