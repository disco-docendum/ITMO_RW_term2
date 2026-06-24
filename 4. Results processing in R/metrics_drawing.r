# =========================================================================
# FINAL ENSEMBLE APPLICATION & VISUALIZATION PIPELINE (BIFURCATED)
# =========================================================================

library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(scales)

cat("Loading Targets, Database, and Best Parameters...\n")

# 1. Load best parameters dynamically from bifurcated grid searches
del_grid <- read.csv("Grid_Search_Optimal_Deletions.csv", stringsAsFactors = FALSE)
best_del <- del_grid %>% arrange(desc(Ens_F1)) %>% head(1)

dup_grid <- read.csv("Grid_Search_Optimal_Duplications.csv", stringsAsFactors = FALSE)
best_dup <- dup_grid %>% arrange(desc(Ens_F1)) %>% head(1)

cat(sprintf("=> Best Deletion Params: Len >= %d, Weight >= %s, F1 = %.3f\n", 
            best_del$Rel_Len, best_del$Rel_Weight, best_del$Ens_F1))
cat(sprintf("=> Best Duplication Params: Len >= %d, Weight >= %s, F1 = %.3f\n", 
            best_dup$Rel_Len, best_dup$Rel_Weight, best_dup$Ens_F1))

# 2. Load targets
targets_df <- read.table("targets_hg19_padded.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, 1:3]
colnames(targets_df) <- c("chromosome", "start", "end")
targets_df$chromosome <- ifelse(grepl("^chr", targets_df$chromosome), targets_df$chromosome, paste0("chr", targets_df$chromosome))
targets_gr <- GenomicRanges::reduce(makeGRangesFromDataFrame(targets_df))

# 3. Load and filter raw data
df_raw <- read.csv("Unified_Cohort_CNVs.csv", stringsAsFactors = FALSE) %>% 
  filter(patient != "WGS023") %>%
  filter(!chromosome %in% c("chrX", "chrY", "X", "Y", "chrM", "MT", "M"))

df_raw_gr <- makeGRangesFromDataFrame(df_raw, keep.extra.columns = TRUE)
hits <- findOverlaps(df_raw_gr, targets_gr, ignore.strand = TRUE)
df <- df_raw[unique(queryHits(hits)), ] %>% 
  mutate(
    seq_id = paste(patient, cnv_type, chromosome, sep="___"),
    Weight_Per_Probe = ifelse(probes > 0, weight / probes, 0)
  )

# 4. Fused ID Hack (universal seqinfo)
all_seq_ids <- unique(df$seq_id)
seq_info_univ <- Seqinfo(seqnames = all_seq_ids)

safe_make_gr <- function(df_sub) {
  if (nrow(df_sub) == 0) return(GRanges(seqinfo = seq_info_univ))
  
  # Structural safety check
  reserved_keywords <- c("width", "strand", "seqnames", "ranges", "element")
  df_sub <- df_sub[, !(colnames(df_sub) %in% reserved_keywords), drop = FALSE]
  
  df_sub$seq_id <- factor(df_sub$seq_id, levels = all_seq_ids)
  gr <- makeGRangesFromDataFrame(df_sub, seqnames.field = "seq_id", keep.extra.columns = TRUE)
  seqinfo(gr) <- seq_info_univ
  return(gr)
}

# =========================================================================
# 5. EXTRACT NAIVE BASELINE (raw CBS vs raw WGS)
# =========================================================================
cat("Extracting Naive Baseline Data...\n")

wgs_naive_df <- df %>% filter(source == "WGS-HMM", abs(log2) >= 0.25)
cbs_naive_df <- df %>% filter(source == "vdWGS-CBS")

wgs_naive_gr <- safe_make_gr(wgs_naive_df)
cbs_naive_gr <- safe_make_gr(cbs_naive_df)

# =========================================================================
# 6. ENSEMBLE PRE-PROCESSING (stitching and artifact masking)
# =========================================================================
cat("Applying Ensemble Pre-Processing (Stitching & Systemic Artifact Masking)...\n")

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

wgs_df  <- df %>% filter(source == "WGS-HMM", abs(log2) >= 0.25)
cbs_raw <- df %>% filter(source == "vdWGS-CBS")
hmm_raw <- stitch_hmm_fragments(df %>% filter(source == "vdWGS-HMM"))

wgs_ens_df <- mask_artifacts(wgs_df)
cbs_ens_df <- mask_artifacts(cbs_raw)
hmm_ens_df <- mask_artifacts(hmm_raw)

wgs_ens_gr <- safe_make_gr(wgs_ens_df)

# =========================================================================
# 7. APPLY BIFURCATED ENSEMBLE FILTERS
# =========================================================================
cat("Building Final Ensemble...\n")

vdwgs_max_dup_log2 <- 1.25

# Deletions
c_del_rel <- cbs_ens_df %>% filter(cnv_type == "Deletion", cnv_length >= best_del$Rel_Len, probes >= best_del$Rel_Prb, Weight_Per_Probe >= best_del$Rel_Weight, log2 <= best_del$Rel_Log2)
h_del_rel <- hmm_ens_df %>% filter(cnv_type == "Deletion", cnv_length >= best_del$Rel_Len, probes >= best_del$Rel_Prb, Weight_Per_Probe >= best_del$Rel_Weight, log2 <= best_del$Rel_Log2)
c_del_high <- cbs_ens_df %>% filter(cnv_type == "Deletion", cnv_length >= best_del$High_Len, probes >= best_del$High_Prb, Weight_Per_Probe >= best_del$High_Weight, log2 <= best_del$High_Log2)
h_del_high <- hmm_ens_df %>% filter(cnv_type == "Deletion", cnv_length >= best_del$High_Len, probes >= best_del$High_Prb, Weight_Per_Probe >= best_del$High_Weight, log2 <= best_del$High_Log2)

# Duplications
c_dup_rel <- cbs_ens_df %>% filter(cnv_type == "Duplication", cnv_length >= best_dup$Rel_Len, probes >= best_dup$Rel_Prb, Weight_Per_Probe >= best_dup$Rel_Weight, log2 >= best_dup$Rel_Log2 & log2 <= vdwgs_max_dup_log2)
h_dup_rel <- hmm_ens_df %>% filter(cnv_type == "Duplication", cnv_length >= best_dup$Rel_Len, probes >= best_dup$Rel_Prb, Weight_Per_Probe >= best_dup$Rel_Weight, log2 >= best_dup$Rel_Log2 & log2 <= vdwgs_max_dup_log2)
c_dup_high <- cbs_ens_df %>% filter(cnv_type == "Duplication", cnv_length >= best_dup$High_Len, probes >= best_dup$High_Prb, Weight_Per_Probe >= best_dup$High_Weight, log2 >= best_dup$High_Log2 & log2 <= vdwgs_max_dup_log2)
h_dup_high <- hmm_ens_df %>% filter(cnv_type == "Duplication", cnv_length >= best_dup$High_Len, probes >= best_dup$High_Prb, Weight_Per_Probe >= best_dup$High_Weight, log2 >= best_dup$High_Log2 & log2 <= vdwgs_max_dup_log2)

# Bind components
cbs_rel_gr <- safe_make_gr(bind_rows(c_del_rel, c_dup_rel))
hmm_rel_gr <- safe_make_gr(bind_rows(h_del_rel, h_dup_rel))
cbs_high_gr <- safe_make_gr(bind_rows(c_del_high, c_dup_high))
hmm_high_gr <- safe_make_gr(bind_rows(h_del_high, h_dup_high))

# Geometric logic
cons_gr <- GenomicRanges::intersect(cbs_rel_gr, hmm_rel_gr, ignore.strand = TRUE)
ensemble_gr <- GenomicRanges::union(cons_gr, cbs_high_gr)
ensemble_gr <- GenomicRanges::union(ensemble_gr, hmm_high_gr)

# =========================================================================
# 8. METRICS EVALUATION & LABELING (GROUPED BY CNV TYPE)
# =========================================================================
cat("Calculating Metrics & Generating Classifications...\n")

evaluate_segments <- function(test_gr, truth_gr, method_name) {
  test_df <- as.data.frame(test_gr) %>% mutate(patient = sapply(strsplit(as.character(seqnames), "___"), `[`, 1), cnv_type = sapply(strsplit(as.character(seqnames), "___"), `[`, 2), Method = method_name)
  truth_df <- as.data.frame(truth_gr) %>% mutate(patient = sapply(strsplit(as.character(seqnames), "___"), `[`, 1), cnv_type = sapply(strsplit(as.character(seqnames), "___"), `[`, 2), Method = method_name)
  
  ov <- findOverlaps(test_gr, truth_gr)
  is_tp <- rep(FALSE, length(test_gr))
  if(length(ov) > 0) {
    ov_int <- pintersect(test_gr[queryHits(ov)], truth_gr[subjectHits(ov)])
    ov_agg <- data.frame(qHit = queryHits(ov), ov_len = width(ov_int)) %>% group_by(qHit) %>% summarise(tot = sum(ov_len), .groups="drop")
    is_tp[ov_agg$qHit] <- (ov_agg$tot / width(test_gr)[ov_agg$qHit]) >= 0.10
  }
  
  test_df$Classification <- ifelse(is_tp, "TP", "FP")
  
  hit_truths <- unique(subjectHits(ov)[queryHits(ov) %in% which(is_tp)])
  truth_df$Classification <- "FN"
  truth_df$Classification[hit_truths] <- "Found"
  
  return(list(Test = test_df, Truth = truth_df))
}

# Evaluate naive baseline (raw CBS vs raw WGS)
eval_naive <- evaluate_segments(cbs_naive_gr, wgs_naive_gr, "Naive (Raw CBS)")

# Evaluate final ensemble (ensemble vs clean WGS)
eval_ens <- evaluate_segments(ensemble_gr, wgs_ens_gr, "Final Ensemble")

# Aggregate for per-patient whisker plots (grouped by CNV type)
aggregate_patient_metrics <- function(eval_list) {
  test_df <- eval_list$Test; truth_df <- eval_list$Truth
  
  tp_fp <- test_df %>% group_by(patient, cnv_type) %>% summarise(TP = sum(Classification=="TP"), FP = sum(Classification=="FP"), .groups="drop")
  fn <- truth_df %>% group_by(patient, cnv_type) %>% summarise(FN = sum(Classification=="FN"), .groups="drop")
  
  combined <- bind_rows(test_df %>% select(patient, cnv_type), truth_df %>% select(patient, cnv_type)) %>% distinct()
  res <- combined %>% left_join(tp_fp, by=c("patient", "cnv_type")) %>% left_join(fn, by=c("patient", "cnv_type")) %>% replace_na(list(TP=0, FP=0, FN=0))
  
  res <- res %>% mutate(
    Precision = ifelse(TP+FP==0, 0, TP/(TP+FP)),
    Recall = ifelse(TP+FN==0, 0, TP/(TP+FN)),
    F1 = ifelse(Precision+Recall==0, 0, 2*(Precision*Recall)/(Precision+Recall)),
    Method = test_df$Method[1]
  )
  return(res)
}

pat_metrics <- bind_rows(
  aggregate_patient_metrics(eval_naive), 
  aggregate_patient_metrics(eval_ens)
)

plot_metrics <- pat_metrics %>%
  group_by(Method, cnv_type) %>%
  summarise(
    Mean_Precision = mean(Precision), SD_Precision = sd(Precision),
    Mean_Recall = mean(Recall), SD_Recall = sd(Recall),
    Mean_F1 = mean(F1), SD_F1 = sd(F1),
    .groups="drop"
  ) %>%
  pivot_longer(cols = starts_with("Mean"), names_to = "Metric", values_to = "Mean") %>%
  mutate(Metric = gsub("Mean_", "", Metric)) %>%
  left_join(
    pat_metrics %>% group_by(Method, cnv_type) %>% summarise(SD_Precision=sd(Precision), SD_Recall=sd(Recall), SD_F1=sd(F1), .groups="drop") %>%
      pivot_longer(cols = starts_with("SD"), names_to = "Metric", values_to = "SD") %>% mutate(Metric = gsub("SD_", "", Metric)),
    by = c("Method", "cnv_type", "Metric")
  )

plot_metrics$Method <- factor(plot_metrics$Method, levels = c("Naive (Raw CBS)", "Final Ensemble"))
plot_metrics$Metric <- factor(plot_metrics$Metric, levels = c("Precision", "Recall", "F1"))

# =========================================================================
# 9. VISUALIZATIONS
# =========================================================================
cat("Generating Plots...\n")

# GRAPH 1: Naive vs final ensemble bar chart (faceted by deletion/duplication)
p1 <- ggplot(plot_metrics, aes(x = Metric, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color="black", alpha=0.85) +
  geom_errorbar(aes(ymin = pmax(0, Mean-SD), ymax = Mean+SD), position = position_dodge(width = 0.8), width = 0.25) +
  scale_fill_manual(values = c("Naive (Raw CBS)" = "#E41A1C", "Final Ensemble" = "#377EB8")) +
  facet_wrap(~ cnv_type) +  
  theme_bw(base_size = 14) +
  labs(title = "Performance: naive pipeline vs. final ensemble",
       subtitle = "Metrics separated by deletions and duplications (error bars = std. dev)",
       x = "Metric", y = "Score") +
  theme(legend.position = "bottom", 
        plot.title = element_text(face="bold", hjust=0.5), 
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_rect(fill="grey90"),
        strip.text = element_text(face="bold"))

ggsave("Performance_Naive_vs_Ensemble.png", plot = p1, width = 10, height = 6, dpi = 300)

# GRAPH 2: length distributions with TRUE NEGATIVES
# Build combined dataframe using the final ensemble evaluation
ens_test <- eval_ens$Test %>% select(cnv_length = width, cnv_type, Classification)
ens_truth_fn <- eval_ens$Truth %>% filter(Classification == "FN") %>% select(cnv_length = width, cnv_type, Classification)
gt_full <- eval_ens$Truth %>% select(cnv_length = width, cnv_type) %>% mutate(Classification = "Ground Truth")

# Calculate true negatives (raw artifacts that were successfully filtered out)
cbs_all_gr <- safe_make_gr(cbs_ens_df)
hmm_all_gr <- safe_make_gr(hmm_ens_df)
raw_pool_gr <- GenomicRanges::union(cbs_all_gr, hmm_all_gr)

raw_eval <- evaluate_segments(raw_pool_gr, wgs_ens_gr, "Raw Pool")
raw_fp_gr <- raw_pool_gr[raw_eval$Test$Classification == "FP"]
raw_fp_df <- raw_eval$Test %>% filter(Classification == "FP")

ov_ens <- findOverlaps(raw_fp_gr, ensemble_gr)
is_retained <- rep(FALSE, length(raw_fp_gr))
if(length(ov_ens) > 0) {
  ov_int <- pintersect(raw_fp_gr[queryHits(ov_ens)], ensemble_gr[subjectHits(ov_ens)])
  ov_agg <- data.frame(qHit = queryHits(ov_ens), ov_len = width(ov_int)) %>% group_by(qHit) %>% summarise(tot = sum(ov_len), .groups="drop")
  is_retained[ov_agg$qHit] <- (ov_agg$tot / width(raw_fp_gr)[ov_agg$qHit]) >= 0.10
}

ens_truth_tn <- raw_fp_df[!is_retained, ] %>% mutate(Classification = "TN") %>% select(cnv_length = width, cnv_type, Classification)

# Bind all 5 classifications together
length_df <- bind_rows(ens_test, ens_truth_fn, gt_full, ens_truth_tn)

length_df <- length_df %>%
  mutate(Classification = case_when(
    Classification == "Ground Truth" ~ "1. WGS Ground Truth",
    Classification == "TP" ~ "2. vdWGS ensemble true positives",
    Classification == "FP" ~ "3. vdWGS ensemble false positives",
    Classification == "FN" ~ "4. Filtered-out vdWGS true positives\n(false negatives of filtration)",
    Classification == "TN" ~ "5. Correctly filtered vdWGS artifacts\n(true negatives)"
  ))

length_df$Classification <- factor(length_df$Classification, levels = c(
  "1. WGS Ground Truth",
  "2. vdWGS ensemble true positives",
  "3. vdWGS ensemble false positives",
  "4. Filtered-out vdWGS true positives\n(false negatives of filtration)",
  "5. Correctly filtered vdWGS artifacts\n(true negatives)"
))

p2 <- ggplot(length_df, aes(x = cnv_length, fill = cnv_type)) +
  geom_histogram(bins = 40, color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("Deletion" = "#2C7BB6", "Duplication" = "#D7191C")) +
  facet_grid(Classification ~ cnv_type) +  
  scale_x_log10(labels = comma) +
  theme_bw(base_size = 14) +
  labs(
    title = "Final ensemble CNV length evaluation",
    subtitle = "Geometric consensus vs clean WGS ground truth",
    x = "CNV length (bp) [Log10 Scale]",
    y = "Absolute count (frequency)"
  ) +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0), 
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

ggsave("Ensemble_Length_Distributions.png", plot = p2, width = 14, height = 12, dpi = 300)

# =========================================================================
# 10. EXPORT FINAL BED FILE
# =========================================================================
cat("Exporting Final Ensemble BED File...\n")

final_bed <- as.data.frame(ensemble_gr) %>%
  mutate(
    chromosome = sapply(strsplit(as.character(seqnames), "___"), `[`, 3),
    patient = sapply(strsplit(as.character(seqnames), "___"), `[`, 1),
    cnv_type = sapply(strsplit(as.character(seqnames), "___"), `[`, 2)
  ) %>%
  select(chromosome, start, end, patient, cnv_type, width) %>%
  rename(cnv_length = width) %>%
  arrange(patient, chromosome, start)

write.table(final_bed, "Final_Ensemble_CNVs.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

cat("[SUCCESS] Pipeline Complete. Plots and BED file generated.\n")
