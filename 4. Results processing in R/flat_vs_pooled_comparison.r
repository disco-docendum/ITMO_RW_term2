library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(igraph)
library(tidyverse)

library(scales)

# Graph logic (with nested length calculations)

get_cluster_sizes <- function(df, reciprocal_threshold, condition_label) {
  
  process_cnv_type <- function(sub_df) {
    # if empty, return a properly formatted empty dataframe to prevent errors
    if (nrow(sub_df) == 0) return(data.frame(unique_patients = numeric(0), avg_len = numeric(0)))
    
    gr <- makeGRangesFromDataFrame(sub_df, seqnames.field = "chromosome", keep.extra.columns = TRUE)
    
    overlaps <- findOverlaps(gr, gr, ignore.strand = TRUE)
    pairs <- as.data.frame(overlaps) %>% filter(queryHits < subjectHits)
    
    if (nrow(pairs) > 0) {
      q_samples <- sub_df$sample_id[pairs$queryHits]
      s_samples <- sub_df$sample_id[pairs$subjectHits]
      pairs <- pairs[q_samples != s_samples, ]
    }
    
    if (nrow(pairs) > 0) {
      gr_q <- gr[pairs$queryHits]
      gr_s <- gr[pairs$subjectHits]
      inter <- pintersect(gr_q, gr_s, ignore.strand = TRUE)
      
      cov_q <- width(inter) / width(gr_q)
      cov_s <- width(inter) / width(gr_s)
      
      valid_edges <- pairs[cov_q >= reciprocal_threshold & cov_s >= reciprocal_threshold, ]
    } else {
      valid_edges <- data.frame(queryHits = integer(), subjectHits = integer())
    }
    
    g <- make_empty_graph(n = length(gr), directed = FALSE)
    if (nrow(valid_edges) > 0) {
      edge_vector <- c(rbind(valid_edges$queryHits, valid_edges$subjectHits))
      g <- add_edges(g, edge_vector)
    }
    
    # Extract graph membership
    cluster_df <- data.frame(
      comp_id = components(g)$membership,
      sample_id = sub_df$sample_id,
      cnv_length = sub_df$cnv_length
    )
    
    # Calculate total length of CNVs for each patient in the cluster
    patient_totals <- cluster_df %>%
      group_by(comp_id, sample_id) %>%
      summarise(tot_len = sum(cnv_length), .groups = 'drop')
    
    # Find the average of those patient totals, and count unique patients
    cluster_summary <- patient_totals %>%
      group_by(comp_id) %>%
      summarise(
        unique_patients = n_distinct(sample_id),
        avg_len = mean(tot_len),
        .groups = 'drop'
      )
    
    return(cluster_summary %>% select(unique_patients, avg_len))
  }
  
  # Process and tag
  df_dup <- process_cnv_type(df %>% filter(log2 > 0.25))
  df_del <- process_cnv_type(df %>% filter(log2 < -0.25))
  
  res_dup <- df_dup %>% mutate(
    condition = condition_label, 
    run_type = df$run_type[1], 
    cnv_category = "duplications only"
  )
  
  res_del <- df_del %>% mutate(
    condition = condition_label, 
    run_type = df$run_type[1], 
    cnv_category = "deletions only"
  )
  
  return(bind_rows(res_dup, res_del))
}

# Execute data generation
cat("calculating unique patient clusters and lengths...\n")

# Data loading function
load_wgs_run <- function(dir_name, run_label) {
  sample_ids <- sprintf("WGS%03d", 1:23)
  df_list <- list()
  
  for (sid in sample_ids) {
    file_path <- file.path(dir_name, paste0(sid, ".call.cns"))
    
    # Check if file exists to prevent errors if a sample failed
    if (file.exists(file_path)) {
      df <- read_tsv(file_path, show_col_types = FALSE) %>%
        mutate(
          sample_id = sid,
          run_type = run_label,
          cnv_length = end - start,
          cnv_type = ifelse(log2 > 0, "Duplication", "Deletion")
        ) %>%
        filter(abs(log2) > 0.25, cnv_length > 0) # Apply significance filter
      
      df_list[[sid]] <- df
    }
  }
  return(bind_rows(df_list))
}

# Load both runs
cat("Loading WGS data...\n")
df_pooled <- load_wgs_run("WGS_pooled", "Pooled Reference")
df_flat   <- load_wgs_run("WGS_flat", "Flat Reference")

# Combine for initial visualization
df_all <- bind_rows(df_pooled, df_flat)

p_1bp <- get_cluster_sizes(df_pooled, 0.00, "1bp mutual overlap")
p_25p <- get_cluster_sizes(df_pooled, 0.25, "25% mutual overlap")
p_50p <- get_cluster_sizes(df_pooled, 0.50, "50% mutual overlap")

f_1bp <- get_cluster_sizes(df_flat, 0.00, "1bp mutual overlap")
f_25p <- get_cluster_sizes(df_flat, 0.25, "25% mutual overlap")
f_50p <- get_cluster_sizes(df_flat, 0.50, "50% mutual overlap")

cluster_data <- bind_rows(p_1bp, p_25p, p_50p, f_1bp, f_25p, f_50p)

cluster_data$condition <- factor(cluster_data$condition, 
                                 levels = c("1bp mutual overlap", "25% mutual overlap", "50% mutual overlap"))

# Dynamic plotting function (with vertical median length labels)

plot_artifact_condition <- function(data, condition_name, category_name) {
  
  if (category_name == "all cnvs") {
    plot_data <- data %>% filter(condition == condition_name)
  } else {
    plot_data <- data %>% filter(condition == condition_name, cnv_category == category_name)
  }
  
  label_data <- plot_data %>%
    group_by(unique_patients, run_type) %>%
    summarise(
      n_clusters = n(),
      median_len = median(avg_len),
      .groups = 'drop'
    ) %>%
    mutate(
      label_text = ifelse(median_len >= 1000, 
                          paste0(round(median_len / 1000, 1), "k"), 
                          as.character(round(median_len, 0)))
    )
  
  p <- ggplot(plot_data, aes(x = as.factor(unique_patients), fill = run_type)) +
    geom_bar(color = "black", alpha = 0.8) +
    geom_text(data = label_data, aes(y = n_clusters, label = label_text), 
              angle = 90, hjust = -0.2, vjust = 0.5, size = 3.5, color = "black") +
    facet_wrap(~ run_type) + 
    scale_y_log10(labels = comma, expand = expansion(mult = c(0, 0.25))) + 
    scale_fill_manual(values = c("Pooled Reference" = "#2E8B57", "Flat Reference" = "#CD5C5C")) +
    theme_bw() +
    labs(
      title = paste("cnv artifacts:", condition_name, "-", category_name),
      subtitle = "x-axis: number of unique patients sharing the cnv. (max = 23)\nlabel: median of the average cnv lengths in the cluster",
      x = "number of unique patients in cluster",
      y = "number of clusters (log10 scale)"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

# Directory creation and plotting export
cat("Creating directories and exporting 9 annotated plots...\n")

conditions <- c("1bp mutual overlap", "25% mutual overlap", "50% mutual overlap")
categories <- c("all cnvs", "duplications only", "deletions only")

for (cat_name in categories) {
  dir_name <- gsub(" ", "_", cat_name)
  dir.create(dir_name, showWarnings = FALSE)
  
  for (cond_name in conditions) {
    p <- plot_artifact_condition(cluster_data, cond_name, cat_name)
    
    safe_file_name <- gsub("%", "p", cond_name)
    safe_file_name <- gsub(" ", "_", safe_file_name)
    
    file_path <- file.path(dir_name, paste0(safe_file_name, ".png"))
    
    ggsave(filename = file_path, plot = p, width = 12, height = 7.5, units = "in", dpi = 300, bg = "white")
    
    cat("saved:", file_path, "\n")
  }
}
cat("Export complete!\n")








# =======================================================================
# Pre-extraction: Combine datasets into 'all_data'
# =======================================================================

# NOTE: Adjust 'pooled_filtered' and 'flat_filtered' to exactly match 
# the names of your loaded and log2-filtered dataframes in your environment.

cat("\n--- Combining datasets into all_data for extraction ---\n")

all_data <- bind_rows(
  df_pooled %>% mutate(run_type = "Pooled Reference"),
  df_flat %>% mutate(run_type = "Flat Reference")
)

# Optional safety check to ensure required columns exist
required_cols <- c("chromosome", "start", "end", "sample_id", "log2", "run_type")
if(!all(required_cols %in% colnames(all_data))) {
  stop("ERROR: all_data is missing required columns. Please check input dataframe column names.")
}






# =======================================================================
# 5. EXTRACTION: Persistent Pooled Reference Artifacts (>= 6 Patients)
# =======================================================================

cat("\n--- Extracting high-frequency artifacts (pooled reference, >= 6 patients, 50% overlap) ---\n")

extract_artifact_coords <- function(df, target_run = "Pooled Reference", min_patients = 6) {
  
  # Filter for the specific pipeline run
  df_target <- df %>% filter(run_type == target_run)
  
  if(nrow(df_target) == 0) return(NULL)
  
  # Helper function to process duplications and deletions safely without mixing them
  process_artifacts <- function(sub_df, type_label) {
    if(nrow(sub_df) == 0) return(NULL)
    
    # Convert to GRanges
    gr <- GRanges(
      seqnames = sub_df$chromosome,
      ranges = IRanges(start = sub_df$start, end = sub_df$end),
      sample_id = sub_df$sample_id,
      log2 = sub_df$log2
    )
    
    # Find Overlaps
    overlaps <- findOverlaps(gr, gr)
    
    # Apply 50% Mutual overlap math using pintersect
    intersections <- pintersect(gr[queryHits(overlaps)], gr[subjectHits(overlaps)])
    overlap_widths <- width(intersections)
    query_widths <- width(gr[queryHits(overlaps)])
    subject_widths <- width(gr[subjectHits(overlaps)])
    
    valid_hits <- overlap_widths >= (0.5 * query_widths) & overlap_widths >= (0.5 * subject_widths)
    filtered_overlaps <- overlaps[valid_hits]
    
    # Build graph and find connected components
    edges <- data.frame(
      from = queryHits(filtered_overlaps),
      to = subjectHits(filtered_overlaps)
    ) %>% filter(from != to) # remove self-loops to keep the graph clean
    
    g <- graph_from_data_frame(edges, directed = FALSE, vertices = 1:nrow(sub_df))
    comps <- components(g)
    
    # Append cluster ID and filter by unique patient count
    sub_df$cluster_id <- comps$membership
    
    artifact_clusters <- sub_df %>%
      group_by(cluster_id) %>%
      mutate(unique_patients = n_distinct(sample_id)) %>%
      filter(unique_patients >= min_patients) %>%
      ungroup() %>%
      arrange(desc(unique_patients), cluster_id, chromosome, start)
    
    if(nrow(artifact_clusters) > 0) {
      artifact_clusters$cnv_type <- type_label
      return(artifact_clusters)
    }
    return(NULL)
  }
  
  # Process independently to maintain the biological sense
  dups <- process_artifacts(df_target %>% filter(log2 > 0.25), "Duplication")
  dels <- process_artifacts(df_target %>% filter(log2 < -0.25), "Deletion")
  
  combined <- bind_rows(dups, dels)
  return(combined)
}

# Execute the extraction
persistent_artifacts <- extract_artifact_coords(all_data, target_run = "Pooled Reference", min_patients = 3)

if(!is.null(persistent_artifacts) && nrow(persistent_artifacts) > 0) {
  
  # Summarize the clusters to remove redundancy
  cluster_summary <- persistent_artifacts %>%
    group_by(cluster_id, cnv_type, chromosome, unique_patients) %>%
    summarise(
      consensus_start = format(round(median(start)), big.mark = ",", scientific = FALSE),
      consensus_end = format(round(median(end)), big.mark = ",", scientific = FALSE),
      median_log2 = round(median(log2), 2),
      .groups = "drop"
    ) %>%
    arrange(desc(unique_patients), chromosome, consensus_start)
  
  cat(sprintf("Found %d high-frequency artifact clusters.\n\n", nrow(cluster_summary)))
  
  # Print the collapsed, clean table to the console
  print(as.data.frame(cluster_summary))
  
  # Automatically save full coordinates to a TSV for IGV pileup inspection
  out_file <- "pooled_reference_persistent_artifacts_50pct.tsv"
  write.table(persistent_artifacts, out_file, sep="\t", row.names=FALSE, quote=FALSE)
  cat(sprintf("\nSummary printed above. Full, per-patient coordinates successfully saved to '%s'\n", out_file))
  
} else {
  cat("No clusters with >= 3 unique patients.\n")
}
