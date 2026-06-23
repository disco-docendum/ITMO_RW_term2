#!/usr/bin/env Rscript

# ==============================================================================
# QDNAseq Pipeline (Optimized for low-coverage WGS 0.5x-1x)
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) stop("Usage: Rscript qdnaseq_wgs_pipeline.R <out_dir> <bin_size> <exclude_bed> <chromosome> <threads> <alpha> <bam1> ...")

output_dir  <- args[1]
bin_size_kb <- as.numeric(args[2])
exclude_bed <- args[3]
chromosome  <- args[4]
threads     <- as.numeric(args[5])
alpha_val   <- as.numeric(args[6])
tumor_bams  <- args[7:length(args)]

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(QDNAseq)
  library(QDNAseq.hg19) 
  library(Biobase)
  library(GenomicRanges)
  library(future)
})

# --- DYNAMIC MULTICORE CONFIGURATION ---
if (!is.na(threads) && threads > 1) {
  future::plan("multisession", workers = threads)
  message(paste0("[INFO] Multicore processing enabled: ", threads, " workers"))
} else {
  future::plan("sequential")
  message("[INFO] Sequential processing enabled (1 worker)")
}

message(paste0("[INFO] Mode: WGS Pipeline (Bin Size: ", bin_size_kb, "kb)"))

# --- HELPER FUNCTIONS ---

exportManualResults <- function(qdna_object, sample_name, output_dir) {
  
  all_segs    <- assayDataElement(qdna_object, "segmented")
  sample_segs <- all_segs[, 1]
  chroms      <- fData(qdna_object)$chromosome
  starts      <- fData(qdna_object)$start
  ends        <- fData(qdna_object)$end
  
  segment_list <- list()
  
  for (chr in unique(chroms)) {
    idx <- which(chroms == chr)
    if (length(idx) == 0) next
    
    chr_segs   <- sample_segs[idx]
    chr_starts <- starts[idx]
    chr_ends   <- ends[idx]
    
    # FIX: Find valid (non-NA) bins to prevent CBS gap bridging
    valid_idx <- which(!is.na(chr_segs))
    if (length(valid_idx) == 0) next
    
    # Group contiguous valid bins to naturally jump over unmappable/blacklisted regions
    breaks <- c(0, which(diff(valid_idx) > 1), length(valid_idx))
    
    for (b in 1:(length(breaks) - 1)) {
      block_valid <- valid_idx[(breaks[b] + 1):breaks[b + 1]]
      
      block_segs   <- chr_segs[block_valid]
      block_starts <- chr_starts[block_valid]
      block_ends   <- chr_ends[block_valid]
      
      rle_res <- rle(block_segs)
      end_idx   <- cumsum(rle_res$lengths)
      start_idx <- c(1, head(end_idx, -1) + 1)
      
      # Log2 transform the linear segmented values for thresholding
      log2_means <- log2(rle_res$values + 1e-6) 
      
      df <- data.frame(
        ID    = sample_name,
        chrom = chr,
        start = block_starts[start_idx],
        end   = block_ends[end_idx],
        num_bins = rle_res$lengths,
        seg_mean_log2 = log2_means,
        stringsAsFactors = FALSE
      )
      
      segment_list[[length(segment_list) + 1]] <- df
    }
  }
  
  if (length(segment_list) == 0) return(NULL)
  final_df <- do.call(rbind, segment_list)
  
  write.table(final_df, file = file.path(output_dir, paste0(sample_name, ".seg")), 
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  # FIX: static thresholding for duplications and deletions
  final_df$call <- "Normal"
  final_df$call[final_df$seg_mean_log2 > 0.25]  <- "Duplication"            
  final_df$call[final_df$seg_mean_log2 < -0.25] <- "Deletion"            
  
  bed_df <- final_df[final_df$call != "Normal", ]
  
  if (nrow(bed_df) > 0) {
    write.table(bed_df[, c("chrom", "start", "end", "call", "seg_mean_log2")], 
                file = file.path(output_dir, paste0(sample_name, "_segments.bed")), 
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    message(paste0("[INFO] Detected ", nrow(bed_df), " CNV segments. (Static Threshold: +/- 0.25 log2)"))
  } else {
    message("[WARN] No CNVs detected.")
  }
}

# --- PREPARE ANNOTATIONS ---
message("[INFO] Loading annotations (keeping native QDNAseq filters)...")
bins <- getBinAnnotations(binSize = bin_size_kb, genome = "hg19")

# --- CHROMOSOME FILTERING LOGIC ---
if (tolower(chromosome) != "all") {
  clean_chr <- gsub("^chr", "", chromosome) 
  message(paste0("[INFO] Filtering annotations for chromosome: ", clean_chr))
  # FIX: Use fData() instead of pData() to filter features
  bins <- bins[fData(bins)$chromosome == clean_chr, ]
} else {
  message("[INFO] Running analysis on ALL chromosomes.")
}

# --- APPLY EXCLUSION BED (ROI Masking) ---
message(paste0("[INFO] Processing ROI exclusion mask: ", exclude_bed))

ex_df <- tryCatch({
  read.table(exclude_bed, header=FALSE, sep="\t", stringsAsFactors=FALSE)
}, error=function(e) { stop("Could not read exclusion BED file.") })

colnames(ex_df)[1:3] <- c("chrom", "start", "end")
ex_df$chrom <- gsub("^chr", "", ex_df$chrom)

ex_gr <- makeGRangesFromDataFrame(ex_df, keep.extra.columns=FALSE)
# FIX: Use fData() to correctly extract bin coordinates for GRanges
bin_gr <- makeGRangesFromDataFrame(fData(bins), keep.extra.columns=FALSE)

hits <- findOverlaps(bin_gr, ex_gr)

if (length(hits) > 0) {
  masked_indices <- unique(queryHits(hits))
  # FIX: write the mask status strictly to the feature data (fData)
  fData(bins)$use[masked_indices] <- FALSE
  message(paste0("[INFO] Masked ", length(masked_indices), " bins overlapping with high-coverage ROI regions."))
} else {
  message("[INFO] No bins overlapped with exclusion regions.")
}

# --- COUNTING & CORRECTION ---
message("[STEP 1/3] Counting Reads...")
counts <- binReadCounts(bins, bamfiles = tumor_bams)

message("[STEP 2/3] Applying GC and Mappability Corrections...")
counts <- applyFilters(counts, residual=TRUE, blacklist=TRUE)
counts <- estimateCorrection(counts)
cn <- correctBins(counts)

message("[INFO] Normalizing & Imputing...")
cn <- normalizeBins(cn)
cn <- smoothOutlierBins(cn)

# --- SEGMENTATION ---
message("[STEP 3/3] Segmenting (Parameters reverted to biological reality)...")
# FIX: relaxed alpha and disabled destructive SD undo to allow single-copy CNV detection
cn_segmented <- segmentBins(cn, transformFun="log2", 
                            alpha = alpha_val, 
                            undo.splits = "none", 
                            min.width = 2)

# --- OUTPUT ---
message("[INFO] Writing output files...")

for (i in 1:ncol(cn_segmented)) {
  sample_name <- sampleNames(cn_segmented)[i]
  sample_name <- gsub("\\.bam$", "", basename(sample_name))
  single_sample <- cn_segmented[, i]
  
  exportManualResults(single_sample, sample_name, output_dir)
  
  pdf(file.path(output_dir, paste0(sample_name, "_report.pdf")), width = 12, height = 6)
  plot(single_sample)
  dev.off()
}

message("[SUCCESS] Analysis Complete.")
