#!/usr/bin/env python3

import argparse
import os
import math
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def main():
    parser = argparse.ArgumentParser(description="Filter cohort to the most stable 2/3 using mosdepth PCA.")
    parser.add_argument("-m", "--mosdepth", nargs='+', required=True, help="List of mosdepth .regions.bed.gz files")
    parser.add_argument("-b", "--bams", nargs='+', required=True, help="List of original BAM files")
    parser.add_argument("-o", "--output", required=True, help="Output file to write selected BAM paths")
    args = parser.parse_args()

    # Create a mapping of base names to full BAM paths for reliable output
    bam_map = {os.path.basename(b).replace('.bam', ''): b for b in args.bams}
    
    depth_data = {}
    
    print(">>> [PCA Triage] Loading mosdepth region arrays...")
    for f in args.mosdepth:
        sample_name = os.path.basename(f).replace('.regions.bed.gz', '')
        
        try:
            # Read chrom, start, end, and depth
            df = pd.read_csv(f, sep='\t', header=None, usecols=[0, 1, 2, 3])
            # Create a unique genomic coordinate index (e.g., "chr1_0_1000000")
            df['idx'] = df[0].astype(str) + '_' + df[1].astype(str) + '_' + df[2].astype(str)
            df.set_index('idx', inplace=True)
            
            # Store as a pandas series so it aligns automatically by index
            depth_data[sample_name] = df[3]
        except Exception as e:
            print(f"Error reading {f}: {e}")
            exit(1)

    # Construct dataFrame: rows = samples, columns = 1Mb bins
    df_depth = pd.DataFrame(depth_data).T
    
    # Safely drop any bins (columns) that weren't perfectly uniform across all BAMs
    df_depth.dropna(axis=1, inplace=True)

    if df_depth.shape[1] == 0:
        print("ERROR: [PCA Triage] All genomic bins were dropped due to misalignment across samples. Check BAM headers/contigs.")
        exit(1)

    print(f">>> [PCA Triage] Running analysis on {df_depth.shape[0]} samples across {df_depth.shape[1]} genomic bins...")
    
    n_samples = df_depth.shape[0]
    
    # Bulletproof bypass for small test cohorts
    if n_samples < 3:
        print(">>> [PCA Triage] Cohort too small for PCA (< 3 samples). Bypassing triage and selecting all samples.")
        best_samples = df_depth.index.tolist()
    else:
        # Standardize data and run PCA
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(df_depth)
        
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(scaled_data)

        # Calculate cluster centroid
        centroid = np.mean(pcs, axis=0)
        
        # Calculate Euclidean distance from each sample to the centroid
        distances = np.linalg.norm(pcs - centroid, axis=1)

        # Create results dataframe
        results = pd.DataFrame({
            'Sample': df_depth.index,
            'Distance': distances
        }).sort_values('Distance')

        # Determine dynamic cutoff using ceiling
        cutoff = math.floor(n_samples * 0.8)
        
        print(f">>> [PCA Triage] Selecting top {cutoff} out of {n_samples} samples.")
        best_samples = results.head(cutoff)['Sample'].tolist()

    # Map back to original BAM paths and write to output
    with open(args.output, 'w') as out_file:
        for sample in best_samples:
            if sample in bam_map:
                out_file.write(f"{bam_map[sample]}\n")
            else:
                print(f"Warning: Could not map {sample} back to original BAM.")

    print(f">>> [PCA Triage] Success. Selected BAMs written to {args.output}")

if __name__ == "__main__":
    main()
