#!/bin/bash

# ==============================================================================
# HYBRID CNV PIPELINE (split environments)
# Usage: ./run_hybrid_cnv.sh -i <bam> -r <roi.bed> -c <chrom> [-g ref] [-s bin|auto] [-q bin|skip] [-t threads] [-m method] [-p p-val] [-b blacklist]
# ==============================================================================

# Stop on error
set -e

# --- 1. CONFIGURATION & DEFAULTS ---

# Default values
REF_GENOME="hg19"
CNVKIT_BIN_SIZE="auto"
QDNA_BIN_SIZE_KB="50"
THREADS="4"
CHROMOSOME="all"
METHOD="hmm-germline"
THRESH=0.0001
BLACKLIST=1
ALPHA_VAL=0.01

# Environment names
ENV_CNVKIT="cnvkit_env"
ENV_QDNA="qdnaseq_env"

# Settings
RSCRIPT_NAME="qdnaseq_wes_pipeline.R"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
RSCRIPT_PATH="$SCRIPT_DIR/$RSCRIPT_NAME"

# --- 2. ARGUMENT PARSING ---

usage() {
    echo "Usage: $0 -i <input_pattern> -r <roi_bed> -c <chromosome> [-p <pooled_ref.cnn>] [OPTIONS]"
    echo ""
    echo "Arguments:"
    echo "  -i  Input BAM file. Use \"*\" (quoted) to process all .bam files in the current directory. (Required)"
    echo "  -r  ROI BED file (Required)"
    echo "  -c  Chromosome (e.g., \"22\" or \"all\") (Required)"
    echo "  -g  Reference genome (default: hg19)"
    echo "  -s  CNVkit bin size in bp or 'auto' (default: auto)"
    echo "  -q  QDNAseq bin size in kb, or 'skip' to bypass QDNAseq (default: 50)"
    echo "  -w  Alpha value for QDNAseq segmentation significance (default: 0.01)"
    echo "  -@  Threads (default: 4)"
    echo "  -m  Method of segmentation (default: cbs)"
    echo "  -t  Threshold/p-value of statistical significance (default: 0.0001)"
    echo "  -p  Path to a pooled reference (.cnn) file. If omitted, pipeline uses Flat Reference."
    echo "  -b  Blacklist exclusion: 1 for On, 0 for Off (default: 1)"
    echo "  -h  Show this help message"
    exit 1
}

INPUT_PATTERN=""
POOLED_REF_FILE=""

while getopts "i:r:c:g:s:q:@:m:t:p:b:h" opt; do
    case $opt in
        i) INPUT_PATTERN="$OPTARG" ;;
        r) ROI_BED="$OPTARG" ;;
        c) CHROMOSOME="$OPTARG" ;;
        g) REF_GENOME="$OPTARG" ;;
        s) CNVKIT_BIN_SIZE="$OPTARG" ;;
        q) QDNA_BIN_SIZE_KB="$OPTARG" ;;
        w ) ALPHA_VAL=$OPTARG ;;
        @) THREADS="$OPTARG" ;;
        m) METHOD="$OPTARG" ;;
        t) THRESH="$OPTARG" ;;
        p) POOLED_REF_FILE="$OPTARG" ;;
        b) BLACKLIST="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check mandatory arguments
if [[ -z "$INPUT_PATTERN" || -z "$ROI_BED" || -z "$CHROMOSOME" ]]; then
    echo "ERROR: Missing mandatory arguments -i, -r, or -c."
    usage
fi

# Evaluate the optional -p pooled reference flag
IS_POOLED_REF=0
if [[ -n "$POOLED_REF_FILE" ]]; then
    if [[ -f "$POOLED_REF_FILE" && "$POOLED_REF_FILE" == *.cnn ]]; then
        echo ">>> Valid pooled reference detected: $POOLED_REF_FILE"
        IS_POOLED_REF=1
    else
        echo "ERROR: Provided -p argument ('$POOLED_REF_FILE') is not a valid existing .cnn file."
        exit 1
    fi
else
    echo ">>> No pooled reference provided. Defaulting to Flat Reference mode."
fi

# Early verification: ensure R script exists if QDNAseq is not being skipped
if [[ "$QDNA_BIN_SIZE_KB" != "skip" ]]; then
    if [[ ! -f "$RSCRIPT_PATH" ]]; then
        echo "ERROR: R script '$RSCRIPT_NAME' not found in $SCRIPT_DIR."
        echo "Execution aborted before starting."
        exit 1
    fi
fi

# --- GATHER INPUT BAM(S) ---
if [[ "$INPUT_PATTERN" == "*" ]]; then
    shopt -s nullglob
    FILES=(*.bam)
    shopt -u nullglob
else
    FILES=("$INPUT_PATTERN")
fi

if [ ${#FILES[@]} -eq 0 ]; then
    echo "ERROR: No input BAM files found to process."
    exit 1
fi

# ==============================================================================
# 3. CONDA INITIALIZATION
# ==============================================================================

if command -v conda &> /dev/null; then
    CONDA_BASE=$(conda info --base)
    source "$CONDA_BASE/etc/profile.d/conda.sh"
else
    # Fallback paths
    source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || \
    source "$HOME/anaconda3/etc/profile.d/conda.sh" 2>/dev/null
fi

# Directory setup
OUT_CNVKIT="./cnvkit_results"
OUT_QDNA="./qdnaseq_results"
REF_DIR="./reference_${REF_GENOME}"
mkdir -p "$OUT_CNVKIT" "$OUT_QDNA" "$REF_DIR"

# ==============================================================================
# 4. PHASE 1: CNVkit (Exome/Targeted Mode)
# ==============================================================================

echo ">>> Activating Environment: $ENV_CNVKIT"
conda activate "$ENV_CNVKIT"

echo "=========================================="
echo " RUNNING CNVkit (Exome/Hybrid Capture)"
echo " Inputs: BAM=$INPUT_PATTERN | ROI=$ROI_BED"
echo " Config: Ref=$REF_GENOME | Chr=$CHROMOSOME | Bin=$CNVKIT_BIN_SIZE"
echo " Params: Method=$METHOD | Threshold=$THRESH | Blacklist=$BLACKLIST"
echo "=========================================="

# --- 4.1 Reference preparation (chromosome specific) ---

if [ "$IS_POOLED_REF" -eq 0 ]; then
    if [ "$REF_GENOME" == "hs37d5" ]; then
        REF_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
        REF_RAW="${REF_DIR}/hs37d5_raw.fa.gz"

        if [ "$CHROMOSOME" == "all" ]; then
            REF_CLEAN="${REF_DIR}/hs37d5.fa"
        else
            REF_CLEAN="${REF_DIR}/hs37d5_chr${CHROMOSOME}.fa"
        fi

        if [ ! -f "$REF_CLEAN" ]; then
            if [ ! -f "$REF_RAW" ]; then
                echo "   Downloading HS37D5 Raw Reference..."
                wget -q -O "$REF_RAW" "$REF_URL"
            fi

            echo "   Extracting reference sequence..."
            if [ "$CHROMOSOME" == "all" ]; then
                zcat "$REF_RAW" > "$REF_CLEAN"
            else
                zcat "$REF_RAW" | awk -v chr="$CHROMOSOME" 'BEGIN {P=0} /^>/ {if ($1 == ">"chr) P=1; else P=0} {if (P) print $0}' > "$REF_CLEAN"
            fi
            
            echo "   Indexing Reference..."
            samtools faidx "$REF_CLEAN"
        fi

    elif [ "$CHROMOSOME" == "all" ]; then
        REF_URL="https://hgdownload.cse.ucsc.edu/goldenpath/${REF_GENOME}/bigZips/${REF_GENOME}.fa.gz"
        REF_RAW="${REF_DIR}/${REF_GENOME}_raw.fa.gz"
        REF_CLEAN="${REF_DIR}/${REF_GENOME}.fa"
        
        if [ ! -f "$REF_CLEAN" ]; then
            echo "   Downloading Reference ($REF_GENOME)..."
            wget -q -O "$REF_RAW" "$REF_URL"
            zcat "$REF_RAW" > "$REF_CLEAN"
            rm "$REF_RAW"
            samtools faidx "$REF_CLEAN"
        fi
    else
        REF_URL="https://hgdownload.cse.ucsc.edu/goldenpath/${REF_GENOME}/chromosomes/chr${CHROMOSOME}.fa.gz"
        REF_RAW="${REF_DIR}/chr${CHROMOSOME}_raw.fa.gz"
        REF_CLEAN="${REF_DIR}/chr${CHROMOSOME}.fa"

        if [ ! -f "$REF_CLEAN" ]; then
            echo "   Downloading Reference ($REF_GENOME chr$CHROMOSOME)..."
            wget -q -O "$REF_RAW" "$REF_URL"
            zcat "$REF_RAW" > "$REF_CLEAN"
            rm "$REF_RAW"
            samtools faidx "$REF_CLEAN"
        fi
    fi
fi

# Determine REF_STYLE (Needed for Blacklist and ROI parsing)
if [ "$IS_POOLED_REF" -eq 1 ]; then
    FIRST_CONTIG=$(awk 'NR==2 {print $1}' "$POOLED_REF_FILE")
    if [[ "$FIRST_CONTIG" == chr* ]]; then REF_STYLE="chr"; else REF_STYLE="plain"; fi
else
    REF_HEADER=$(head -n 1 "$REF_CLEAN")
    CLEAN_HEADER="${REF_HEADER#>}" 
    if [[ "$CLEAN_HEADER" == chr* ]]; then REF_STYLE="chr"; else REF_STYLE="plain"; fi
fi

# --- 4.2 Blacklist preparation ---

if [ "$BLACKLIST" -eq 1 ]; then
    BLACKLIST_DIR="./blacklist"
    mkdir -p "$BLACKLIST_DIR"
    BL_FETCH_REF=$([ "$REF_GENOME" == "hs37d5" ] && echo "hg19" || echo "$REF_GENOME")
    BLACKLIST_URL="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/${BL_FETCH_REF}-blacklist.v2.bed.gz"
    BLACKLIST_RAW="${BLACKLIST_DIR}/${BL_FETCH_REF}_blacklist.bed.gz"
    BLACKLIST_CLEAN="${BLACKLIST_DIR}/${REF_GENOME}_${REF_STYLE}_blacklist_clean.bed"

    # Download the raw file ONLY if it is not already present
    if [ ! -f "$BLACKLIST_RAW" ]; then
        echo "   Downloading Blacklist (Source: $BL_FETCH_REF)..."
        wget -q -O "$BLACKLIST_RAW" "$BLACKLIST_URL"
    fi

    # Format the clean file ONLY if it is not already present
    if [ ! -f "$BLACKLIST_CLEAN" ]; then
        echo "   Formatting Blacklist for $REF_GENOME ($REF_STYLE)..."
        if [ "$REF_STYLE" == "plain" ]; then
            zcat "$BLACKLIST_RAW" | sed 's/^chr//g' > "$BLACKLIST_CLEAN"
        else
            zcat "$BLACKLIST_RAW" > "$BLACKLIST_CLEAN"
        fi
    fi
fi

if [ "$BLACKLIST" -eq 1 ] && [ -f "$BLACKLIST_CLEAN" ]; then
    ACCESS_BED="${REF_DIR}/access_${REF_GENOME}_chr${CHROMOSOME}.bed"
    if [ "$IS_POOLED_REF" -eq 0 ]; then
        if [ ! -f "$ACCESS_BED" ]; then
            echo "   Generating CNVkit Access regions (excluding Blacklist)..."
            cnvkit.py access "$REF_CLEAN" -x "$BLACKLIST_CLEAN" -o "$ACCESS_BED"
        fi
        ACCESS_ARG="-g $ACCESS_BED"
    else
        ACCESS_ARG=""
    fi
else
    ACCESS_ARG=""
fi

# --- 4.3 Target, antitarget & ROI preparation ---

FILTERED_BED="${REF_DIR}/roi_filtered_${CHROMOSOME}.bed"
SHARED_TARGETS="${REF_DIR}/targets_${REF_GENOME}_chr${CHROMOSOME}_${CNVKIT_BIN_SIZE}bp.bed"
SHARED_ANTITARGETS="${REF_DIR}/antitargets_${REF_GENOME}_chr${CHROMOSOME}_${CNVKIT_BIN_SIZE}bp.bed"
SHARED_FLAT_REF="${REF_DIR}/flat_ref_${REF_GENOME}_chr${CHROMOSOME}_${CNVKIT_BIN_SIZE}bp.cnn"

# REF_STYLE is determined above.

# Filter BED for specific chromosome
echo "   Filtering BED file for Chr: $CHROMOSOME..."
awk -v style="$REF_STYLE" -v target="$CHROMOSOME" '
{
    gsub(/\r/, "", $0)
    c=$1; s=$2; e=$3;
    sub(/^chr/, "", c)
    target_clean = target
    sub(/^chr/, "", target_clean)
    
    if (target == "all" || c == target_clean) {
        if (style == "chr") { print "chr" c "\t" s "\t" e } 
        else { print c "\t" s "\t" e }
    }
}' "$ROI_BED" > "$FILTERED_BED"

if [ "$BLACKLIST" -eq 1 ] && [ -f "$BLACKLIST_CLEAN" ]; then
    echo "   Subtracting Blacklist regions from Target BED..."
    bedtools subtract -a "$FILTERED_BED" -b "$BLACKLIST_CLEAN" > "${FILTERED_BED}_clean"
    mv "${FILTERED_BED}_clean" "$FILTERED_BED"
fi

if [ ! -s "$FILTERED_BED" ]; then
    echo "ERROR: No matching regions found in BED for chromosome $CHROMOSOME."
    exit 1
fi

if [ "$IS_POOLED_REF" -eq 1 ]; then
    REF_TO_USE="$POOLED_REF_FILE"
    SHARED_TARGETS="${REF_DIR}/extracted_targets_${CHROMOSOME}.bed"
    SHARED_ANTITARGETS="${REF_DIR}/extracted_antitargets_${CHROMOSOME}.bed"
    
    if [ ! -f "$SHARED_TARGETS" ] || [ ! -f "$SHARED_ANTITARGETS" ]; then
        echo "   Extracting target/antitarget BED from pooled reference: $POOLED_REF_FILE..."
        if [ "$CHROMOSOME" == "all" ]; then
            awk 'BEGIN {FS="\t"; OFS="\t"} NR>1 {if ($4 != "Antitarget") print $1, $2, $3, $4}' "$POOLED_REF_FILE" > "${SHARED_TARGETS}_raw"
            awk 'BEGIN {FS="\t"; OFS="\t"} NR>1 {if ($4 == "Antitarget") print $1, $2, $3, $4}' "$POOLED_REF_FILE" > "${SHARED_ANTITARGETS}_raw"
        else
            awk -v chr="$CHROMOSOME" 'BEGIN {FS="\t"; OFS="\t"} NR>1 {
                if (($1 == chr || $1 == "chr"chr) && $4 != "Antitarget") print $1, $2, $3, $4
            }' "$POOLED_REF_FILE" > "${SHARED_TARGETS}_raw"
            awk -v chr="$CHROMOSOME" 'BEGIN {FS="\t"; OFS="\t"} NR>1 {
                if (($1 == chr || $1 == "chr"chr) && $4 == "Antitarget") print $1, $2, $3, $4
            }' "$POOLED_REF_FILE" > "${SHARED_ANTITARGETS}_raw"
        fi

        echo "   Applying ROI and Blacklist filters to extracted pooled targets..."
        bedtools intersect -u -a "${SHARED_TARGETS}_raw" -b "$FILTERED_BED" > "$SHARED_TARGETS"
        
        if [ "$BLACKLIST" -eq 1 ] && [ -f "$BLACKLIST_CLEAN" ]; then
            bedtools subtract -A -a "${SHARED_ANTITARGETS}_raw" -b "$BLACKLIST_CLEAN" > "$SHARED_ANTITARGETS"
        else
            cp "${SHARED_ANTITARGETS}_raw" "$SHARED_ANTITARGETS"
        fi
        
        rm -f "${SHARED_TARGETS}_raw" "${SHARED_ANTITARGETS}_raw"
    fi
else
    REF_TO_USE="$SHARED_FLAT_REF"

    if [ ! -f "$SHARED_TARGETS" ] || [ ! -f "$SHARED_ANTITARGETS" ]; then
        if [[ "${CNVKIT_BIN_SIZE,,}" == "auto" ]]; then
            echo "   Calculating optimal bin size via autobin..."
            cnvkit.py autobin "${FILES[0]}" -t "$FILTERED_BED" $ACCESS_ARG
            
            FIRST_BASE=$(basename "${FILES[0]}" .bam)
            if [ -f "${FIRST_BASE}.target.bed" ] && [ -f "${FIRST_BASE}.antitarget.bed" ]; then
                mv "${FIRST_BASE}.target.bed" "$SHARED_TARGETS"
                mv "${FIRST_BASE}.antitarget.bed" "$SHARED_ANTITARGETS"
                echo "   Autobinning complete. Targets and Antitargets saved."
            else
                echo "ERROR: autobin failed to generate targets/antitargets."
                exit 1
            fi
        else
            echo "   Creating CNVkit Targets with fixed bin size: $CNVKIT_BIN_SIZE..."
            cnvkit.py target "$FILTERED_BED" --split -a "$CNVKIT_BIN_SIZE" -o "$SHARED_TARGETS"
            
            echo "   Creating CNVkit Antitargets for Exome Mode..."
            cnvkit.py antitarget "$SHARED_TARGETS" $ACCESS_ARG -o "$SHARED_ANTITARGETS"
        fi
    fi

    if [ ! -f "$SHARED_FLAT_REF" ]; then
        echo "   Creating CNVkit Flat Exome Reference (Targets + Antitargets)..."
        cnvkit.py reference -o "$SHARED_FLAT_REF" -f "$REF_CLEAN" -t "$SHARED_TARGETS" -a "$SHARED_ANTITARGETS"
    fi
fi

# --- 4.4 Execution ---

for INPUT_BAM in "${FILES[@]}"; do
    BASENAME=$(basename "$INPUT_BAM" .bam)
    SAMPLE_DIR="${OUT_CNVKIT}/${BASENAME}"
    mkdir -p "$SAMPLE_DIR"

    case "$BASENAME" in
        WGS001|WGS002|WGS004|WGS005|WGS006|WGS009|WGS010|WGS015|WGS016|WGS021)
            SAMPLE_SEX="female"
            if [ -n "$POOLED_REF" ]; then
                ACTIVE_REF="${POOLED_REF%.cnn}_female.cnn"
            else
                ACTIVE_REF="$FLAT_REF"
            fi
            ;;
        WGS003|WGS007|WGS008|WGS011|WGS012|WGS013|WGS014|WGS017|WGS018|WGS019|WGS020|WGS022|WGS023)
            SAMPLE_SEX="male"
            if [ -n "$POOLED_REF" ]; then
                ACTIVE_REF="${POOLED_REF%.cnn}_male.cnn"
            else
                ACTIVE_REF="$FLAT_REF"
            fi
            ;;
        *)
            echo "   [WARNING] Unknown sex for $BASENAME. Assuming female baseline."
            SAMPLE_SEX="female"
            if [ -n "$POOLED_REF" ]; then
                ACTIVE_REF="${POOLED_REF%.cnn}_female.cnn"
            else
                ACTIVE_REF="$FLAT_REF"
            fi
            ;;
    esac

    # Fail-fast if the assigned pooled reference is missing
    if [ -n "$POOLED_REF" ] && [ ! -f "$ACTIVE_REF" ]; then
        echo "ERROR: Sex-specific pooled reference not found at $ACTIVE_REF"
        exit 1
    fi

    # -----------------------------------------------------
    # CNVKIT EXECUTION (With Sex Flags & Noise Filters)
    # -----------------------------------------------------
    
    echo "   [CNVkit] Calculating Target/Antitarget coverage..."
    cnvkit.py coverage "$BAM" "$SHARED_TARGETS" -p "$THREADS" -o "${SAMPLE_OUT_DIR}/${BASENAME}.targetcoverage.cnn"
    
    if [ "$QDNA_BIN_SIZE_KB" != "skip" ]; then
        cnvkit.py coverage "$BAM" "$SHARED_ANTITARGETS" -p "$THREADS" -o "${SAMPLE_OUT_DIR}/${BASENAME}.antitargetcoverage.cnn"
        
        echo "   [CNVkit] Fixing log2 ratios (Hybrid Mode) against $SAMPLE_SEX reference..."
        cnvkit.py fix "${SAMPLE_OUT_DIR}/${BASENAME}.targetcoverage.cnn" \
                      "${SAMPLE_OUT_DIR}/${BASENAME}.antitargetcoverage.cnn" \
                      "$ACTIVE_REF" -o "${SAMPLE_OUT_DIR}/${BASENAME}.cnr"
    else
        echo "   [CNVkit] Fixing log2 ratios (Targeted Mode) against $SAMPLE_SEX reference..."
        cnvkit.py fix "${SAMPLE_OUT_DIR}/${BASENAME}.targetcoverage.cnn" \
                      "$ACTIVE_REF" -o "${SAMPLE_OUT_DIR}/${BASENAME}.cnr"
    fi

    echo "   [CNVkit] Running segmentation ($METHOD)..."
    if [[ "$METHOD" == *"hmm"* ]]; then
        # HMM Branch (Excludes smooth-cbs and threshold)
        cnvkit.py segment "${SAMPLE_OUT_DIR}/${BASENAME}.cnr" \
            -m "$METHOD" \
            --drop-low-coverage \
            --drop-outliers 10 \
            -p "$THREADS" \
            -o "${SAMPLE_OUT_DIR}/${BASENAME}_${METHOD}.cns"
    else
        # CBS Branch
        cnvkit.py segment "${SAMPLE_OUT_DIR}/${BASENAME}.cnr" \
            -m "$METHOD" \
            --threshold "$THRESH" \
            --smooth-cbs \
            --drop-low-coverage \
            --drop-outliers 10 \
            -p "$THREADS" \
            -o "${SAMPLE_OUT_DIR}/${BASENAME}_${METHOD}.cns"
    fi

    echo "   [CNVkit] Calling discrete copy numbers and merging bins..."
    # CRITICAL FIX: Evaluates segment states against proper sex baseline
    cnvkit.py call "${SAMPLE_OUT_DIR}/${BASENAME}_${METHOD}.cns" \
        -x "$SAMPLE_SEX" \
        --filter cn \
        -o "${SAMPLE_OUT_DIR}/${BASENAME}_${METHOD}.call.cns"

    echo "   [CNVkit] Filtering out hyper-segmented noise (<10kb)..."
    awk -F'\t' 'NR==1 {print $0; next} ($3 - $2 >= 10000)' \
        "${SAMPLE_OUT_DIR}/${BASENAME}_${METHOD}.call.cns" \
        > "${SAMPLE_OUT_DIR}/${BASENAME}_${METHOD}.call.filtered.cns"

    echo "   [CNVkit] Exporting to BED..."
    CNVKIT_RESULTS_BED="${OUT_CNVKIT}/${BASENAME}_${METHOD}_cnvkit_results.bed"

    # Export to temporary raw bed keeping original contig names for correct bedtools subtract
    tail -n +2 "${SAMPLE_OUT_DIR}/${BASENAME}_${METHOD}.call.filtered.cns" | \
    awk '{print $1 "\t" $2 "\t" $3 "\tcnvkit\t" $5}' \
    > "${CNVKIT_RESULTS_BED}_temp1"

    if [ "$BLACKLIST" -eq 1 ] && [ -f "$BLACKLIST_CLEAN" ]; then
        echo "   [bedtools] Removing blacklist overlaps from final BED..."
        bedtools subtract -a "${CNVKIT_RESULTS_BED}_temp1" -b "$BLACKLIST_CLEAN" > "${CNVKIT_RESULTS_BED}_temp2"
    else
        cp "${CNVKIT_RESULTS_BED}_temp1" "${CNVKIT_RESULTS_BED}_temp2"
    fi      

    # Write directly to the main directory, stripping 'chr' for normalization
    awk 'BEGIN{OFS="\t"} {
        c=$1; 
        sub(/^chr/, "", c); 
        print c, $2, $3, $4, $5
    }' "${CNVKIT_RESULTS_BED}_temp2" > "$CNVKIT_RESULTS_BED"

    rm -f "${CNVKIT_RESULTS_BED}_temp1" "${CNVKIT_RESULTS_BED}_temp2"
done

echo ">>> Deactivating $ENV_CNVKIT"
conda deactivate

# ==============================================================================
# 5. PHASE 2: QDNAseq (WGS)
# ==============================================================================

if [[ "$QDNA_BIN_SIZE_KB" != "skip" ]]; then

    echo ">>> Generating 1kb padded exclusion mask for QDNAseq..."
    PADDED_EXCLUDE_BED="${OUT_QDNA}/targets_padded_1kb.bed"
    
    # Use awk instead of bedtools to avoid .genome file dependencies.
    # Expands the target by 1000bp on both sides (floor of 0 to prevent negative coords).
    awk -v pad=1000 'BEGIN{OFS="\t"} {
        start = $2 - pad;
        if (start < 0) start = 0;
        print $1, start, $3 + pad
    }' "$ROI_BED" > "$PADDED_EXCLUDE_BED"

    echo ">>> Activating Environment: $ENV_QDNA"
    conda activate "$ENV_QDNA"

    echo "=========================================="
    echo " RUNNING QDNAseq (WGS) IN BATCH MODE"
    echo " Bin Size: ${QDNA_BIN_SIZE_KB}kb | Chromosome: ${CHROMOSOME}"
    echo " Threads: ${THREADS:-4} | Total Samples: ${#FILES[@]}"
    echo "=========================================="

    echo "   [QDNAseq] Running R Pipeline..."
    
    # --- Pass arguments: OutDir, BinSize, ROI, Chromosome, Threads, and ALL BAMs ---
    Rscript "$RSCRIPT_PATH" \
        "$OUT_QDNA" \
        "$QDNA_BIN_SIZE_KB" \
        "$PADDED_EXCLUDE_BED" \
        "$CHROMOSOME" \
        "$THREADS" \
        "$ALPHA_VAL" \
        "${FILES[@]}"

    echo ">>> Deactivating $ENV_QDNA"
    conda deactivate
else
    echo "=========================================="
    echo " SKIPPING QDNAseq (WGS) Phase (-q skip)"
    echo "=========================================="
fi

# ==============================================================================
# 6. MERGE & POST-PROCESSING
# ==============================================================================

echo "=========================================="
echo " PROCESSING & MERGING"
echo "=========================================="

for INPUT_BAM in "${FILES[@]}"; do
    BASENAME=$(basename "$INPUT_BAM" .bam)
    CNVKIT_RESULTS_BED="${OUT_CNVKIT}/${BASENAME}_${METHOD}_cnvkit_results.bed"

    if [[ "$QDNA_BIN_SIZE_KB" != "skip" ]]; then
        # FIX: Explicitly target the deterministic file path
        QDNA_RAW_BED="${OUT_QDNA}/${BASENAME}_segments.bed"

        if [ ! -f "$QDNA_RAW_BED" ]; then
            echo "Error: QDNAseq output BED file not found for $BASENAME at $QDNA_RAW_BED."
            exit 1
        fi

        TEMP_BED="${QDNA_RAW_BED%.bed}_temp.bed"
        FINAL_QDNA_BED="${OUT_QDNA}/${BASENAME}_qdna.bed"

        echo "   Formatting QDNAseq output for $BASENAME..."
        awk 'BEGIN{OFS="\t"} {
            c=$1; 
            sub(/^chr/, "", c);
            print c, $2, $3, "qdnaseq", $5
        }' "$QDNA_RAW_BED" > "$TEMP_BED"

        mv "$TEMP_BED" "$FINAL_QDNA_BED"

        TOTAL_BED="./${BASENAME}_total.bed"

        echo "   Merging into $TOTAL_BED..."
        cat "$CNVKIT_RESULTS_BED" "$FINAL_QDNA_BED" > "$TOTAL_BED"

        # Simple sort
        sort -k1,1 -k2,2n "$TOTAL_BED" > "${TOTAL_BED}.tmp" && mv "${TOTAL_BED}.tmp" "$TOTAL_BED"
    else
        echo "   QDNAseq pipeline was bypassed for $BASENAME."
        echo "   Copying CNVkit output to the main directory..."
        cp "$CNVKIT_RESULTS_BED" "./${BASENAME}_${METHOD}_cnvkit_results.bed"
        echo "   Final output saved to ./${BASENAME}_${METHOD}_cnvkit_results.bed"
    fi
done

echo "Pipeline Complete."
