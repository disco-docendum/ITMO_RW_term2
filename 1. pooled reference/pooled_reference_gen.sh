#!/bin/bash

# ==============================================================================
# CNVkit Pooled Reference Generator (with Blacklist, WGS fallback & Auto-FASTA)
# Usage: bash ./pooled_ref_gen.sh -o <out_ref.cnn> [-r <roi.bed>] [-g <genome>] [-d <bam_dir>] [BAM_FILES...]
# ==============================================================================

set -e

# CONFIGURATION & DEFAULTS
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
BIN_SIZE="auto"
ANTITARGET_BIN_SIZE=""
THREADS="4"
ENV_CNVKIT="cnvkit_env"
ENV_MOSDEPTH="mosdepth_env"
OUT_DIR="./pooled_reference_build"

REF_GENOME="hg19"
BLACKLIST=1
BAM_DIR=""
ROI_BED=""

usage() {
    echo "Usage: $0 -o <output_ref.cnn> [OPTIONS] [file1.bam file2.bam ...]"
    echo ""
    echo "Arguments:"
    echo "  -o  Output reference filename (e.g., pooled_reference.cnn) (required)"
    echo "  -r  ROI BED file. If OMITTED, script defaults to WGS mode (optional)"
    echo "  -d  Directory containing normal/control BAM files (optional if BAM files are listed)"
    echo "  -g  Reference genome to download (e.g., hg19, hg38, hs37d5) (default: hg19)"
    echo "  -b  Blacklist exclusion: 1 for On, 0 for Off (default: 1)"
    echo "  -s  CNVkit target bin size in bp or 'auto' (default: auto). For vdWGS mode, must be manual."
    echo "  -a  CNVkit antitarget bin size in bp. Required for vdWGS hybrid mode."
    echo "  -@  Threads (default: 4)"
    echo "  -h  Show this help message"
    exit 1
}

while getopts "d:r:o:g:b:s:a:@:h" opt; do
    case $opt in
        d) BAM_DIR="$OPTARG" ;;
        r) ROI_BED="$OPTARG" ;;
        o) OUT_REF="$OPTARG" ;;
        g) REF_GENOME="$OPTARG" ;;
        b) BLACKLIST="$OPTARG" ;;
        s) BIN_SIZE="$OPTARG" ;;
        a) ANTITARGET_BIN_SIZE="$OPTARG" ;;
        @) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

shift $((OPTIND-1))

# Fail-fast: If ROI_BED is provided (vdWGS mode) but ANTITARGET_BIN_SIZE is empty
if [ -n "$ROI_BED" ] && [ -z "$ANTITARGET_BIN_SIZE" ]; then
    echo "ERROR: ROI_BED was provided (vdWGS mode triggered), but ANTITARGET_BIN_SIZE is empty."
    echo "Please provide a valid integer for the antitarget bin size (e.g., 50000)."
    exit 1
fi

# Check mandatory structural arguments
if [[ -z "$OUT_REF" ]]; then
    echo "ERROR: Missing mandatory argument (-o)."
    usage
fi

# Determine pipeline mode (WGS vs Targeted)
IS_WGS=0
if [[ -z "$ROI_BED" ]]; then
    echo ">>> No ROI BED provided. Defaulting to Whole Genome Sequencing (WGS) mode."
    IS_WGS=1
else
    echo ">>> ROI BED provided. Defaulting to vdWGS (Hybrid) mode."
    
    if [[ "${BIN_SIZE,,}" == "auto" ]] || [[ -z "$ANTITARGET_BIN_SIZE" ]]; then
        echo "WARNING: vdWGS hybrid pipeline is active, but manual target (-s) and antitarget (-a) bin sizes are missing or set to auto."
        echo "For this mode, only manual bin size of targets and antitargets should be accepted."
        echo "Exiting script."
        exit 1
    fi
fi

# ==========================================
# 1. COMPILE INPUT BAM LIST
# ==========================================

BAM_LIST=()

if [[ -n "$BAM_DIR" ]]; then
    shopt -s nullglob
    for b in "$BAM_DIR"/*.bam; do
        BAM_LIST+=("$b")
    done
    shopt -u nullglob
fi

for b in "$@"; do
    BAM_LIST+=("$b")
done

if [ ${#BAM_LIST[@]} -eq 0 ]; then
    echo "ERROR: No BAM files found."
    usage
fi

# ==========================================
# 2. CONDA INITIALIZATION & ENV CHECKS
# ==========================================

# Locate and source conda
if command -v conda &> /dev/null; then
    CONDA_BASE=$(conda info --base)
    source "$CONDA_BASE/etc/profile.d/conda.sh"
else
    source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || \
    source "$HOME/anaconda3/etc/profile.d/conda.sh" 2>/dev/null
fi

# Verify conda is actually available
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda is not installed or not in PATH."
    exit 1
fi

# Verify both environments exist before starting the pipeline
if ! conda info --envs | grep -q "^${ENV_MOSDEPTH} "; then
    echo "ERROR: Conda environment '$ENV_MOSDEPTH' not found. Please create it first."
    exit 1
fi

if ! conda info --envs | grep -q "^${ENV_CNVKIT} "; then
    echo "ERROR: Conda environment '$ENV_CNVKIT' not found. Please create it first."
    exit 1
fi

# Fast dry-run to see if conda activation is actually necessary
NEEDS_INDEX=0
for BAM in "${BAM_LIST[@]}"; do
    if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
        NEEDS_INDEX=1
        break # The moment we find one missing index, stop checking and set the flag
    fi
done

# Only load the environment if work needs to be done
if [ "$NEEDS_INDEX" -eq 1 ]; then
    echo ">>> Activating Environment: $ENV_CNVKIT (for SAMtools pre-flight)"
    conda activate "$ENV_CNVKIT"

    for BAM in "${BAM_LIST[@]}"; do
        if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
            echo "   [SAMtools] Index missing for $(basename "$BAM"). Generating .bai index..."
            samtools index -@ "$THREADS" "$BAM"
        fi
    done
    
    conda deactivate
else
    echo "   [Pre-flight] All BAM indices found. Skipping SAMtools activation."
fi
echo ">>> Activating Environment: $ENV_MOSDEPTH (for PCA Triage)"
conda activate "$ENV_MOSDEPTH"

# Verify required tools for triage are present
if ! command -v mosdepth &> /dev/null; then echo "ERROR: 'mosdepth' not found in $ENV_MOSDEPTH."; exit 1; fi
if ! command -v python3 &> /dev/null; then echo "ERROR: 'python3' not found in $ENV_MOSDEPTH."; exit 1; fi

mkdir -p "$OUT_DIR"

# ==========================================
# 2.5. COHORT TRIAGE (MOSDEPTH + PCA)
# ==========================================

echo "=========================================="
echo " STAGE 0.5: Selecting optimal cohort (PCA)"
echo "=========================================="

TRIAGE_DIR="${OUT_DIR}/triage"
mkdir -p "$TRIAGE_DIR"
MOSDEPTH_FILES=()

# Run mosdepth rapidly in 1Mb bins
for BAM in "${BAM_LIST[@]}"; do
    BASENAME=$(basename "$BAM" .bam)
    REGIONS_FILE="${TRIAGE_DIR}/${BASENAME}.regions.bed.gz"
    
    if [ ! -f "$REGIONS_FILE" ]; then
        echo "   [mosdepth] Profiling $BASENAME..."
        # --fast-mode assumes no overlapping mates, perfectly fine for 1Mb binning
        mosdepth -n --fast-mode -b 1000000 "${TRIAGE_DIR}/${BASENAME}" "$BAM"
    fi
    MOSDEPTH_FILES+=("$REGIONS_FILE")
done

# Execute python PCA triage
BEST_BAMS_LIST="${TRIAGE_DIR}/selected_bams.txt"
python3 "$SCRIPT_DIR/pca_triage.py" \
    -m "${MOSDEPTH_FILES[@]}" \
    -b "${BAM_LIST[@]}" \
    -o "$BEST_BAMS_LIST"

if [ ! -f "$BEST_BAMS_LIST" ]; then
    echo "ERROR: PCA Triage failed."
    exit 1
fi

# Overwrite the BAM_LIST array with the filtered selection
BAM_LIST=()
while IFS= read -r line || [ -n "$line" ]; do
    BAM_LIST+=("$line")
done < "$BEST_BAMS_LIST"

# FAIL-FAST: Ensure samples survived
if [ ${#BAM_LIST[@]} -eq 0 ]; then
    echo "ERROR: Zero BAM files survived PCA triage or file mapping failed."
    exit 1
fi

echo "   [Triage Complete] Proceeding to reference generation with ${#BAM_LIST[@]} samples."

# Switch environments for the remainder of the pipeline
echo ">>> Switching Environment: $ENV_CNVKIT (for Reference Building)"
conda deactivate
conda activate "$ENV_CNVKIT"

# Verify tools for the main pipeline

if ! command -v samtools &> /dev/null; then echo "ERROR: 'samtools' not found in $ENV_CNVKIT."; exit 1; fi
if ! command -v cnvkit.py &> /dev/null; then echo "ERROR: 'cnvkit.py' not found in $ENV_CNVKIT."; exit 1; fi
if [ "$BLACKLIST" -eq 1 ] && ! command -v bedtools &> /dev/null; then echo "ERROR: 'bedtools' not found in $ENV_CNVKIT. Required for blacklist subtraction."; exit 1; fi


# ==========================================
# 3. FASTA PREPARATION (STAGE 0)
# ==========================================

echo "=========================================="
echo " STAGE 0: Reference FASTA preparation"
echo "=========================================="

REF_DIR="${OUT_DIR}/reference"
mkdir -p "$REF_DIR"

if [ "$REF_GENOME" == "hs37d5" ]; then
    REF_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
    REF_RAW="${REF_DIR}/hs37d5_raw.fa.gz"
    REF_FASTA="${REF_DIR}/hs37d5.fa"

    if [ ! -f "$REF_FASTA" ]; then
        if [ ! -f "$REF_RAW" ]; then
            echo "   Downloading HS37D5 Raw Reference..."
            wget -q -O "$REF_RAW" "$REF_URL"
        fi
        echo "   Extracting reference sequence..."
        zcat "$REF_RAW" > "$REF_FASTA"
        rm -f "$REF_RAW"
    fi
else
    REF_URL="https://hgdownload.cse.ucsc.edu/goldenpath/${REF_GENOME}/bigZips/${REF_GENOME}.fa.gz"
    REF_RAW="${REF_DIR}/${REF_GENOME}_raw.fa.gz"
    REF_FASTA="${REF_DIR}/${REF_GENOME}.fa"

    if [ ! -f "$REF_FASTA" ]; then
        echo "   Downloading Reference ($REF_GENOME)..."
        wget -q -O "$REF_RAW" "$REF_URL"
        echo "   Extracting reference sequence..."
        zcat "$REF_RAW" > "$REF_FASTA"
        rm -f "$REF_RAW"
    fi
fi

if [ ! -f "${REF_FASTA}.fai" ]; then
    echo "   Indexing Reference..."
    samtools faidx "$REF_FASTA"
fi

# ==========================================
# 4. BLACKLIST & ACCESS PREPARATION
# ==========================================

ACCESS_BED="${OUT_DIR}/access_${REF_GENOME}.bed"
ACCESS_ARG="-g $ACCESS_BED"
BLACKLIST_CLEAN=""

echo "=========================================="
echo " STAGE 1: Blacklist & access preparation"
echo "=========================================="

# 4A. Handle blacklist download if enabled
if [ "$BLACKLIST" -eq 1 ]; then
    BLACKLIST_DIR="${OUT_DIR}/blacklist"
    mkdir -p "$BLACKLIST_DIR"
    
    BL_FETCH_REF=$([ "$REF_GENOME" == "hs37d5" ] && echo "hg19" || echo "$REF_GENOME")
    BLACKLIST_URL="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/${BL_FETCH_REF}-blacklist.v2.bed.gz"
    BLACKLIST_RAW="${BLACKLIST_DIR}/${BL_FETCH_REF}_blacklist.bed.gz"
    BLACKLIST_CLEAN="${BLACKLIST_DIR}/${REF_GENOME}_blacklist_clean.bed"

    if [ ! -f "$BLACKLIST_CLEAN" ]; then
        echo "   Downloading Blacklist (Source: $BL_FETCH_REF)..."
        wget -q -O "$BLACKLIST_RAW" "$BLACKLIST_URL"
        echo "   Formatting Blacklist for $REF_GENOME..."
        if [ "$REF_GENOME" == "hs37d5" ]; then
            zcat "$BLACKLIST_RAW" | sed 's/^chr//g' > "$BLACKLIST_CLEAN"
        else
            zcat "$BLACKLIST_RAW" > "$BLACKLIST_CLEAN"
        fi
        rm -f "$BLACKLIST_RAW"
    fi
fi

# 4B. ALWAYS generate access BED (with or without blacklist)
if [ ! -f "$ACCESS_BED" ]; then
    if [ "$BLACKLIST" -eq 1 ] && [ -f "$BLACKLIST_CLEAN" ]; then
        echo "   Calculating accessible regions from FASTA (Excluding Blacklist)..."
        cnvkit.py access "$REF_FASTA" -x "$BLACKLIST_CLEAN" -o "$ACCESS_BED"
    else
        echo "   Calculating accessible regions from FASTA (No Blacklist applied)..."
        cnvkit.py access "$REF_FASTA" -o "$ACCESS_BED"
    fi
fi

# 4C. Subtract Blacklist from ROI (targeted mode only)
if [ "$IS_WGS" -eq 0 ] && [ "$BLACKLIST" -eq 1 ] && [ -f "$BLACKLIST_CLEAN" ]; then
    CLEAN_ROI="${OUT_DIR}/roi_clean.bed"
    echo "   Subtracting Blacklist regions from targeted ROI BED..."
    bedtools subtract -a "$ROI_BED" -b "$BLACKLIST_CLEAN" > "$CLEAN_ROI"
    ROI_BED="$CLEAN_ROI"
fi

# ==========================================
# 5. TARGET / ANTITARGET GENERATION
# ==========================================

TARGET_BED="${OUT_DIR}/targets.bed"
ANTITARGET_BED="${OUT_DIR}/antitargets.bed"

echo "=========================================="
echo " STAGE 2: Generating targets & antitargets"
echo "=========================================="

if [ "$IS_WGS" -eq 1 ]; then
    echo "   [WGS Mode] Generating targets (antitargets are skipped for WGS)..."
    if [[ "${BIN_SIZE,,}" == "auto" ]]; then
        # Safely grab the exact name of the autobin output
        FIRST_BAM="${BAM_LIST[0]}"
        FIRST_BASE=$(basename "$FIRST_BAM" .bam)
        
        cnvkit.py autobin "${BAM_LIST[@]}" $ACCESS_ARG -m wgs
        mv "${FIRST_BASE}.target.bed" "$TARGET_BED"
    else
        cnvkit.py target "$ACCESS_BED" --split -a "$BIN_SIZE" -o "$TARGET_BED"
    fi
else
    echo "   [vdWGS Hybrid Mode] Generating targets and antitargets..."
    echo "   Creating targets with bin size: $BIN_SIZE bp"
    cnvkit.py target "$ROI_BED" --split -a "$BIN_SIZE" -o "$TARGET_BED"
    
    echo "   Creating antitargets with bin size: $ANTITARGET_BIN_SIZE bp"
    cnvkit.py antitarget "$TARGET_BED" $ACCESS_ARG -a "$ANTITARGET_BIN_SIZE" -o "$ANTITARGET_BED"
fi

# ==========================================
# 6. COVERAGE CALCULATION & SEX SEGREGATION
# ==========================================

echo "=========================================="
echo " STAGE 3: Calculating coverage & segregating by sex"
echo "=========================================="

FEMALE_TARGETS=()
FEMALE_ANTITARGETS=()
MALE_TARGETS=()
MALE_ANTITARGETS=()

for BAM in "${BAM_LIST[@]}"; do
    if [ ! -f "$BAM" ]; then
        echo "ERROR: BAM file not found: $BAM"
        exit 1
    fi
    
    BASENAME=$(basename "$BAM" .bam)
    echo ">>> Processing $BASENAME..."
    
    if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
        echo "   [SAMtools] Index not found. Generating .bai index using $THREADS threads..."
        samtools index -@ "$THREADS" "$BAM"
    fi
    
    T_COV="${OUT_DIR}/${BASENAME}.targetcoverage.cnn"
    echo "   [CNVKit] Calculating target coverage..."
    cnvkit.py coverage "$BAM" "$TARGET_BED" -p "$THREADS" -o "$T_COV"
    
    # Conditionally process antitargets ONLY for vdWGS
    if [ "$IS_WGS" -eq 0 ]; then
        A_COV="${OUT_DIR}/${BASENAME}.antitargetcoverage.cnn"
        echo "   [CNVKit] Calculating antitarget coverage..."
        cnvkit.py coverage "$BAM" "$ANTITARGET_BED" -p "$THREADS" -o "$A_COV"
    fi

    # ----------------------------------------
    # Hardcoded sex segregation
    # ----------------------------------------
    case "$BASENAME" in
        WGS001|WGS002|WGS004|WGS005|WGS006|WGS009|WGS010|WGS015|WGS016|WGS021)
            echo "   [Info] Assigned to FEMALE cohort."
            FEMALE_TARGETS+=("$T_COV")
            [ "$IS_WGS" -eq 0 ] && FEMALE_ANTITARGETS+=("$A_COV")
            ;;
        WGS003|WGS007|WGS008|WGS011|WGS012|WGS013|WGS014|WGS017|WGS018|WGS019|WGS020|WGS022|WGS023)
            echo "   [Info] Assigned to MALE cohort."
            MALE_TARGETS+=("$T_COV")
            [ "$IS_WGS" -eq 0 ] && MALE_ANTITARGETS+=("$A_COV")
            ;;
        *)
            echo "   [WARNING] Unknown sex for $BASENAME. Excluding from final pooled reference."
            ;;
    esac
done

# ==========================================
# 7. SEX-SPECIFIC REFERENCE COMPILATION
# ==========================================

echo "=========================================="
echo " STAGE 4: Building sex-specific pooled references"
echo "=========================================="

OUT_REF_FEMALE="${OUT_REF%.cnn}_female.cnn"
OUT_REF_MALE="${OUT_REF%.cnn}_male.cnn"

# Compile female reference
if [ ${#FEMALE_TARGETS[@]} -gt 0 ]; then
    echo ">>> Compiling Female Reference (${#FEMALE_TARGETS[@]} samples)..."
    if [ "$IS_WGS" -eq 1 ]; then
        cnvkit.py reference "${FEMALE_TARGETS[@]}" \
            -f "$REF_FASTA" -x female -o "$OUT_REF_FEMALE"
    else
        cnvkit.py reference "${FEMALE_TARGETS[@]}" "${FEMALE_ANTITARGETS[@]}" \
            -f "$REF_FASTA" -x female -o "$OUT_REF_FEMALE"
    fi
    echo "   Saved: $OUT_REF_FEMALE"
else
    echo "   [Skip] No female samples survived triage."
fi

# Compile male reference
if [ ${#MALE_TARGETS[@]} -gt 0 ]; then
    echo ">>> Compiling Male Reference (${#MALE_TARGETS[@]} samples)..."
    if [ "$IS_WGS" -eq 1 ]; then
        cnvkit.py reference "${MALE_TARGETS[@]}" \
            -f "$REF_FASTA" -x male -o "$OUT_REF_MALE"
    else
        cnvkit.py reference "${MALE_TARGETS[@]}" "${MALE_ANTITARGETS[@]}" \
            -f "$REF_FASTA" -x male -o "$OUT_REF_MALE"
    fi
    echo "   Saved: $OUT_REF_MALE"
else
    echo "   [Skip] No male samples survived triage."
fi

echo "=========================================="
echo " SUCCESS: Pipeline Execution Complete."
echo "=========================================="

conda deactivate
