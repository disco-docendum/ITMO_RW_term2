#!/bin/bash

# Stop immediately if any command fails
set -e

# ==========================================
# USAGE & HELP
# ==========================================

usage() {
    echo "Usage: $0 -i <input_bam> -r <ref_genome> -c <chromosome> [-f <ref_fasta>] [-g <access_bed>] [-p <pooled_ref.cnn>] [-s <bin_size>|auto] [-m <method>] [-b <blacklist_flag>] [-@ <threads>]"
    echo ""
    echo "Arguments:"
    echo "  -i  Input BAM file. Use \"*\" (quoted) to process all .bam files in the current directory."
    echo "  -r  Reference genome build (hs37d5, hg19, hg38). REQUIRED for auto-downloads & blacklists."
    echo "  -c  Chromosome to process (e.g., \"20\" or \"chr20\"). Use \"all\" for all chromosomes (WGS)."
    echo "  -f  Reference FASTA file. Optional: If omitted, the script will auto-download based on -r."
    echo "  -g  Access BED file. Optional: Generated on-the-fly from FASTA if omitted."
    echo "  -p  Path to a pooled reference (.cnn) file. If omitted, pipeline builds a Flat Reference."
    echo "  -s  Bin size (default: auto). Use 'auto' for autobinning or a specific integer (e.g., 50000)."
    echo "  -m  Method of segmentation (default is cbs, supports hmm, hmm-tumor, hmm-germline)."
    echo "  -t  Threshold of statistical significance (default: 0.00001)."
    echo "  -b  Blacklist exclusion: 1 for On, 0 for Off (default: 1)."
    echo "  -@  Number of threads (default: 4)."
    echo "  -h  Show this help message."
    exit 1
}

# Default values
BIN_SIZE="auto"
METHOD="cbs"
THRESH=0.00001
BLACKLIST=1
THREADS=4
POOLED_REF=""
REF_GENOME=""
REF_FASTA=""
ACCESS_BED=""
CHROMOSOME="all"

while getopts "i:f:r:g:p:c:s:m:t:b:@:h" opt; do
    case $opt in
        i) INPUT_BAM="$OPTARG" ;;
        f) REF_FASTA="$OPTARG" ;;
        r) REF_GENOME="$OPTARG" ;;
        g) ACCESS_BED="$OPTARG" ;;
        p) POOLED_REF="$OPTARG" ;;
        c) CHROMOSOME="$OPTARG" ;;
        s) BIN_SIZE="$OPTARG" ;;
        m) METHOD="$OPTARG" ;;
        t) THRESH="$OPTARG" ;;
        b) BLACKLIST="$OPTARG" ;;
        @) THREADS="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [ -z "$INPUT_BAM" ] || [ -z "$REF_GENOME" ]; then
    echo "Error: Missing required arguments. -i and -r are mandatory."
    usage
fi

OUT_CNVKIT="./wgs_cnvkit_out"
mkdir -p "$OUT_CNVKIT"

# ==========================================
# 1. INPUT RESOLUTION
# ==========================================
echo "=========================================="
echo " PHASE 1: Input Resolution"
echo "=========================================="

BAM_FILES=()
for file in $INPUT_BAM; do
    if [ -f "$file" ] && [[ "$file" == *.bam ]]; then
        BAM_FILES+=("$file")
    fi
done

if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "Error: No valid BAM files found matching the input."
    exit 1
fi
echo "Found ${#BAM_FILES[@]} BAM files for processing."

# ==========================================
# 2. BLACKLIST SETUP
# ==========================================
echo "=========================================="
echo " PHASE 2: Blacklist Setup"
echo "=========================================="

if [ "$BLACKLIST" -eq 1 ]; then
    BLACKLIST_DIR="./blacklist_${REF_GENOME}"
    mkdir -p "$BLACKLIST_DIR"
    
    if [ "$REF_GENOME" == "hg19" ]; then
        BLACKLIST_URL="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz"
        BLACKLIST_RAW="${BLACKLIST_DIR}/hg19-blacklist.v2.bed.gz"
        BLACKLIST_CLEAN="${BLACKLIST_DIR}/hg19-blacklist_clean.bed"
    elif [ "$REF_GENOME" == "hg38" ]; then
        BLACKLIST_URL="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"
        BLACKLIST_RAW="${BLACKLIST_DIR}/hg38-blacklist.v2.bed.gz"
        BLACKLIST_CLEAN="${BLACKLIST_DIR}/hg38-blacklist_clean.bed"
    else
        echo "Warning: No standard blacklist available for $REF_GENOME. Skipping blacklist."
        BLACKLIST=0
    fi

    if [ "$BLACKLIST" -eq 1 ]; then
        if [ ! -f "$BLACKLIST_CLEAN" ]; then
            echo ">>> Downloading Boyle-Lab Blacklist for $REF_GENOME..."
            wget -q -O "$BLACKLIST_RAW" "$BLACKLIST_URL"
            zcat "$BLACKLIST_RAW" | sed 's/^chr//' | sort -k1,1 -k2,2n > "$BLACKLIST_CLEAN"
            echo "   Saved clean blacklist to $BLACKLIST_CLEAN"
        else
            echo ">>> Blacklist already exists at $BLACKLIST_CLEAN"
        fi
    fi
else
    echo ">>> Blacklist exclusion disabled."
fi

# ==========================================
# 2.5 FASTA AUTO-DOWNLOADER
# ==========================================
echo "=========================================="
echo " PHASE 2.5: FASTA Reference Setup"
echo "=========================================="

REF_DIR="./reference_genomes"
mkdir -p "$REF_DIR"

if [ -z "$REF_FASTA" ]; then
    echo ">>> No FASTA provided. Checking auto-download configuration for $REF_GENOME..."
    if [ "$REF_GENOME" == "hs37d5" ]; then
        REF_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
        REF_FASTA="${REF_DIR}/hs37d5.fa"
    elif [ "$REF_GENOME" == "hg19" ]; then
        REF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
        REF_FASTA="${REF_DIR}/hg19.fa"
    elif [ "$REF_GENOME" == "hg38" ]; then
        REF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
        REF_FASTA="${REF_DIR}/hg38.fa"
    else
        echo "Error: Auto-download not supported for $REF_GENOME. Please provide -f <fasta_file> manually."
        exit 1
    fi

    if [ ! -f "$REF_FASTA" ]; then
        if [[ "$REF_URL" == *.gz ]]; then
            REF_RAW="${REF_FASTA}.gz"
        else
            REF_RAW="${REF_FASTA}.raw"
        fi

        if [ ! -f "$REF_RAW" ]; then
            echo "   Downloading $REF_GENOME Raw Reference (This may take a while)..."
            wget -c -O "$REF_RAW" "$REF_URL"
        fi
        
        echo "   Processing reference sequence..."
        if [[ "$REF_RAW" == *.gz ]]; then
            zcat "$REF_RAW" > "$REF_FASTA"
        else
            cp "$REF_RAW" "$REF_FASTA"
        fi
    else
        echo "   [FASTA Configured] Reference already exists at $REF_FASTA"
    fi
else
    echo ">>> Using user-supplied FASTA: $REF_FASTA"
fi


# ==========================================
# 3. REFERENCE & TARGET PREPARATION
# ==========================================
echo "=========================================="
echo " PHASE 3: Target Preparation"
echo "=========================================="

EMPTY_ANTI="${OUT_CNVKIT}/empty_antitarget.bed"
touch "$EMPTY_ANTI"

if [ -z "$ACCESS_BED" ]; then
    echo ">>> No Access BED provided. Generating from FASTA (this may take a few minutes)..."
    ACCESS_BED="${OUT_CNVKIT}/access_${REF_GENOME}.bed"
    if [ ! -f "$ACCESS_BED" ]; then
        cnvkit.py access "$REF_FASTA" -o "$ACCESS_BED"
    fi
fi

if [ "$BIN_SIZE" == "auto" ]; then
    echo ">>> Running CNVkit autobin to determine optimal bin size..."
    FIRST_BAM="${BAM_FILES[0]}"
    cnvkit.py autobin "$FIRST_BAM" -m wgs -g "$ACCESS_BED"
    TARGET_BED="$(basename "${FIRST_BAM%.bam}.target.bed")"
    echo "   Autobin created target BED: $TARGET_BED"
else
    echo ">>> Using fixed bin size: $BIN_SIZE"
    TARGET_BED="${OUT_CNVKIT}/targets_${BIN_SIZE}.bed"
    if [ ! -f "$TARGET_BED" ]; then
        echo ">>> Generating uniform WGS targets using fixed bin size..."
        cnvkit.py target "$ACCESS_BED" --split -b "$BIN_SIZE" -o "$TARGET_BED"
    fi
fi

if [ -z "$POOLED_REF" ]; then
    echo ">>> No pooled reference provided. Generating Flat Reference..."
    FLAT_REF="${OUT_CNVKIT}/flat_reference_all.cnn"
    if [ ! -f "$FLAT_REF" ]; then
        cnvkit.py reference -f "$REF_FASTA" -o "$FLAT_REF" -t "$TARGET_BED" -a "$EMPTY_ANTI"
    fi
else
    echo ">>> User supplied pooled reference baseline mapping: $POOLED_REF"
fi


# ==========================================
# 4. BATCH PROCESSING
# ==========================================

for BAM in "${BAM_FILES[@]}"; do
    BASENAME=$(basename "$BAM" .bam)
    SAMPLE_OUT_DIR="${OUT_CNVKIT}/${BASENAME}"
    mkdir -p "$SAMPLE_OUT_DIR"

    echo "=========================================="
    echo " PHASE 4: Processing $BASENAME"
    echo "=========================================="
    
    # -----------------------------------------------------
    # DYNAMIC SEX ASSIGNMENT & REFERENCE ROUTING
    # -----------------------------------------------------
    case "$BASENAME" in
        WGS001|WGS002|WGS004|WGS005|WGS006|WGS009|WGS010|WGS015|WGS016|WGS021)
            SAMPLE_SEX="female"
            ;;
        WGS003|WGS007|WGS008|WGS011|WGS012|WGS013|WGS014|WGS017|WGS018|WGS019|WGS020|WGS022|WGS023)
            SAMPLE_SEX="male"
            ;;
        *)
            echo "   [WARNING] Unknown sex for $BASENAME. Assuming female."
            SAMPLE_SEX="female"
            ;;
    esac

    if [ -n "$POOLED_REF" ]; then
        if [ "$SAMPLE_SEX" == "female" ] && [ -f "${POOLED_REF%.cnn}_female.cnn" ]; then
            ACTIVE_REF="${POOLED_REF%.cnn}_female.cnn"
        elif [ "$SAMPLE_SEX" == "male" ] && [ -f "${POOLED_REF%.cnn}_male.cnn" ]; then
            ACTIVE_REF="${POOLED_REF%.cnn}_male.cnn"
        else
            ACTIVE_REF="$POOLED_REF" 
        fi
    else
        ACTIVE_REF="$FLAT_REF"
    fi

    if [ ! -f "$ACTIVE_REF" ]; then
        echo "ERROR: Reference not found at $ACTIVE_REF"
        exit 1
    fi

    # -----------------------------------------------------
    # CNVKIT CALLING (With HMM/CBS Agnostic Flags)
    # -----------------------------------------------------
    echo "   [CNVkit] Calculating Target coverage..."
    cnvkit.py coverage "$BAM" "$TARGET_BED" -p "$THREADS" -o "$SAMPLE_OUT_DIR/${BASENAME}.targetcoverage.cnn"

    echo "   [CNVkit] Calculating (Empty) Antitarget coverage..."
    cnvkit.py coverage "$BAM" "$EMPTY_ANTI" -p "$THREADS" -o "$SAMPLE_OUT_DIR/${BASENAME}.antitargetcoverage.cnn"

    echo "   [CNVkit] Fixing log2 ratios against $SAMPLE_SEX reference..."
    cnvkit.py fix "$SAMPLE_OUT_DIR/${BASENAME}.targetcoverage.cnn" \
                  "$SAMPLE_OUT_DIR/${BASENAME}.antitargetcoverage.cnn" \
                  "$ACTIVE_REF" \
                  -x "$SAMPLE_SEX" \
                  -o "$SAMPLE_OUT_DIR/${BASENAME}.cnr"

    echo "   [CNVkit] Segmenting using method: $METHOD ..."
    
    if [[ "$METHOD" == "cbs" ]]; then
        SEG_ARGS="--smooth-cbs --drop-outliers 10 -t $THRESH"
    elif [[ "$METHOD" == hmm* ]]; then
        SEG_ARGS=""
    else
        SEG_ARGS="--drop-outliers 10" # Fallback for haar, etc.
    fi

    # Conditionally add the chromosome flag if not set to "all"
    if [ "$CHROMOSOME" != "all" ]; then
        SEG_ARGS="$SEG_ARGS -c $CHROMOSOME"
        echo "   [CNVkit] Restricting segmentation to chromosome: $CHROMOSOME"
    fi

    cnvkit.py segment "$SAMPLE_OUT_DIR/${BASENAME}.cnr" \
        -m "$METHOD" \
        -p "$THREADS" \
        $SEG_ARGS \
        -o "$SAMPLE_OUT_DIR/${BASENAME}_${METHOD}.cns"

    echo "   [CNVkit] Calling Copy Number..."
    cnvkit.py call "$SAMPLE_OUT_DIR/${BASENAME}_${METHOD}.cns" \
        -x "$SAMPLE_SEX" \
        --filter cn \
        --min-variant-bins 3 \
        -o "$SAMPLE_OUT_DIR/${BASENAME}_${METHOD}.call.cns"

    # -----------------------------------------------------
    # EXPORT & BLACKLIST
    # -----------------------------------------------------
    echo "   [CNVkit] Exporting to BED..."
    tail -n +2 "$SAMPLE_OUT_DIR/${BASENAME}_${METHOD}.call.cns" | \
    awk -v id="$BASENAME" '{print $1 "\t" $2 "\t" $3 "\t" id "\t" $9}' \
    > "${SAMPLE_OUT_DIR}/${BASENAME}_${METHOD}_cnvkit_results_raw.bed"

    if [ "$BLACKLIST" -eq 1 ] && [ -f "$BLACKLIST_CLEAN" ]; then
        echo "   [bedtools] Removing blacklist overlaps from final BED..."
        bedtools subtract -a "${SAMPLE_OUT_DIR}/${BASENAME}_${METHOD}_cnvkit_results_raw.bed" \
                          -b "$BLACKLIST_CLEAN" \
                          > "${SAMPLE_OUT_DIR}/${BASENAME}_${METHOD}_cnvkit_results.bed"
        rm "${SAMPLE_OUT_DIR}/${BASENAME}_${METHOD}_cnvkit_results_raw.bed"
    else
        mv "${SAMPLE_OUT_DIR}/${BASENAME}_${METHOD}_cnvkit_results_raw.bed" "${SAMPLE_OUT_DIR}/${BASENAME}_${METHOD}_cnvkit_results.bed"
    fi

    echo "Finished $BASENAME"
done

echo "=========================================="
echo " Pipeline Complete."
echo "=========================================="
