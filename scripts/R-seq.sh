#!/bin/bash

# ═══════════════════════════════════════════════════════════════════
#   R-SEQ-Flow  v3  —  Production RNA-seq Pipeline 
# ═══════════════════════════════════════════════════════════════════

# ── Strict mode ──────────────────────────────────────────────────
set -euo pipefail

# ── Colour output ─────────────────────────────────────────────────
GRN="\033[0;32m"; YLW="\033[1;33m"; RED="\033[0;31m"; BLU="\033[0;34m"; RST="\033[0m"
info()  { echo -e "${BLU}[INFO]${RST}  $*"; }
ok()    { echo -e "${GRN}[OK]${RST}    $*"; }
warn()  { echo -e "${YLW}[SKIP]${RST}  $*"; }
err()   { echo -e "${RED}[ERROR]${RST} $*" >&2; }
die()   { err "$*"; exit 1; }

echo ""
echo "═══════════════════════════════════════════════════"
echo "    R-SEQ-Flow v3 — Complete RNA-seq Pipeline      "
echo "═══════════════════════════════════════════════════"
echo ""

# ═══════════════════════════════════════════════════════════════════
# CONFIGURATION  ── edit these to match your machine
# ═══════════════════════════════════════════════════════════════════

THREADS=4       # CPU threads used by fastp / fastqc / kallisto
BOOTSTRAP=100   # Kallisto bootstrap resamples

# Shared folders — live next to this script, reused across experiments
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARED_DIR="/mnt/d/Bioinformatics/RNAseq/SHARED"
SRA_DIR="${SHARED_DIR}/sra"              # prefetch .sra files stored here
FASTQ_DIR="${SHARED_DIR}/fastq"          # converted FASTQs cached here
REF_DIR="${SHARED_DIR}/reference"        # transcriptome FASTA
IDX_DIR="${SHARED_DIR}/kallisto_index"   # kallisto index
GLOBAL_CKPT="${SHARED_DIR}/checkpoints"  # tools / ref / index stamps

# ═══════════════════════════════════════════════════════════════════
# STEP 1 — USER INPUT
# ═══════════════════════════════════════════════════════════════════

echo "────────────────────────────────────────────"
echo "  STEP 1 : Experiment Setup"
echo "────────────────────────────────────────────"
echo ""
read -rp "  Experiment name : " EXPERIMENT_NAME
echo ""
read -rp "  CONTROL SRA ID(s)  (space-separated) : " CONTROL_IDS_INPUT
read -rp "  TREATED SRA ID(s)  (space-separated)  : " TREATED_IDS_INPUT

CONTROL_IDS=($CONTROL_IDS_INPUT)
TREATED_IDS=($TREATED_IDS_INPUT)
ALL_IDS=("${CONTROL_IDS[@]}" "${TREATED_IDS[@]}")

[[ -z "${EXPERIMENT_NAME}" ]]   && die "Experiment name cannot be empty."
[[ ${#CONTROL_IDS[@]} -eq 0 ]] && die "No CONTROL IDs entered."
[[ ${#TREATED_IDS[@]} -eq 0 ]] && die "No TREATED IDs entered."

echo ""
echo "  Experiment  : ${EXPERIMENT_NAME}"
echo "  Control IDs : ${CONTROL_IDS[*]}"
echo "  Treated IDs : ${TREATED_IDS[*]}"
echo ""

# ═══════════════════════════════════════════════════════════════════
# STEP 2 — CREATE DIRECTORY STRUCTURE (only if new)
# ═══════════════════════════════════════════════════════════════════

echo "────────────────────────────────────────────"
echo "  STEP 2 : Directory Setup"
echo "────────────────────────────────────────────"

PROJECT_DIR="${SCRIPT_DIR}/${EXPERIMENT_NAME}"
EXP_CKPT="${PROJECT_DIR}/checkpoints"

# Check if this is a RESUMED run
RESUME_MODE=0
if [[ -d "${PROJECT_DIR}" && -d "${EXP_CKPT}" ]]; then
    RESUME_MODE=1
    info "RESUME MODE DETECTED: Using existing experiment folder"
fi

# Create only missing directories (preserves existing data)
mkdir -p "${SRA_DIR}" "${FASTQ_DIR}" "${REF_DIR}" "${IDX_DIR}" "${GLOBAL_CKPT}"
mkdir -p "${PROJECT_DIR}"
mkdir -p "${EXP_CKPT}"

# Only create experiment subdirectories if new (don't overwrite on resume)
if [[ ${RESUME_MODE} -eq 0 ]]; then
    mkdir -p "${PROJECT_DIR}"/{trimmed_data,qc_results,kallisto_output,metadata,logs}
    info "NEW EXPERIMENT: Created all directories"
else
    # On resume, ensure log directory exists but don't touch others
    mkdir -p "${PROJECT_DIR}/logs"
    warn "Experiment directory already exists — preserving all data"
fi

echo "  Shared cache   : ${SHARED_DIR}"
echo "  Experiment dir : ${PROJECT_DIR}"
[[ ${RESUME_MODE} -eq 1 ]] && echo "  Mode           : RESUME (skipping completed steps)"
echo ""

# ═══════════════════════════════════════════════════════════════════
# CHECKPOINT HELPERS
# global_done/mark → SHARED/checkpoints/  (tools, ref, index)
# exp_done/mark    → <experiment>/checkpoints/  (per-sample steps)
# ═══════════════════════════════════════════════════════════════════

global_done() { [[ -f "${GLOBAL_CKPT}/${1}.done" ]]; }
global_mark() { touch "${GLOBAL_CKPT}/${1}.done";    }
exp_done()    { [[ -f "${EXP_CKPT}/${1}.done" ]];    }
exp_mark()    { touch "${EXP_CKPT}/${1}.done";        }

# ═══════════════════════════════════════════════════════════════════
# STEP 3 — INSTALL TOOLS  (once per machine)
# ═══════════════════════════════════════════════════════════════════

echo "────────────────────────────────────────────"
echo "  STEP 3 : Tool Installation"
echo "────────────────────────────────────────────"

if ! global_done "tools_installed"; then
    info "Installing tools (requires sudo)..."
    sudo apt-get update -y -qq
    sudo apt-get install -y -qq \
        fastqc \
        fastp \
        kallisto \
        sra-toolkit \
        wget \
        pigz \
        unzip
    global_mark "tools_installed"
    ok "All tools installed."
else
    warn "Tools already installed — skipping."
fi
echo ""

# ═══════════════════════════════════════════════════════════════════
# STEP 4 — DOWNLOAD + CONVERT  (smart 3-level cache)
#
#  Per sample, checked in order:
#
#  LEVEL 1  FASTQ files already in SHARED/fastq/
#           → skip entirely
#
#  LEVEL 2  .sra file already in SHARED/sra/<ID>/
#           → skip prefetch, run fasterq-dump on local file only
#           (local conversion = very fast, just CPU, no network)
#
#  LEVEL 3  Nothing exists
#           → prefetch (resumable) then fasterq-dump locally
#
#  After fasterq-dump succeeds, .sra file is deleted to save disk.
# ═══════════════════════════════════════════════════════════════════

echo "────────────────────────────────────────────"
echo "  STEP 4 : SRA Download & FASTQ Conversion"
echo "────────────────────────────────────────────"
echo ""

AVAIL_KB=$(df -k "${SHARED_DIR}" | awk 'NR==2{print $4}')
AVAIL_GB=$(( AVAIL_KB / 1024 / 1024 ))
if [[ ${AVAIL_GB} -lt 10 ]]; then
    warn "Only ${AVAIL_GB} GB free — pipeline may run out of disk space."
else
    info "Disk space available: ${AVAIL_GB} GB"
fi
echo ""

for ID in "${ALL_IDS[@]}"; do

    SRA_FILE="${SRA_DIR}/${ID}/${ID}.sra"

    # LEVEL 1 — FASTQs already cached
    if [[ -f "${FASTQ_DIR}/${ID}_1.fastq" && -f "${FASTQ_DIR}/${ID}_2.fastq" ]]; then
        ok "[CACHE HIT - PE]  ${ID} → FASTQs already present."
        continue
    fi
    if [[ -f "${FASTQ_DIR}/${ID}.fastq" ]]; then
        ok "[CACHE HIT - SE]  ${ID} → FASTQ already present."
        continue
    fi

    # LEVEL 2 — .sra already cached, skip prefetch
    if [[ ! -f "${SRA_FILE}" ]]; then
        # LEVEL 3 — download fresh
        info "[DOWNLOAD]  ${ID} — running prefetch (resumable)..."
        prefetch "${ID}" \
            --output-directory "${SRA_DIR}" \
            --progress \
            --resume yes \
            2>> "${PROJECT_DIR}/logs/${ID}_prefetch.log" \
        || die "${ID} prefetch failed. Check logs/${ID}_prefetch.log"
        ok "${ID} prefetch complete."
    else
        info "[CACHE HIT - SRA]  ${ID}.sra found — skipping prefetch."
    fi

    # Convert local .sra → FASTQ (fast, no network needed)
    info "  Converting ${ID} → FASTQ  (${THREADS} threads)..."
    fasterq-dump "${SRA_FILE}" \
        --split-files \
        --threads "${THREADS}" \
        --progress \
        --outdir "${FASTQ_DIR}" \
        2>> "${PROJECT_DIR}/logs/${ID}_fasterq.log" \
    || die "${ID} fasterq-dump failed. Check logs/${ID}_fasterq.log"

    # Remove .sra to free disk space
    rm -rf "${SRA_DIR:?}/${ID}"
    ok "${ID} converted successfully. .sra removed."

done

echo ""
info "FASTQ cache:"
ls -lh "${FASTQ_DIR}/"
echo ""

# ═══════════════════════════════════════════════════════════════════
# PE / SE DETECTION  (checks FASTQ cache after real conversion)
# ═══════════════════════════════════════════════════════════════════

check_reads() {
    local SAMPLE=$1
    if [[ -f "${FASTQ_DIR}/${SAMPLE}_1.fastq" && \
          -f "${FASTQ_DIR}/${SAMPLE}_2.fastq" ]]; then
        echo "PE"
    elif [[ -f "${FASTQ_DIR}/${SAMPLE}.fastq" ]]; then
        echo "SE"
    else
        echo "UNKNOWN"
    fi
}

echo "  Detected library types:"
for ID in "${ALL_IDS[@]}"; do
    TYPE=$(check_reads "${ID}")
    echo "    ${ID}  →  ${TYPE}"
    [[ "${TYPE}" == "UNKNOWN" ]] && die "Cannot detect read type for ${ID}. Files missing."
done
echo ""

# ═══════════════════════════════════════════════════════════════════
# STEP 5 — FASTQC ON RAW DATA  [CHECKPOINT-PROTECTED]
# ═══════════════════════════════════════════════════════════════════

echo "────────────────────────────────────────────"
echo "  STEP 5 : FastQC (Raw)"
echo "────────────────────────────────────────────"

if ! exp_done "fastqc_raw"; then
    info "Running FastQC on raw FASTQs..."

    RAW_FASTQS=()
    for ID in "${ALL_IDS[@]}"; do
        TYPE=$(check_reads "${ID}")
        if [[ "${TYPE}" == "PE" ]]; then
            RAW_FASTQS+=("${FASTQ_DIR}/${ID}_1.fastq" "${FASTQ_DIR}/${ID}_2.fastq")
        else
            RAW_FASTQS+=("${FASTQ_DIR}/${ID}.fastq")
        fi
    done

    fastqc "${RAW_FASTQS[@]}" \
        -o "${PROJECT_DIR}/qc_results/" \
        -t "${THREADS}" \
        2>> "${PROJECT_DIR}/logs/fastqc_raw.log" \
    || die "FastQC raw failed. Check logs/fastqc_raw.log"

    exp_mark "fastqc_raw"
    ok "FastQC raw complete."
else
    warn "FastQC raw already done — skipping."
fi
echo ""

# ═══════════════════════════════════════════════════════════════════
# STEP 6 — TRIMMING  (fastp, per-sample checkpoint)
# ═══════════════════════════════════════════════════════════════════

echo "────────────────────────────────────────────"
echo "  STEP 6 : Trimming (fastp)"
echo "────────────────────────────────────────────"
echo ""

trim_sample() {
    local SAMPLE=$1
    local CKPT="trim_${SAMPLE}"

    if exp_done "${CKPT}"; then
        warn "${SAMPLE} already trimmed — skipping."
        return 0
    fi

    local TYPE
    TYPE=$(check_reads "${SAMPLE}")
    info "Trimming ${SAMPLE}  [${TYPE}]..."

    if [[ "${TYPE}" == "PE" ]]; then
        fastp \
            -i  "${FASTQ_DIR}/${SAMPLE}_1.fastq" \
            -I  "${FASTQ_DIR}/${SAMPLE}_2.fastq" \
            -o  "${PROJECT_DIR}/trimmed_data/${SAMPLE}_1.trimmed.fastq" \
            -O  "${PROJECT_DIR}/trimmed_data/${SAMPLE}_2.trimmed.fastq" \
            -h  "${PROJECT_DIR}/trimmed_data/${SAMPLE}_fastp.html" \
            -j  "${PROJECT_DIR}/trimmed_data/${SAMPLE}_fastp.json" \
            -w  "${THREADS}" \
            --detect_adapter_for_pe \
            2>> "${PROJECT_DIR}/logs/${SAMPLE}_fastp.log" \
        || die "fastp failed for ${SAMPLE}. Check logs/${SAMPLE}_fastp.log"

    elif [[ "${TYPE}" == "SE" ]]; then
        fastp \
            -i  "${FASTQ_DIR}/${SAMPLE}.fastq" \
            -o  "${PROJECT_DIR}/trimmed_data/${SAMPLE}.trimmed.fastq" \
            -h  "${PROJECT_DIR}/trimmed_data/${SAMPLE}_fastp.html" \
            -j  "${PROJECT_DIR}/trimmed_data/${SAMPLE}_fastp.json" \
            -w  "${THREADS}" \
            2>> "${PROJECT_DIR}/logs/${SAMPLE}_fastp.log" \
        || die "fastp failed for ${SAMPLE}. Check logs/${SAMPLE}_fastp.log"
    fi

    exp_mark "${CKPT}"
    ok "${SAMPLE} trimmed."
}

for ID in "${CONTROL_IDS[@]}"; do trim_sample "${ID}"; done
for ID in "${TREATED_IDS[@]}";  do trim_sample "${ID}"; done
echo ""

# ═══════════════════════════════════════════════════════════════════
# STEP 7 — FASTQC AFTER TRIMMING  [CHECKPOINT-PROTECTED]
# ═══════════════════════════════════════════════════════════════════

echo "────────────────────────────────────────────"
echo "  STEP 7 : FastQC (Trimmed)"
echo "────────────────────────────────────────────"

if ! exp_done "fastqc_trimmed"; then
    info "Running FastQC on trimmed FASTQs..."

    TRIMMED_FASTQS=()
    for ID in "${ALL_IDS[@]}"; do
        TYPE=$(check_reads "${ID}")
        if [[ "${TYPE}" == "PE" ]]; then
            TRIMMED_FASTQS+=(
                "${PROJECT_DIR}/trimmed_data/${ID}_1.trimmed.fastq"
                "${PROJECT_DIR}/trimmed_data/${ID}_2.trimmed.fastq"
            )
        else
            TRIMMED_FASTQS+=("${PROJECT_DIR}/trimmed_data/${ID}.trimmed.fastq")
        fi
    done

    fastqc "${TRIMMED_FASTQS[@]}" \
        -o "${PROJECT_DIR}/qc_results/" \
        -t "${THREADS}" \
        2>> "${PROJECT_DIR}/logs/fastqc_trimmed.log" \
    || die "FastQC trimmed failed. Check logs/fastqc_trimmed.log"

    exp_mark "fastqc_trimmed"
    ok "FastQC trimmed complete."
else
    warn "FastQC trimmed already done — skipping."
fi
echo ""

# ═══════════════════════════════════════════════════════════════════
# STEP 8 — DOWNLOAD TRANSCRIPTOME  (once ever, shared)
# ═══════════════════════════════════════════════════════════════════

echo "────────────────────────────────────────────"
echo "  STEP 8 : Reference Transcriptome"
echo "────────────────────────────────────────────"

REF_FA="${REF_DIR}/gencode.v29.transcripts.fa"
REF_GZ="${REF_FA}.gz"

if ! global_done "reference_downloaded"; then
    info "Downloading GENCODE v29 transcriptome (~750 MB)..."
    wget -c \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.transcripts.fa.gz" \
        -O "${REF_GZ}" \
        2>> "${PROJECT_DIR}/logs/wget_transcriptome.log" \
    || die "Transcriptome download failed. Check logs/wget_transcriptome.log"

    info "Decompressing..."
    gunzip -f "${REF_GZ}"
    global_mark "reference_downloaded"
    ok "Transcriptome ready."
else
    warn "Transcriptome already downloaded (shared) — skipping."
fi
echo ""

# ═══════════════════════════════════════════════════════════════════
# STEP 9 — KALLISTO INDEX  (once ever, shared)
# ═══════════════════════════════════════════════════════════════════

echo "────────────────────────────────────────────"
echo "  STEP 9 : Kallisto Index"
echo "────────────────────────────────────────────"

IDX="${IDX_DIR}/gencode_v29.idx"

if ! global_done "kallisto_index_built"; then
    info "Building Kallisto index (done only once, ~5 min)..."
    [[ ! -f "${REF_FA}" ]] && die "Transcriptome FASTA missing: ${REF_FA}"

    kallisto index \
        -i "${IDX}" \
        "${REF_FA}" \
        2>> "${PROJECT_DIR}/logs/kallisto_index.log" \
    || die "Kallisto index failed. Check logs/kallisto_index.log"

    global_mark "kallisto_index_built"
    ok "Kallisto index built."
else
    warn "Kallisto index already exists (shared) — skipping."
fi
echo ""

# ═══════════════════════════════════════════════════════════════════
# STEP 10 — KALLISTO QUANTIFICATION  (per-sample checkpoint)
# ═══════════════════════════════════════════════════════════════════

echo "────────────────────────────────────────────"
echo "  STEP 10 : Kallisto Quantification"
echo "────────────────────────────────────────────"
echo ""

quant_sample() {
    local SAMPLE=$1
    local CKPT="quant_${SAMPLE}"

    if exp_done "${CKPT}"; then
        warn "${SAMPLE} already quantified — skipping."
        return 0
    fi

    local TYPE
    TYPE=$(check_reads "${SAMPLE}")
    info "Quantifying ${SAMPLE}  [${TYPE}]  -b ${BOOTSTRAP}  -t ${THREADS}..."

    mkdir -p "${PROJECT_DIR}/kallisto_output/${SAMPLE}"

    if [[ "${TYPE}" == "PE" ]]; then
        kallisto quant \
            -i "${IDX}" \
            -o "${PROJECT_DIR}/kallisto_output/${SAMPLE}" \
            -b "${BOOTSTRAP}" \
            -t "${THREADS}" \
            "${PROJECT_DIR}/trimmed_data/${SAMPLE}_1.trimmed.fastq" \
            "${PROJECT_DIR}/trimmed_data/${SAMPLE}_2.trimmed.fastq" \
            2>> "${PROJECT_DIR}/logs/${SAMPLE}_kallisto.log" \
        || die "Kallisto quant failed for ${SAMPLE}. Check logs/${SAMPLE}_kallisto.log"

    elif [[ "${TYPE}" == "SE" ]]; then
        kallisto quant \
            -i "${IDX}" \
            -o "${PROJECT_DIR}/kallisto_output/${SAMPLE}" \
            -b "${BOOTSTRAP}" \
            -t "${THREADS}" \
            --single -l 200 -s 20 \
            "${PROJECT_DIR}/trimmed_data/${SAMPLE}.trimmed.fastq" \
            2>> "${PROJECT_DIR}/logs/${SAMPLE}_kallisto.log" \
        || die "Kallisto quant failed for ${SAMPLE}. Check logs/${SAMPLE}_kallisto.log"
    fi

    exp_mark "${CKPT}"
    ok "${SAMPLE} quantified → kallisto_output/${SAMPLE}/abundance.tsv"
}

for ID in "${CONTROL_IDS[@]}"; do quant_sample "${ID}"; done
for ID in "${TREATED_IDS[@]}";  do quant_sample "${ID}"; done
echo ""

# ═══════════════════════════════════════════════════════════════════
# STEP 11 — DESEQ2 METADATA FILE  [CHECKPOINT-PROTECTED]
# ═══════════════════════════════════════════════════════════════════

echo "────────────────────────────────────────────"
echo "  STEP 11 : DESeq2 Metadata"
echo "────────────────────────────────────────────"

METADATA="${PROJECT_DIR}/metadata/samples.txt"

if ! exp_done "metadata_written"; then
    mkdir -p "${PROJECT_DIR}/metadata"
    echo -e "sample\tcondition\tpath" > "${METADATA}"

    for ID in "${CONTROL_IDS[@]}"; do
        echo -e "${ID}\tcontrol\t${PROJECT_DIR}/kallisto_output/${ID}/abundance.tsv" >> "${METADATA}"
    done
    for ID in "${TREATED_IDS[@]}"; do
        echo -e "${ID}\ttreated\t${PROJECT_DIR}/kallisto_output/${ID}/abundance.tsv" >> "${METADATA}"
    done

    exp_mark "metadata_written"
    ok "Metadata written."
else
    warn "Metadata already written — skipping."
fi

echo ""
echo "Metadata file:"
cat "${METADATA}"
echo ""

exp_mark "EXPERIMENT_COMPLETE"

echo "═══════════════════════════════════════════════════"
echo "    PIPELINE COMPLETED SUCCESSFULLY                "
echo "═══════════════════════════════════════════════════"
echo ""
echo "  Experiment     : ${EXPERIMENT_NAME}"
echo "  Project dir    : ${PROJECT_DIR}"
echo "  QC reports     : ${PROJECT_DIR}/qc_results/"
echo "  Kallisto out   : ${PROJECT_DIR}/kallisto_output/"
echo "  DESeq2 meta    : ${METADATA}"
echo "  Logs           : ${PROJECT_DIR}/logs/"
echo ""
echo "  Shared (reused across all future experiments):"
echo "    FASTQ cache  : ${FASTQ_DIR}/"
echo "    Reference    : ${REF_FA}"
echo "    Index        : ${IDX}"
echo ""
echo "  To resume after a crash:"
echo "    → Re-run this script with the SAME experiment name + IDs"
echo "    → Completed steps detected via checkpoints & skipped"
echo "════════════════════════════════════════════════════════"
echo ""
echo ""
read -p "Continue to DESeq2 analysis? (yes/no): " RUN_DESEQ2

if [[ "$RUN_DESEQ2" =~ ^([Yy][Ee][Ss]|[Yy])$ ]]; then

    echo ""
    echo "Starting DESeq2 analysis..."

    Rscript deseq2_analysis.R \
        "${PROJECT_DIR}/metadata/samples.txt" \
        "${PROJECT_DIR}/deseq2_results"

else

    echo ""
    echo "Skipping DESeq2 analysis."
    echo "Kallisto results available in:"
    echo "${PROJECT_DIR}/kallisto_output"

fi
