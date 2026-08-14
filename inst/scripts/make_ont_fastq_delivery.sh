#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/make_ont_fastq_delivery.sh \
    --input /path/to/demux_fastq \
    --output /path/to/delivery_root \
    --project PROJECT_ID \
    [--sample-sheet /path/to/sample_sheet.xlsx|csv|tsv] \
    [--threads 8] \
    [--overwrite] \
    [--no-reuse-nanoplot] \
    [--skip-nanoplot] \
    [--skip-multiqc]

Description:
  Build a customer-deliverable package for demultiplexed ONT FASTQ files.
  The script does not filter, rename, or modify FASTQ contents.

Outputs:
  <output>/<project>_delivery/
    00_sample_sheet/
    01_fastq/
    02_qc_report/
      fastq_stats.tsv
      nanoplot/
      multiqc_report.html
      multiqc_data/
    03_md5/
      md5.txt
    README.txt
  <output>/<project>_delivery.tar.gz

Required:
  seqkit

Recommended:
  NanoPlot
  multiqc
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

has_cmd() {
  command -v "$1" >/dev/null 2>&1
}

md5_write() {
  local target_dir="$1"
  local md5_file="$2"

  if has_cmd md5sum; then
    (
      cd "$target_dir"
      find . -type f | sort | while IFS= read -r file; do
        md5sum "${file#./}"
      done
    ) > "$md5_file"
  elif has_cmd md5; then
    (
      cd "$target_dir"
      find . -type f | sort | while IFS= read -r file; do
        local rel="${file#./}"
        local sum
        sum="$(md5 -q "$rel")"
        printf '%s  %s\n' "$sum" "$rel"
      done
    ) > "$md5_file"
  else
    die "Neither md5sum nor md5 was found."
  fi
}

INPUT_DIR=""
OUTPUT_ROOT=""
PROJECT_ID=""
SAMPLE_SHEET=""
THREADS=8
RUN_NANOPLOT=1
RUN_MULTIQC=1
OVERWRITE=0
REUSE_NANOPLOT=1
FASTQ_LIST=""
PRESERVED_NANOPLOT_DIR=""

cleanup() {
  [[ -n "${FASTQ_LIST:-}" ]] && rm -f "$FASTQ_LIST"
  [[ -n "${PRESERVED_NANOPLOT_DIR:-}" ]] && rm -rf "$PRESERVED_NANOPLOT_DIR"
}
trap cleanup EXIT

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)
      INPUT_DIR="${2:-}"
      shift 2
      ;;
    --output)
      OUTPUT_ROOT="${2:-}"
      shift 2
      ;;
    --project)
      PROJECT_ID="${2:-}"
      shift 2
      ;;
    --sample-sheet)
      SAMPLE_SHEET="${2:-}"
      shift 2
      ;;
    --threads)
      THREADS="${2:-}"
      shift 2
      ;;
    --overwrite)
      OVERWRITE=1
      shift
      ;;
    --no-reuse-nanoplot)
      REUSE_NANOPLOT=0
      shift
      ;;
    --skip-nanoplot)
      RUN_NANOPLOT=0
      shift
      ;;
    --skip-multiqc)
      RUN_MULTIQC=0
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "Unknown argument: $1"
      ;;
  esac
done

[[ -n "$INPUT_DIR" ]] || die "--input is required."
[[ -n "$OUTPUT_ROOT" ]] || die "--output is required."
[[ -n "$PROJECT_ID" ]] || die "--project is required."
[[ -d "$INPUT_DIR" ]] || die "Input directory does not exist: $INPUT_DIR"
[[ "$THREADS" =~ ^[0-9]+$ ]] || die "--threads must be a positive integer."
[[ "$THREADS" -gt 0 ]] || die "--threads must be a positive integer."

if [[ -n "$SAMPLE_SHEET" && ! -f "$SAMPLE_SHEET" ]]; then
  die "Sample sheet does not exist: $SAMPLE_SHEET"
fi

require_cmd seqkit
require_cmd tar
require_cmd gzip

if [[ "$RUN_NANOPLOT" -eq 1 ]] && ! has_cmd NanoPlot; then
  die "NanoPlot not found. Install it or rerun with --skip-nanoplot."
fi

if [[ "$RUN_MULTIQC" -eq 1 ]] && ! has_cmd multiqc; then
  die "multiqc not found. Install it or rerun with --skip-multiqc."
fi

mkdir -p "$OUTPUT_ROOT"

DELIVERY_DIR="$OUTPUT_ROOT/${PROJECT_ID}_delivery"
ARCHIVE="$OUTPUT_ROOT/${PROJECT_ID}_delivery.tar.gz"

if [[ -e "$DELIVERY_DIR" || -e "$ARCHIVE" ]]; then
  if [[ "$OVERWRITE" -eq 1 ]]; then
    if [[ "$RUN_NANOPLOT" -eq 1 && "$REUSE_NANOPLOT" -eq 1 && -d "$DELIVERY_DIR/02_qc_report/nanoplot" ]]; then
      PRESERVED_NANOPLOT_DIR="$(mktemp -d "${TMPDIR:-/tmp}/make_ont_fastq_delivery_nanoplot.XXXXXX")"
      log "Preserving existing NanoPlot reports for reuse."
      cp -pR "$DELIVERY_DIR/02_qc_report/nanoplot/." "$PRESERVED_NANOPLOT_DIR/"
    fi
    log "Removing existing output because --overwrite was supplied."
    rm -rf "$DELIVERY_DIR" "$ARCHIVE"
  else
    die "Output already exists. Remove it first or choose another --project/--output: $DELIVERY_DIR"
  fi
fi

mkdir -p \
  "$DELIVERY_DIR/00_sample_sheet" \
  "$DELIVERY_DIR/01_fastq" \
  "$DELIVERY_DIR/02_qc_report/nanoplot" \
  "$DELIVERY_DIR/03_md5"

if [[ -n "$PRESERVED_NANOPLOT_DIR" ]]; then
  cp -pR "$PRESERVED_NANOPLOT_DIR/." "$DELIVERY_DIR/02_qc_report/nanoplot/"
fi

log "Collecting FASTQ files from: $INPUT_DIR"
FASTQ_LIST="$(mktemp "${TMPDIR:-/tmp}/make_ont_fastq_delivery_fastq_list.XXXXXX")"
find "$INPUT_DIR" -type f \( \
  -name '*.fastq' -o -name '*.fq' -o \
  -name '*.fastq.gz' -o -name '*.fq.gz' \
\) | sort > "$FASTQ_LIST"

FASTQ_COUNT="$(wc -l < "$FASTQ_LIST" | tr -d ' ')"
[[ "$FASTQ_COUNT" -gt 0 ]] || die "No FASTQ files were found in: $INPUT_DIR"

log "Found $FASTQ_COUNT FASTQ file(s). Copying into delivery package."
while IFS= read -r fastq; do
  base="$(basename "$fastq")"
  if [[ -e "$DELIVERY_DIR/01_fastq/$base" ]]; then
    die "Duplicate FASTQ basename detected: $base. Please make filenames unique before delivery."
  fi
  cp -p "$fastq" "$DELIVERY_DIR/01_fastq/"
done < "$FASTQ_LIST"

if [[ -n "$SAMPLE_SHEET" ]]; then
  log "Copying sample sheet."
  cp -p "$SAMPLE_SHEET" "$DELIVERY_DIR/00_sample_sheet/"
fi

log "Testing gzip integrity for compressed FASTQ files."
while IFS= read -r fastq; do
  case "$fastq" in
    *.gz)
      gzip -t "$fastq"
      ;;
  esac
done < "$FASTQ_LIST"

log "Generating seqkit statistics."
(
  cd "$DELIVERY_DIR"
  seqkit stats -a -T 01_fastq/*.f*q* > 02_qc_report/fastq_stats.tsv
)

if [[ "$RUN_NANOPLOT" -eq 1 ]]; then
  log "Generating NanoPlot reports."
  while IFS= read -r delivered_fastq; do
    sample="$(basename "$delivered_fastq")"
    sample="${sample%.gz}"
    sample="${sample%.fastq}"
    sample="${sample%.fq}"
    sample_out="$DELIVERY_DIR/02_qc_report/nanoplot/$sample"
    mkdir -p "$sample_out"
    if [[ "$REUSE_NANOPLOT" -eq 1 ]] && find "$sample_out" -maxdepth 1 -type f -name '*NanoStats.txt' | grep -q .; then
      log "Reusing existing NanoPlot report for: $sample"
      continue
    fi
    NanoPlot \
      --fastq "$delivered_fastq" \
      --outdir "$sample_out" \
      --threads "$THREADS" \
      --prefix "${sample}_"
  done < <(find "$DELIVERY_DIR/01_fastq" -type f \( \
    -name '*.fastq' -o -name '*.fq' -o \
    -name '*.fastq.gz' -o -name '*.fq.gz' \
  \) | sort)
fi

if [[ "$RUN_MULTIQC" -eq 1 ]]; then
  log "Generating MultiQC report."
  MULTIQC_INPUT="$DELIVERY_DIR/02_qc_report"

  if [[ "$RUN_NANOPLOT" -eq 1 ]]; then
    MULTIQC_INPUT="$DELIVERY_DIR/02_qc_report/nanoplot_for_multiqc"
    mkdir -p "$MULTIQC_INPUT"

    while IFS= read -r nanoplot_stats; do
      sample="$(basename "$(dirname "$nanoplot_stats")")"
      cp -p "$nanoplot_stats" "$MULTIQC_INPUT/${sample}_NanoStats.txt"
    done < <(find "$DELIVERY_DIR/02_qc_report/nanoplot" \
      -mindepth 2 \
      -maxdepth 2 \
      -type f \
      -name '*NanoStats.txt' \
      | sort)

    MULTIQC_NANOSTATS_COUNT="$(find "$MULTIQC_INPUT" -type f -name '*_NanoStats.txt' | wc -l | tr -d ' ')"
    [[ "$MULTIQC_NANOSTATS_COUNT" -gt 0 ]] || die "No NanoPlot NanoStats.txt files were found for MultiQC."
  fi

  multiqc \
    "$MULTIQC_INPUT" \
    --outdir "$DELIVERY_DIR/02_qc_report" \
    --filename multiqc_report.html \
    --title "${PROJECT_ID} ONT FASTQ Delivery QC" \
    --force

  if [[ "$RUN_NANOPLOT" -eq 1 ]]; then
    rm -rf "$MULTIQC_INPUT"
  fi
fi

log "Writing README."
cat > "$DELIVERY_DIR/README.txt" <<EOF
Project: ${PROJECT_ID}
Package generated: $(date '+%Y-%m-%d %H:%M:%S %Z')

This package contains demultiplexed Oxford Nanopore FASTQ files and delivery-level QC results.

Directory contents:
  00_sample_sheet/  Sample sheet or barcode/sample metadata, if provided.
  01_fastq/         Demultiplexed FASTQ files. Files were copied without filtering, renaming, or sequence modification.
  02_qc_report/     FASTQ statistics and QC reports.
  03_md5/           MD5 checksums for delivered files.

QC contents:
  fastq_stats.tsv       Summary generated by seqkit stats, including read count, total bases, read length, and quality statistics.
  nanoplot/             Per-sample NanoPlot reports, if NanoPlot was enabled.
  multiqc_report.html   Aggregated MultiQC report, if MultiQC was enabled. NanoPlot statistics are renamed by sample before MultiQC aggregation.
  multiqc_data/         Tables used by MultiQC, if MultiQC was enabled.

Notes:
  No downstream alignment, taxonomic classification, consensus generation, or biological interpretation was performed.
  MD5 checksums can be used to verify file integrity after transfer.
EOF

cat > "$DELIVERY_DIR/README.zh-CN.txt" <<EOF
项目：${PROJECT_ID}
交付包生成时间：$(date '+%Y-%m-%d %H:%M:%S %Z')

本交付包包含已拆分样本的 Oxford Nanopore FASTQ 文件和交付级质量控制结果。

目录内容：
  00_sample_sheet/  样本表或 barcode/样本对应信息（如提供）。
  01_fastq/         已拆分样本的 FASTQ 文件。文件仅复制到交付目录，未进行过滤、重命名或序列内容修改。
  02_qc_report/     FASTQ 基础统计和质量控制报告。
  03_md5/           交付 FASTQ 文件的 MD5 校验值。

质控内容：
  fastq_stats.tsv       由 seqkit stats 生成的 FASTQ 统计表，包括 reads 数量、总碱基数、reads 长度和质量值等信息。
  nanoplot/             每个样本的 NanoPlot 报告（如果启用了 NanoPlot）。
  multiqc_report.html   汇总的 MultiQC 报告（如果启用了 MultiQC）。NanoPlot 统计文件会先按样本名重命名，再用于 MultiQC 汇总。
  multiqc_data/         MultiQC 报告使用的结果表格（如果启用了 MultiQC）。

说明：
  本交付包未进行下游比对、物种分类、共识序列生成或生物学解释。
  MD5 校验值可用于数据传输后的文件完整性验证。
EOF

log "Writing MD5 checksums."
md5_write "$DELIVERY_DIR/01_fastq" "$DELIVERY_DIR/03_md5/md5.txt"

log "Creating archive: $ARCHIVE"
tar -czf "$ARCHIVE" -C "$OUTPUT_ROOT" "$(basename "$DELIVERY_DIR")"

log "Done."
log "Delivery directory: $DELIVERY_DIR"
log "Delivery archive: $ARCHIVE"
