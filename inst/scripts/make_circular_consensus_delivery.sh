#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash make_circular_consensus_delivery.sh \
    --consensus /path/to/final_consensus.fasta \
    --output /path/to/delivery_root \
    --project PROJECT_ID \
    [--sequence-type circular_dna|plasmid|virus] \
    [--sample-sheet /path/to/sample_sheet.csv|tsv|xlsx] \
    [--bam /path/to/all_reads_vs_final.bam] \
    [--bai /path/to/all_reads_vs_final.bam.bai] \
    [--depth /path/to/all_reads_vs_final.depth.txt] \
    [--assembly-info /path/to/assembly_info.txt] \
    [--flye-log /path/to/flye.log] \
    [--notes /path/to/notes.txt] \
    [--overwrite]

Description:
  Build a delivery package for a finished ONT circular DNA consensus sequence.
  The script does not annotate, polish, rotate, rename the sequence header, or
  otherwise modify sequence contents. It copies the provided final consensus
  FASTA and optional evidence files into a customer-facing package.

Outputs:
  <output>/<project>_consensus_delivery/
    00_metadata/
      sample sheet and notes, if provided
    01_consensus/
      consensus.fasta
      consensus.fasta.fai, if samtools is available
    02_evidence/
      final BAM/BAI, depth, Flye assembly_info/log, if provided
    03_qc/
      consensus_stats.tsv
    04_md5/
      md5.txt
    manifest.tsv
    README.txt
    README.zh-CN.txt
  <output>/<project>_consensus_delivery.tar.gz

Required:
  seqkit
  tar
  gzip
  md5sum or md5

Optional:
  samtools, for FASTA index generation
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

check_file() {
  local path="$1"
  local label="$2"
  [[ -f "$path" ]] || die "$label does not exist or is not a file: $path"
  [[ -s "$path" ]] || die "$label is empty: $path"
}

copy_file() {
  local source="$1"
  local dest="$2"
  local label="$3"

  cp -p "$source" "$dest"
  printf '%s\t%s\t%s\n' "$label" "$source" "${dest#$DELIVERY_DIR/}" >> "$MANIFEST"
}

md5_write() {
  local target_dir="$1"
  local md5_file="$2"

  if has_cmd md5sum; then
    (
      cd "$target_dir"
      find . -type f ! -path './04_md5/*' | sort | while IFS= read -r file; do
        md5sum "${file#./}"
      done
    ) > "$md5_file"
  elif has_cmd md5; then
    (
      cd "$target_dir"
      find . -type f ! -path './04_md5/*' | sort | while IFS= read -r file; do
        rel="${file#./}"
        sum="$(md5 -q "$rel")"
        printf '%s  %s\n' "$sum" "$rel"
      done
    ) > "$md5_file"
  else
    die "Neither md5sum nor md5 was found."
  fi
}

CONSENSUS=""
OUTPUT_ROOT=""
PROJECT_ID=""
SEQUENCE_TYPE="circular_dna"
SAMPLE_SHEET=""
BAM=""
BAI=""
DEPTH=""
ASSEMBLY_INFO=""
FLYE_LOG=""
NOTES=""
OVERWRITE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --consensus)
      CONSENSUS="${2:-}"
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
    --sequence-type)
      SEQUENCE_TYPE="${2:-}"
      shift 2
      ;;
    --sample-sheet)
      SAMPLE_SHEET="${2:-}"
      shift 2
      ;;
    --bam)
      BAM="${2:-}"
      shift 2
      ;;
    --bai)
      BAI="${2:-}"
      shift 2
      ;;
    --depth)
      DEPTH="${2:-}"
      shift 2
      ;;
    --assembly-info)
      ASSEMBLY_INFO="${2:-}"
      shift 2
      ;;
    --flye-log)
      FLYE_LOG="${2:-}"
      shift 2
      ;;
    --notes)
      NOTES="${2:-}"
      shift 2
      ;;
    --overwrite)
      OVERWRITE=1
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

[[ -n "$CONSENSUS" ]] || die "--consensus is required."
[[ -n "$OUTPUT_ROOT" ]] || die "--output is required."
[[ -n "$PROJECT_ID" ]] || die "--project is required."

case "$SEQUENCE_TYPE" in
  circular_dna|plasmid|virus)
    ;;
  *)
    die "--sequence-type must be one of: circular_dna, plasmid, virus"
    ;;
esac

check_file "$CONSENSUS" "Consensus FASTA"
[[ -n "$SAMPLE_SHEET" ]] && check_file "$SAMPLE_SHEET" "Sample sheet"
[[ -n "$BAM" ]] && check_file "$BAM" "BAM"
[[ -n "$BAI" ]] && check_file "$BAI" "BAI"
[[ -n "$DEPTH" ]] && check_file "$DEPTH" "Depth file"
[[ -n "$ASSEMBLY_INFO" ]] && check_file "$ASSEMBLY_INFO" "Flye assembly_info"
[[ -n "$FLYE_LOG" ]] && check_file "$FLYE_LOG" "Flye log"
[[ -n "$NOTES" ]] && check_file "$NOTES" "Notes file"

require_cmd seqkit
require_cmd tar
require_cmd gzip

mkdir -p "$OUTPUT_ROOT"

DELIVERY_DIR="$OUTPUT_ROOT/${PROJECT_ID}_consensus_delivery"
ARCHIVE="$OUTPUT_ROOT/${PROJECT_ID}_consensus_delivery.tar.gz"
MANIFEST="$DELIVERY_DIR/manifest.tsv"

if [[ -e "$DELIVERY_DIR" || -e "$ARCHIVE" ]]; then
  if [[ "$OVERWRITE" -eq 1 ]]; then
    log "Removing existing output because --overwrite was supplied."
    rm -rf "$DELIVERY_DIR" "$ARCHIVE"
  else
    die "Output already exists. Remove it first or rerun with --overwrite: $DELIVERY_DIR"
  fi
fi

mkdir -p \
  "$DELIVERY_DIR/00_metadata" \
  "$DELIVERY_DIR/01_consensus" \
  "$DELIVERY_DIR/02_evidence" \
  "$DELIVERY_DIR/03_qc" \
  "$DELIVERY_DIR/04_md5"

printf 'label\tsource\tdelivery_path\n' > "$MANIFEST"

log "Copying consensus FASTA."
copy_file "$CONSENSUS" "$DELIVERY_DIR/01_consensus/consensus.fasta" "consensus_fasta"

if has_cmd samtools; then
  log "Creating FASTA index."
  samtools faidx "$DELIVERY_DIR/01_consensus/consensus.fasta"
  printf '%s\t%s\t%s\n' \
    "consensus_fasta_index" \
    "generated by samtools faidx" \
    "01_consensus/consensus.fasta.fai" >> "$MANIFEST"
else
  log "samtools not found; skipping FASTA index generation."
fi

if [[ -n "$SAMPLE_SHEET" ]]; then
  log "Copying sample sheet."
  copy_file "$SAMPLE_SHEET" "$DELIVERY_DIR/00_metadata/$(basename "$SAMPLE_SHEET")" "sample_sheet"
fi

if [[ -n "$NOTES" ]]; then
  log "Copying notes."
  copy_file "$NOTES" "$DELIVERY_DIR/00_metadata/$(basename "$NOTES")" "notes"
fi

if [[ -n "$BAM" ]]; then
  log "Copying BAM evidence."
  copy_file "$BAM" "$DELIVERY_DIR/02_evidence/final_alignment.bam" "final_alignment_bam"
fi

if [[ -n "$BAI" ]]; then
  log "Copying BAM index."
  copy_file "$BAI" "$DELIVERY_DIR/02_evidence/final_alignment.bam.bai" "final_alignment_bai"
fi

if [[ -n "$DEPTH" ]]; then
  log "Copying depth evidence."
  copy_file "$DEPTH" "$DELIVERY_DIR/02_evidence/final_depth.txt" "final_depth"
fi

if [[ -n "$ASSEMBLY_INFO" ]]; then
  log "Copying Flye assembly_info."
  copy_file "$ASSEMBLY_INFO" "$DELIVERY_DIR/02_evidence/flye_assembly_info.txt" "flye_assembly_info"
fi

if [[ -n "$FLYE_LOG" ]]; then
  log "Copying Flye log."
  copy_file "$FLYE_LOG" "$DELIVERY_DIR/02_evidence/flye.log" "flye_log"
fi

log "Generating consensus statistics."
seqkit stats -a -T "$DELIVERY_DIR/01_consensus/consensus.fasta" \
  > "$DELIVERY_DIR/03_qc/consensus_stats.tsv"
printf '%s\t%s\t%s\n' \
  "consensus_stats" \
  "generated by seqkit stats -a -T" \
  "03_qc/consensus_stats.tsv" >> "$MANIFEST"

log "Writing README files."
cat > "$DELIVERY_DIR/README.txt" <<EOF
Project: ${PROJECT_ID}
Package generated: $(date '+%Y-%m-%d %H:%M:%S %Z')
Sequence type: ${SEQUENCE_TYPE}

This package contains a finished ONT consensus sequence for a circular DNA
assembly project. No downstream genome annotation is included in this package.

Directory contents:
  00_metadata/      Sample sheet and project notes, if provided.
  01_consensus/     Final consensus FASTA. The sequence content was copied as provided.
  02_evidence/      Optional assembly/alignment evidence files supplied at packaging time.
  03_qc/            Basic consensus FASTA statistics.
  04_md5/           MD5 checksums for delivered files.

Key files:
  01_consensus/consensus.fasta
      Final delivered consensus sequence.
  01_consensus/consensus.fasta.fai
      FASTA index generated by samtools faidx, if samtools was available.
  03_qc/consensus_stats.tsv
      FASTA statistics generated by seqkit stats.
  manifest.tsv
      List of files included in the package and their source paths.
  04_md5/md5.txt
      MD5 checksums for integrity verification after transfer.

Notes:
  This delivery script does not annotate, polish, rotate, rename sequence
  headers, or modify sequence contents.
  Optional evidence files are copied only when supplied through command-line
  arguments.
EOF

cat > "$DELIVERY_DIR/README.zh-CN.txt" <<EOF
项目：${PROJECT_ID}
交付包生成时间：$(date '+%Y-%m-%d %H:%M:%S %Z')
序列类型：${SEQUENCE_TYPE}

本交付包包含一个已经完成的 ONT 环状 DNA 共识序列。本交付包不包含后续
基因组注释结果。

目录内容：
  00_metadata/      样本表和项目说明（如提供）。
  01_consensus/     最终共识序列 FASTA。序列内容按输入文件原样复制。
  02_evidence/      打包时提供的组装/比对证据文件（可选）。
  03_qc/            共识序列 FASTA 基础统计。
  04_md5/           交付文件的 MD5 校验值。

关键文件：
  01_consensus/consensus.fasta
      最终交付的共识序列。
  01_consensus/consensus.fasta.fai
      samtools faidx 生成的 FASTA 索引（如果环境中有 samtools）。
  03_qc/consensus_stats.tsv
      seqkit stats 生成的 FASTA 统计表。
  manifest.tsv
      交付包内文件列表及其来源路径。
  04_md5/md5.txt
      文件完整性校验值，可用于传输后核对。

说明：
  本脚本不会进行注释、polish、序列旋转、header 改名或任何序列内容修改。
  证据文件仅在命令行参数显式提供时复制进交付包。
EOF

log "Writing MD5 checksums."
md5_write "$DELIVERY_DIR" "$DELIVERY_DIR/04_md5/md5.txt"

log "Creating archive: $ARCHIVE"
tar -czf "$ARCHIVE" -C "$OUTPUT_ROOT" "$(basename "$DELIVERY_DIR")"

log "Done."
log "Delivery directory: $DELIVERY_DIR"
log "Delivery archive: $ARCHIVE"
