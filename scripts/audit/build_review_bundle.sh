#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPT_REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$PWD}"
BUNDLE_ID="${1:-review-bundle-$(date -u +%Y%m%dT%H%M%SZ)}"
OUTPUT_DIR="$PROJECT_ROOT/docs/review_bundle"
ZIP_PATH="$OUTPUT_DIR/$BUNDLE_ID.zip"
CHECKSUM_PATH="$OUTPUT_DIR/${BUNDLE_ID}_checksums.tsv"
STAGE_DIR="$(mktemp -d "${TMPDIR:-/tmp}/discover-oncology-bundle-${BUNDLE_ID}-XXXXXX")"

cleanup() {
  rm -rf "$STAGE_DIR"
}

trap cleanup EXIT

mkdir -p "$OUTPUT_DIR"
rm -f "$ZIP_PATH" "$CHECKSUM_PATH"

copy_if_exists() {
  local source_root="$1"
  local rel_path="$2"
  if [[ -e "$source_root/$rel_path" ]]; then
    mkdir -p "$STAGE_DIR/$(dirname "$rel_path")"
    cp -R "$source_root/$rel_path" "$STAGE_DIR/$rel_path"
  fi
}

copy_if_exists "$PROJECT_ROOT" "results"
copy_if_exists "$PROJECT_ROOT" "plots/publication"
copy_if_exists "$PROJECT_ROOT" "logs"
copy_if_exists "$PROJECT_ROOT" "docs/COMPUTE_PLAN.md"
copy_if_exists "$PROJECT_ROOT" "docs/HARMONIZATION_NOTES.md"
copy_if_exists "$PROJECT_ROOT" "docs/audit_runs"
copy_if_exists "$PROJECT_ROOT" "data/manifest.tsv"
copy_if_exists "$PROJECT_ROOT" "cloud_runs"
copy_if_exists "$SCRIPT_REPO_ROOT" "README.md"
copy_if_exists "$SCRIPT_REPO_ROOT" "requirements.txt"
copy_if_exists "$SCRIPT_REPO_ROOT" "scripts"

{
  printf 'relative_path\tsha256\n'
  while IFS= read -r file_path; do
    rel_path="${file_path#$STAGE_DIR/}"
    checksum="$(shasum -a 256 "$file_path" | awk '{print $1}')"
    printf '%s\t%s\n' "$rel_path" "$checksum"
  done < <(find "$STAGE_DIR" -type f | sort)
} > "$CHECKSUM_PATH"

(
  cd "$STAGE_DIR"
  zip -qr "$ZIP_PATH" .
)