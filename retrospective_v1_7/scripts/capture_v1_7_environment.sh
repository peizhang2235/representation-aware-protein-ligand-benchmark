#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PACKAGE_ROOT="$(dirname "$SCRIPT_DIR")"
OUTPUT="$PACKAGE_ROOT/qc/V1_7_COMPUTE_ENVIRONMENT.txt"

mkdir -p "$PACKAGE_ROOT/qc"

{
  printf 'captured_at_utc\t'
  date -u '+%Y-%m-%dT%H:%M:%SZ'
  printf 'operating_system\t'
  uname -a
  printf 'cpu\t'
  sysctl -n machdep.cpu.brand_string 2>/dev/null || printf 'unavailable\n'
  printf 'logical_cores\t'
  sysctl -n hw.logicalcpu 2>/dev/null || printf 'unavailable\n'
  printf 'memory_bytes\t'
  sysctl -n hw.memsize 2>/dev/null || printf 'unavailable\n'
  printf 'Rscript\t'
  command -v Rscript || true
  Rscript --version 2>&1 || true
  printf 'python3\t'
  command -v python3 || true
  python3 --version 2>&1 || true
  printf 'analysis_runtime_observed\tabout 50 seconds on this machine\n'
  printf 'figure_runtime_observed\tabout 30 seconds on this machine\n'
} > "$OUTPUT"

printf 'Wrote %s\n' "$OUTPUT"
