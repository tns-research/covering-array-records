#!/usr/bin/env bash
set -e
V="${1:-./ca_verify}"
for f in arrays/CA_*.txt; do
  IFS=_ read -r _ N t k v <<< "$(basename "$f" .txt)"
  "$V" "$f" "$t" "$v" | tail -1
done
