#!/usr/bin/env bash
set -euo pipefail

DATA_ROOT=/autodl-fs/data
INCOMING="$DATA_ROOT/rheumlens.incoming_20260621"
FINAL="$DATA_ROOT/rheumlens"
SOURCE_MANIFEST="${1:-}"
DEST_MANIFEST="${2:-}"

test -d "$INCOMING"
test -f "$INCOMING/.SOURCE_TRANSFER_STATUS"
grep -qx 'source_exit_code=0' "$INCOMING/.SOURCE_TRANSFER_STATUS"
test -n "$SOURCE_MANIFEST" -a -f "$SOURCE_MANIFEST"
test -n "$DEST_MANIFEST" -a -f "$DEST_MANIFEST"
diff -u "$SOURCE_MANIFEST" "$DEST_MANIFEST"

if [ -e "$FINAL" ]; then
  echo "Refusing promotion: $FINAL already exists; Supervisor decision required" >&2
  exit 1
fi
mv "$INCOMING" "$FINAL"
printf 'promoted_utc=%s\nsource_manifest=%s\ndest_manifest=%s\n' \
  "$(date -u +%FT%TZ)" "$SOURCE_MANIFEST" "$DEST_MANIFEST" > "$FINAL/PROMOTION_RECORD.txt"
echo "$FINAL"
