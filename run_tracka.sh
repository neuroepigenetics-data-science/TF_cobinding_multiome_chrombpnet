#!/usr/bin/env bash
# Host-side wrapper: runs the Track A container with the right mounts.
#   ./run_tracka.sh check      # validate inputs only
#   ./run_tracka.sh all        # full pipeline
#   ./run_tracka.sh shell      # interactive
# Override the image with IMAGE=..., e.g. IMAGE=zamboni-tracka:1.0-archr
set -euo pipefail
cd "$(dirname "$0")"

IMAGE=${IMAGE:-zamboni-tracka:1.0}
OBJ=${OBJ:-multiome_obj_clean_240326.rds}
CMD=${1:-all}

[ -d data ] || { echo "no data/ directory here" >&2; exit 1; }
[ -f "$OBJ" ] || { echo "merged object not found: $OBJ (set OBJ=...)" >&2; exit 1; }

# Created up front so they are owned by the invoking user, not by root inside
# the container.
mkdir -p ML markers logs

exec docker run --rm -it \
  -v "$PWD/data:/work/data:ro" \
  -v "$PWD/$OBJ:/work/$OBJ:ro" \
  -v "$PWD/meta:/work/meta:ro" \
  -v "$PWD/ML:/work/ML" \
  -v "$PWD/markers:/work/markers" \
  -v "$PWD/logs:/work/logs" \
  "$IMAGE" "$CMD"
