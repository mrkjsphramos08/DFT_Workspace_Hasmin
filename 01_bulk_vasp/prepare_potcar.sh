#!/usr/bin/env bash
# OPTIONAL MANUAL STEP in the user's Ubuntu WSL, from this package directory.
# This script is supplied as commands only. It was not executed by Codex.
# It never edits the source datasets and refuses to overwrite an output.
set -euo pipefail
umask 077
# All three literal files were confirmed present in Ubuntu-26.04 on 2026-10-28.
# Their proprietary contents, TITEL values and ENMAX were not read by Codex.
library=/home/mark/vasp_paw_library/PBE64
datasets=(Li_sv Co_pv O)
# Complete ALL existence/readability checks before creating an output file.
for dataset in "${datasets[@]}"; do
    source_file="$library/$dataset/POTCAR"
    if [[ ! -f $source_file || ! -r $source_file || ! -s $source_file ]]; then
        printf 'Missing, unreadable, or empty dataset: %s\n' "$source_file" >&2
        exit 2
    fi
done
if [[ -e POTCAR || -L POTCAR ]]; then
    echo 'POTCAR already exists; refusing to overwrite it.' >&2
    exit 2
fi
# Stage privately on the same filesystem; a hard link publishes only after
# successful concatenation and fails if POTCAR appeared in the meantime.
potcar_tmp=$(mktemp ./.POTCAR.XXXXXXXX)
trap 'rm -f -- "$potcar_tmp"' EXIT
cat -- "$library/Li_sv/POTCAR" "$library/Co_pv/POTCAR" "$library/O/POTCAR" > "$potcar_tmp"
[[ -s $potcar_tmp ]]
ln -T -- "$potcar_tmp" POTCAR
printf 'Created POTCAR in Li_sv, Co_pv, O order. Source files are unchanged.\n'
# Inspect these local metadata fields yourself; do not publish licensed data.
grep 'TITEL' POTCAR
grep 'ENMAX' POTCAR
