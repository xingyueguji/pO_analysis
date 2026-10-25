#!/bin/bash
# merge_rootfile/download_Oct_22.sh
#
# Copy the merged Oct-22 production from EOS to the Mac, straight over
# xrootd (no lxplus hop: rsync through the lxplus EOS FUSE mount ran at
# 0.25 MB/s on 2026-10-16, xrdcp at ~30 MB/s).
#
#   source : root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/zheng/pO_2026_Oct_22/
#            (written by run_all_hadd_Oct_22.sh on lxplus)
#   dest   : ~/pO_2026_Oct_22/   (DEST= to override)
#
# Per file: skipped when a local copy of the same size exists; otherwise
# copied to <name>.part with end-to-end adler32 verification against the
# source (xrdcp -C adler32:source) and renamed only on success, so an
# interrupted transfer never leaves a file that looks complete. Re-running
# the script resumes where it stopped.
#
# Needs a CERN Kerberos ticket in the Homebrew krb5 cache (Homebrew xrootd
# links Homebrew's MIT krb5, not the system Heimdal):
#   /opt/homebrew/opt/krb5/bin/kinit zheng@CERN.CH
#
# Usage:
#   ./download_Oct_22.sh                         # every .root in the source dir
#   ./download_Oct_22.sh Wp_mu Wm_mu DY_mu_Z     # only names containing a word
#   DRYRUN=1 ./download_Oct_22.sh                # list what would be copied
#   STREAMS=8 DEST=/path ./download_Oct_22.sh
#
# Logs: merge_rootfile/logs_download/<file>.log (xrdcp output), summary on
# stdout. bash-3.2 safe (macOS stock /bin/bash).

set -u

HOST=eoscms.cern.ch
SRC_DIR=/eos/cms/store/group/phys_heavyions/zheng/pO_2026_Oct_22
DEST="${DEST:-$HOME/pO_2026_Oct_22}"
STREAMS="${STREAMS:-8}"
DRYRUN="${DRYRUN:-0}"
KRB=/opt/homebrew/opt/krb5/bin

cd "$(dirname "$0")"
mkdir -p logs_download

for tool in xrdfs xrdcp; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: $tool not found (brew install xrootd)"
        exit 2
    fi
done
if ! "$KRB/klist" -s 2>/dev/null; then
    echo "ERROR: no valid Kerberos ticket. Run:"
    echo "    $KRB/kinit zheng@CERN.CH"
    exit 2
fi

# --- remote listing: "<size> <path>" for every .root file ---
LISTING=$(xrdfs "$HOST" ls -l "$SRC_DIR" 2>&1)
if [[ $? -ne 0 ]]; then
    echo "ERROR: cannot list root://$HOST/$SRC_DIR"
    echo "$LISTING"
    exit 2
fi
FILES=$(echo "$LISTING" | awk '$NF ~ /\.root$/ {print $(NF-1), $NF}')
if [[ -z "$FILES" ]]; then
    echo "ERROR: no .root files in root://$HOST/$SRC_DIR (has run_all_hadd_Oct_22.sh run?)"
    exit 2
fi

mkdir -p "$DEST"
echo "==================================================="
echo "Source : root://$HOST/$SRC_DIR"
echo "Dest   : $DEST"
echo "Filter : ${*:-(all)}"
echo "==================================================="

NOK=0; NSKIP=0; FAILED=""; TOTB=0
START=$(date +%s)
while read -r size path; do
    name=$(basename "$path")
    if [[ $# -gt 0 ]]; then
        keep=0
        for w in "$@"; do [[ "$name" == *"$w"* ]] && keep=1; done
        [[ "$keep" -eq 1 ]] || continue
    fi
    gb=$(awk -v s="$size" 'BEGIN {printf "%.2f", s/1024/1024/1024}')
    local_f="$DEST/$name"

    if [[ -f "$local_f" ]] && [[ "$(stat -f %z "$local_f")" == "$size" ]]; then
        printf '  %-36s %8s GB  already here (same size) -- skipped\n' "$name" "$gb"
        NSKIP=$((NSKIP + 1))
        continue
    fi
    if [[ "$DRYRUN" == "1" ]]; then
        printf '  %-36s %8s GB  would copy\n' "$name" "$gb"
        continue
    fi

    printf '  %-36s %8s GB  copying ... ' "$name" "$gb"
    t0=$(date +%s)
    if xrdcp -f -S "$STREAMS" -C adler32:source --nopbar \
             "root://$HOST/$path" "$local_f.part" > "logs_download/$name.log" 2>&1; then
        mv -f "$local_f.part" "$local_f"
        dt=$(( $(date +%s) - t0 ))
        rate=$(awk -v s="$size" -v t="$dt" 'BEGIN {if (t < 1) t = 1; printf "%.1f", s/1024/1024/t}')
        echo "ok (${dt}s, $rate MB/s, adler32 verified)"
        NOK=$((NOK + 1))
        TOTB=$((TOTB + size))
    else
        echo "FAILED (see logs_download/$name.log)"
        rm -f "$local_f.part"
        FAILED="$FAILED $name"
    fi
done <<< "$FILES"

echo "==================================================="
echo "Done: $NOK copied ($(awk -v s="$TOTB" 'BEGIN {printf "%.1f", s/1024/1024/1024}') GB), $NSKIP already there, in $(( $(date +%s) - START ))s -> $DEST"
if [[ -n "$FAILED" ]]; then
    echo "FAILED:$FAILED   (re-run to retry; a failed copy leaves no file behind)"
    exit 1
fi
