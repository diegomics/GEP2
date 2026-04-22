#!/usr/bin/env bash
#
# get_files_for_yahs.sh
# (this was massively improved using Claude!)
#
# Prepare YaHS scaffolding input from a chromap --pairs file.
# Creates a self-contained output folder with:
#   - the decompressed assembly FASTA (or a symlink, if the input was plain)
#   - the samtools faidx index
#   - a PA5 file derived from the pairs input
#
# Usage:
#   get_files_for_yahs.sh <assembly.fa[.gz]> <read_length> <input.pairs.gz> <output_dir> [threads]
#
# Tools used (samtools and bgzip are execed through the GEP2 hic_analysis
# container when run outside Snakemake; the rest are standard Unix utilities):
#   samtools, bgzip, zcat, awk, grep, sort, comm, cut, wc, head, ln
#
# Notes:
#   - chromap --pairs columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2
#   - YaHS PA5 columns:        readID chrom1 pos1 chrom2 pos2 [mapq1] [mapq2]
#   - Strands in cols 6/7 are non-numeric, so YaHS ignores them.
#     MAPQ filtering must therefore happen upstream at chromap (-q).

set -euo pipefail

# --- argument parsing ---------------------------------------------------------
if [ $# -lt 4 ] || [ $# -gt 5 ]; then
    cat >&2 <<EOF
Usage: $0 <assembly.fa[.gz]> <read_length> <input.pairs.gz> <output_dir> [threads]

Arguments:
  assembly.fa[.gz]   Contig assembly FASTA (.fa/.fasta/.fna, optionally .gz).
  read_length        Hi-C read length (positive integer). This will only be used
                     in the suggested YaHS command printed at the end.
  input.pairs.gz     (b)gzipped pairs file from chromap --pairs.
  output_dir         Directory for all derived files. Created if absent.
                     Contents:
                       <asm_basename>.fa      (decompressed copy or symlink)
                       <asm_basename>.fa.fai  (samtools faidx index)
                       <pairs_basename>.pa5   (YaHS input)
  threads            Optional. Threads for bgzip decompression (default: 1).
                     Only the pairs decompression benefits from threads;
                     samtools faidx is single-threaded regardless.
EOF
    exit 1
fi

ASM_IN="$1"
READ_LEN="$2"
PAIRS="$3"
OUTDIR="$4"
THREADS="${5:-1}"

# --- input validation ---------------------------------------------------------
[ -s "$ASM_IN" ] || { echo "[GEP2][ERROR] Assembly not found or empty: $ASM_IN" >&2; exit 1; }
[ -s "$PAIRS"  ] || { echo "[GEP2][ERROR] Pairs file not found or empty: $PAIRS" >&2; exit 1; }
case "$READ_LEN" in
    ''|*[!0-9]*) echo "[GEP2][ERROR] read_length must be a positive integer, got: '$READ_LEN'" >&2; exit 1 ;;
esac
[ "$READ_LEN" -gt 0 ] || { echo "[GEP2][ERROR] read_length must be > 0" >&2; exit 1; }
case "$THREADS" in
    ''|*[!0-9]*) echo "[GEP2][ERROR] threads must be a positive integer, got: '$THREADS'" >&2; exit 1 ;;
esac
[ "$THREADS" -gt 0 ] || { echo "[GEP2][ERROR] threads must be > 0" >&2; exit 1; }

mkdir -p "$OUTDIR"

# --- MAPQ warning (printed loudly at the start) -------------------------------
cat <<'EOF'

==============================================================
  IMPORTANT!!! MAPQ filtering
  Chromap --pairs output contains NO MAPQ column, so YaHS
  cannot filter on MAPQ. For Hi-C scaffolding, MAPQ should
  be at least 10. If not, set HIC_MAPQ: 10 (or higher) in the
  GEP2 config and rerun the chromap step before regenerating
  this PA5 file. Running with HIC_MAPQ: 0 may produce noisy
  scaffolding results.
==============================================================

EOF

# --- derive output filenames from inputs --------------------------------------
# Assembly in outdir: <asm_basename without trailing .gz>
ASM_BASE=$(basename "$ASM_IN")
ASM_BASE="${ASM_BASE%.gz}"
ASM_FOR_YAHS="${OUTDIR}/${ASM_BASE}"

# PA5 name: strip .gz then .pairs from pairs basename, add .pa5
PAIRS_BASE=$(basename "$PAIRS")
PAIRS_BASE="${PAIRS_BASE%.gz}"
PAIRS_BASE="${PAIRS_BASE%.pairs}"
OUT_PA5="${OUTDIR}/${PAIRS_BASE}.pa5"

# --- locate script, GEP2 install, and container -------------------------------
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
# Script lives at: <GEP2_FOLDER>/workflow/scripts/get_files_for_yahs.sh
GEP2_FOLDER="$(cd -- "${SCRIPT_DIR}/../.." &> /dev/null && pwd)"

# Resolve the hic_analysis container URI from the canonical config file
# (workflow/envs/containers.yaml) so the version stays in sync automatically.
# Falls back to the known-good URI if the file is missing or unparseable.
CONTAINERS_YAML="${GEP2_FOLDER}/workflow/envs/containers.yaml"
HIC_URI_FALLBACK="docker://diegomics/hic_analysis:0.2"
HIC_URI=""
if [ -f "$CONTAINERS_YAML" ]; then
    # Match a line like:   hic_analysis: "docker://..."   (quotes optional, any indent)
    HIC_URI=$(awk '
        /^[[:space:]]*hic_analysis[[:space:]]*:/ {
            sub(/^[^:]*:[[:space:]]*/, "")   # drop "  hic_analysis: "
            sub(/[[:space:]]*#.*$/, "")      # drop trailing comment
            sub(/[[:space:]]+$/, "")         # trim trailing whitespace
            gsub(/^["'\'']|["'\'']$/, "")    # strip surrounding quotes (after comment removed)
            print
            exit
        }' "$CONTAINERS_YAML")
fi
if [ -z "$HIC_URI" ]; then
    echo "[GEP2][WARNING] Could not parse hic_analysis from $CONTAINERS_YAML" >&2
    echo "[GEP2]   Falling back to: $HIC_URI_FALLBACK" >&2
    HIC_URI="$HIC_URI_FALLBACK"
fi

HIC_VAR=$(echo -n "$HIC_URI" | md5sum | awk '{print $1}')
HIC_CONTAINER="${GEP2_FOLDER}/.snakemake/singularity/${HIC_VAR}.simg"

# Detect whether we are already inside a Snakemake-launched container
inside_container() {
    [ -n "${SINGULARITY_CONTAINER:-}" ] || \
    [ -n "${APPTAINER_CONTAINER:-}"  ] || \
    [ -e /.singularity.d ]
}

# Pick whichever runtime is available
APPTAINER_CMD=""
if ! inside_container; then
    for cmd in apptainer singularity; do
        if command -v "$cmd" >/dev/null 2>&1; then
            APPTAINER_CMD="$cmd"
            break
        fi
    done
fi

# Build dedup'd bind args from all touched dirs (bash 3 compatible, no assoc arrays)
# Note: even when ASM_IN is plain and we symlink it into OUTDIR, the symlink
# target still lives in ASM_IN's directory, so we must bind that dir as well.
APPTAINER_BIND_ARGS=()
_SEEN=""
for p in "$ASM_IN" "$PAIRS" "$OUTDIR"; do
    d="$(cd -- "$(dirname -- "$p")" 2>/dev/null && pwd)" || continue
    # OUTDIR itself is a directory, not a file, so dirname would go one up
    # handle both cases by testing whether the path itself is a directory.
    if [ -d "$p" ]; then
        d="$(cd -- "$p" 2>/dev/null && pwd)"
    fi
    case ":${_SEEN}:" in
        *":${d}:"*) ;;  # already added
        *) APPTAINER_BIND_ARGS+=( -B "${d}:${d}" ); _SEEN="${_SEEN}:${d}" ;;
    esac
done

# Generic wrapper: run a command directly inside container, or via apptainer outside
run_in_container() {
    if inside_container; then
        "$@"
    else
        [ -n "$APPTAINER_CMD" ] || {
            echo "[GEP2][ERROR] Neither 'apptainer' nor 'singularity' on PATH, and not running inside a container." >&2
            exit 1
        }
        [ -s "$HIC_CONTAINER" ] || {
            echo "[GEP2][ERROR] Container image not found: $HIC_CONTAINER" >&2
            echo "[GEP2]   Has the GEP2 pipeline been run at least once to pull containers?" >&2
            exit 1
        }
        "$APPTAINER_CMD" exec "${APPTAINER_BIND_ARGS[@]}" "$HIC_CONTAINER" "$@"
    fi
}

echo "[GEP2] Preparing YaHS input"
echo "[GEP2]   Assembly:    $ASM_IN"
echo "[GEP2]   Read length: $READ_LEN"
echo "[GEP2]   Pairs file:  $PAIRS"
echo "[GEP2]   Output dir:  $OUTDIR"
echo "[GEP2]   Threads:     $THREADS"
if inside_container; then
    echo "[GEP2]   Mode:        inside container (host tools)"
else
    echo "[GEP2]   Mode:        standalone ($APPTAINER_CMD exec)"
    echo "[GEP2]   GEP2_FOLDER: $GEP2_FOLDER"
    echo "[GEP2]   Container:   $HIC_CONTAINER"
    echo "[GEP2]   URI:         $HIC_URI"
fi

# --- assembly: stage into OUTDIR (decompress or symlink), then index ----------
# YaHS reads plain FASTA. We always work from a file inside OUTDIR, either a
# decompressed copy (if input was .gz) or a symlink to the original (if plain).
# This keeps all derived artefacts in one place and avoids polluting the input
# assembly's directory with .fai files.

# Absolute path of the input assembly for use as symlink target
ASM_IN_ABS="$(cd -- "$(dirname -- "$ASM_IN")" && pwd)/$(basename -- "$ASM_IN")"

if [[ "$ASM_IN" == *.gz ]]; then
    if [ ! -s "$ASM_FOR_YAHS" ]; then
        echo "[GEP2] Decompressing assembly to $ASM_FOR_YAHS"
        zcat "$ASM_IN" > "${ASM_FOR_YAHS}.part"
        mv "${ASM_FOR_YAHS}.part" "$ASM_FOR_YAHS"
    else
        echo "[GEP2] Decompressed assembly already exists: $ASM_FOR_YAHS"
    fi
else
    if [ ! -e "$ASM_FOR_YAHS" ]; then
        echo "[GEP2] Symlinking assembly into output dir: $ASM_FOR_YAHS -> $ASM_IN_ABS"
        ln -s "$ASM_IN_ABS" "$ASM_FOR_YAHS"
    else
        echo "[GEP2] Assembly already staged in output dir: $ASM_FOR_YAHS"
    fi
fi

if [ ! -s "${ASM_FOR_YAHS}.fai" ]; then
    echo "[GEP2] Indexing assembly with samtools faidx (single-threaded)"
    run_in_container samtools faidx "$ASM_FOR_YAHS"
else
    echo "[GEP2] FASTA index already exists: ${ASM_FOR_YAHS}.fai"
fi

# --- pairs -> PA5 -------------------------------------------------------------
TMP_PA5="${OUT_PA5}.part"
echo "[GEP2] Stripping header from pairs file (bgzip -@ $THREADS)"
run_in_container bgzip -@ "$THREADS" -dc "$PAIRS" | grep -v '^#' > "$TMP_PA5"

# --- sanity checks ------------------------------------------------------------
N_RECORDS=$(wc -l < "$TMP_PA5")
echo "[GEP2] PA5 records written: $N_RECORDS"
if [ "$N_RECORDS" -eq 0 ]; then
    rm -f "$TMP_PA5"
    echo "[GEP2][ERROR] No records in PA5 -- check input pairs file." >&2
    exit 1
fi

N_COLS=$(awk 'NR==1 {print NF; exit}' "$TMP_PA5")
if [ "$N_COLS" -lt 5 ]; then
    rm -f "$TMP_PA5"
    echo "[GEP2][ERROR] PA5 needs at least 5 columns, got $N_COLS." >&2
    exit 1
fi

# Are columns 6/7 numeric? If yes, YaHS reads them as MAPQ.
if [ "$N_COLS" -ge 7 ]; then
    NUMERIC_67=$(awk 'NR==1 {print ($6 ~ /^[0-9]+$/ && $7 ~ /^[0-9]+$/) ? "yes" : "no"; exit}' "$TMP_PA5")
    if [ "$NUMERIC_67" = "yes" ]; then
        echo "[GEP2] Columns 6/7 numeric -- YaHS will treat them as per-read MAPQ."
    else
        echo "[GEP2] Columns 6/7 non-numeric (strand) -- YaHS ignores them."
        echo "[GEP2]   MAPQ filtering must happen upstream at chromap (HIC_MAPQ)."
    fi
fi

# Cross-check contig names against the .fai. Sampling the first CHECK_LIMIT
# records keeps this fast on large files while still catching the realistic
# failure mode: mapping against a different assembly than the one passed here.
# Contigs referenced only deep in the file will be caught by YaHS itself on
# parse.
CHECK_LIMIT=1000000
echo "[GEP2] Cross-checking PA5 contig names vs assembly index (sampling first $CHECK_LIMIT records)"
MISSING=$(
    comm -23 \
        <(head -n "$CHECK_LIMIT" "$TMP_PA5" | awk '{print $2"\n"$4}' | sort -u) \
        <(cut -f1 "${ASM_FOR_YAHS}.fai" | sort -u)
)
if [ -n "$MISSING" ]; then
    N_MISSING=$(printf '%s\n' "$MISSING" | wc -l)
    echo "[GEP2][WARNING] $N_MISSING contig name(s) in PA5 sample not found in assembly:" >&2
    printf '%s\n' "$MISSING" | head -5 >&2
    [ "$N_MISSING" -gt 5 ] && echo "  ... ($((N_MISSING - 5)) more)" >&2
fi

# Atomic move so partial outputs never look complete
mv "$TMP_PA5" "$OUT_PA5"

# --- summary ------------------------------------------------------------------
echo
echo "[GEP2] Output directory contents:"
echo "[GEP2]   Assembly:    $ASM_FOR_YAHS"
echo "[GEP2]   Index:       ${ASM_FOR_YAHS}.fai"
echo "[GEP2]   PA5:         $OUT_PA5"
echo
echo "Suggested YaHS command (run from any directory):"
echo "  yahs --file-type PA5 --read-length $READ_LEN \\"
echo "       -o ${OUTDIR}/yahs_out \\"
echo "       $ASM_FOR_YAHS \\"
echo "       $OUT_PA5"