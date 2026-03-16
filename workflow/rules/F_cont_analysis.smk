# -------------------------------------------------------------------------------
# GEP2 - Contamination Analysis Rules
# -------------------------------------------------------------------------------

# Note: The following are defined in the main Snakefile:
#   - FCS_GX_DIR: Path to shared FCS-GX database and tools
#   - TAXDUMP_DIR: os.path.join(DB_FOLDER, "taxdump")
#   - UNIPROT_DIR: os.path.join(DB_FOLDER, "uniprot")
#   - DIAMOND_DB_TYPE: "swissprot" or "reference_proteomes" (from DIAMOND_MODE)
#   - DIAMOND_DB_DIR: os.path.join(UNIPROT_DIR, DIAMOND_DB_TYPE)
#   - DIAMOND_DB_PATH: os.path.join(DIAMOND_DB_DIR, DIAMOND_DB_TYPE + ".dmnd")
#   - FAST_TEMP_DIR: Path to fast temp directory (ideally /dev/shm for FCS-GX)
#   - get_assembly_files(): Get all assembly files for species/assembly
#   - get_assembly_basename(): Extract basename from filepath
#   - is_ncbi_assembly_accession(): Check if value is an NCBI accession
#   - is_url(): Check if value is a URL
#   - _should_skip_analysis(): Check if analysis should be skipped
#   - _is_reads_only_entry(): Check if entry has no assembly
#   - get_assembly_input(): Helper to get assembly file path
#   - samples_config: Parsed sample configuration
#   - normalize_read_type()


# -------------------------------------------------------------------------------
# INPUT FUNCTIONS
# -------------------------------------------------------------------------------

def get_fcs_gx_asm_inputs(wildcards):
    """Get resolved assembly file paths for FCS-GX screening.
    
    Handles local paths, NCBI accessions, and URLs by resolving to
    actual file locations (including downloaded paths).
    Returns a sorted list of paths.
    """
    if not _as_bool(config.get("RUN_FCS", False)):
        return []

    if _should_skip_analysis(wildcards.species, wildcards.asm_id, "fcs"):
        return []

    if _is_reads_only_entry(wildcards.species, wildcards.asm_id):
        return []

    asm_files = get_assembly_files(wildcards.species, wildcards.asm_id)
    resolved = []

    for asm_key, asm_path in sorted(asm_files.items()):
        if not asm_path or asm_path == "None":
            continue

        if is_ncbi_assembly_accession(asm_path) or is_url(asm_path):
            asm_basename = get_assembly_basename(asm_path)
            resolved.append(os.path.join(
                config["OUT_FOLDER"], "GEP2_results", "downloaded_data",
                wildcards.species, "assemblies",
                f"{asm_basename}.fna.gz"
            ))
        else:
            resolved.append(asm_path)

    return resolved


def get_fcs_gx_asm_count(wildcards):
    """Get number of assembly files for this asm_id."""
    asm_files = get_assembly_files(wildcards.species, wildcards.asm_id)
    return len([v for v in asm_files.values() if v and v != "None"])


def get_blob_asm_input(wildcards):
    """Get assembly file for blobtools analysis.

    Resolves the assembly path from wildcards, handling both local files
    and downloaded assemblies (accessions/URLs).
    """
    asm_files = get_assembly_files(wildcards.species, wildcards.asm_id)

    for asm_key, asm_path in asm_files.items():
        if not asm_path or asm_path == "None":
            continue

        if get_assembly_basename(asm_path) == wildcards.asm_basename:
            if is_ncbi_assembly_accession(asm_path) or is_url(asm_path):
                return os.path.join(
                    config["OUT_FOLDER"], "GEP2_results", "downloaded_data",
                    wildcards.species, "assemblies",
                    f"{wildcards.asm_basename}.fna.gz"
                )
            else:
                return asm_path

    raise ValueError(
        f"Could not find assembly for {wildcards.species}/{wildcards.asm_id}/{wildcards.asm_basename}"
    )

def _get_best_read_type_for_blob(species, asm_id):
    """Determine the best available read type for blobtools coverage."""
    try:
        asm_data = samples_config["sp_name"][species]["asm_id"][asm_id]
        read_type_dict = asm_data.get("read_type", {})

        available = set()
        for rt_key, rt_data in read_type_dict.items():
            if not rt_key or rt_key == "None" or not rt_data:
                continue
            rt_normalized = normalize_read_type(rt_key)
            read_files = rt_data.get("read_files", {})
            if any(v and v != "None" for v in read_files.values()):
                available.add(rt_normalized)

        for rt in ["hifi", "ont", "illumina", "10x", "hic"]:
            if rt in available:
                return rt

    except (KeyError, TypeError, AttributeError):
        pass

    return None


def _get_blob_read_files(species, asm_id, target_read_type):
    """Get centralized read files for the given read type.

    Returns a list of file paths. For paired-end types, R1 and R2 files are
    returned as consecutive pairs: [r1_a, r2_a, r1_b, r2_b, ...].
    For single-end types: [reads_a, reads_b, ...].
    """
    files = []
    is_paired = target_read_type in ["illumina", "10x", "hic"]

    try:
        asm_data = samples_config["sp_name"][species]["asm_id"][asm_id]
        read_type_dict = asm_data.get("read_type", {})

        for rt_key, rt_data in read_type_dict.items():
            if not rt_key or rt_key == "None" or not rt_data:
                continue

            rt_normalized = normalize_read_type(rt_key)
            if rt_normalized != target_read_type:
                continue

            read_files = rt_data.get("read_files", {})

            for path_key, path_value in sorted(read_files.items()):
                if not path_value or path_value == "None":
                    continue

                if isinstance(path_value, str) and "," in path_value:
                    paths = [p.strip() for p in path_value.split(",")]
                elif isinstance(path_value, list):
                    paths = path_value
                else:
                    paths = [str(path_value)]

                files.extend(paths)

    except (KeyError, TypeError, AttributeError):
        pass

    return files


def get_blob_reads_input(wildcards):
    """Get read files for blobtools coverage mapping.

    Returns the centralized read paths for the best available read type.
    These are already under .../data/{species}/reads/{read_type}/...
    """
    best_rt = _get_best_read_type_for_blob(wildcards.species, wildcards.asm_id)
    if not best_rt:
        return []

    return _get_blob_read_files(wildcards.species, wildcards.asm_id, best_rt)


def get_blob_read_type_param(wildcards):
    """Return the best read type string for use as a shell param."""
    best_rt = _get_best_read_type_for_blob(wildcards.species, wildcards.asm_id)
    return best_rt if best_rt else "none"


def get_minimap2_preset(wildcards):
    """Return the minimap2 preset flag for the best available read type."""
    presets = {
        "hifi":     "map-hifi",
        "ont":      "map-ont",
        "illumina": "sr",
        "10x":      "sr",
        "hic":      "sr",
    }
    best_rt = _get_best_read_type_for_blob(wildcards.species, wildcards.asm_id)
    return presets.get(best_rt, "sr")


# -------------------------------------------------------------------------------
# RULES
# -------------------------------------------------------------------------------

rule F00_download_fcs_gx_db:
    """Download FCS-GX database, runner script, and container image (runs once)."""
    output:
        flag = os.path.join(FCS_GX_DIR, "fcs_gx_db.done")
    params:
        fcs_dir = FCS_GX_DIR,
        db_dir  = os.path.join(FCS_GX_DIR, "gxdb"),
        sif     = os.path.join(FCS_GX_DIR, "fcs-gx.sif"),
        fcs_py  = os.path.join(FCS_GX_DIR, "fcs.py"),
        manifest_url = "https://ncbi-fcs-gx.s3.amazonaws.com/gxdb/latest/all.manifest"
    threads: cpu_func("samtools_index")
    resources:
        mem_mb = mem_func("samtools_index"),
        runtime = time_func("samtools_index")
    log:
        os.path.join(FCS_GX_DIR, "download_fcs_gx_db.log")
    benchmark:
        os.path.join(FCS_GX_DIR, "download_fcs_gx_db_benchmark.txt")
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1

        echo "[GEP2] =========================================="
        echo "[GEP2] Downloading FCS-GX tools and database"
        echo "[GEP2] Target directory: {params.fcs_dir}"
        echo "[GEP2] =========================================="

        mkdir -p {params.fcs_dir}
        mkdir -p {params.db_dir}

        # Ensure singularity/apptainer is available
        if command -v apptainer &>/dev/null; then
            CONTAINER_CMD="apptainer"
            if ! command -v singularity &>/dev/null; then
                echo "[GEP2] Creating singularity -> apptainer wrapper"
                mkdir -p "$HOME/.local/bin"
                cat > "$HOME/.local/bin/singularity" <<'WRAPPER'
#!/bin/bash
exec apptainer "$@"
WRAPPER
                chmod +x "$HOME/.local/bin/singularity"
                export PATH="$HOME/.local/bin:$PATH"
            fi
        elif command -v singularity &>/dev/null; then
            CONTAINER_CMD="singularity"
        else
            echo "[GEP2] ERROR: Neither apptainer nor singularity found on this system."
            echo "[GEP2] Please install apptainer or load the appropriate module."
            exit 1
        fi
        echo "[GEP2] Using container runtime: $CONTAINER_CMD ($(which $CONTAINER_CMD))"

        # Download fcs.py runner script
        if [ ! -f "{params.fcs_py}" ]; then
            echo "[GEP2] Downloading fcs.py runner script..."
            curl -fsSL https://github.com/ncbi/fcs/raw/main/dist/fcs.py -o "{params.fcs_py}"
            chmod +x "{params.fcs_py}"
        else
            echo "[GEP2] fcs.py already present, skipping download"
        fi

        # Download FCS-GX container image
        if [ ! -f "{params.sif}" ]; then
            echo "[GEP2] Downloading fcs-gx.sif container image..."
            curl -fSL https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/latest/fcs-gx.sif \
                -o "{params.sif}.tmp"
            mv "{params.sif}.tmp" "{params.sif}"
        else
            echo "[GEP2] fcs-gx.sif already present, skipping download"
        fi

        # Download FCS-GX database
        export FCS_DEFAULT_IMAGE="{params.sif}"

        if [ -f "{params.db_dir}/all.gxi" ] && \
           [ -f "{params.db_dir}/all.gxs" ] && \
           [ -f "{params.db_dir}/all.manifest" ]; then
            echo "[GEP2] FCS-GX database files found, verifying integrity..."
            if python3 "{params.fcs_py}" db check \
                --mft "{params.manifest_url}" \
                --dir "{params.db_dir}" 2>&1; then
                echo "[GEP2] Database integrity verified, skipping download"
                touch {output.flag}
                echo "[GEP2] FCS-GX setup complete (database already present)"
                exit 0
            else
                echo "[GEP2] Database verification failed, re-downloading..."
            fi
        fi

        echo "[GEP2] Downloading FCS-GX database (this may take a while, ~470 GiB)..."
        python3 "{params.fcs_py}" db get \
            --mft "{params.manifest_url}" \
            --dir "{params.db_dir}"

        echo "[GEP2] Verifying database integrity..."
        python3 "{params.fcs_py}" db check \
            --mft "{params.manifest_url}" \
            --dir "{params.db_dir}"

        touch {output.flag}
        echo "[GEP2] FCS-GX setup complete (database downloaded and verified)"
        """


rule F01_run_fcs_gx:
    """Run FCS-GX contamination screening for all assemblies of an asm_id.
    
    This rule:
    1. Copies the FCS-GX database to RAM (FAST_TEMP_DIR) for performance
    2. Resolves the species tax-id via the GoaT API
    3. Screens each assembly file (1 or 2 per asm_id)
    4. Cleans up RAM copy
    
    Requires 512 GB RAM. Runs without container directive because fcs.py
    manages its own container via apptainer/singularity.
    """
    input:
        db_flag    = os.path.join(FCS_GX_DIR, "fcs_gx_db.done"),
        assemblies = get_fcs_gx_asm_inputs
    output:
        flag = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "decontamination", "fcs-gx", "fcs_gx.done"
        )
    params:
        fcs_dir      = FCS_GX_DIR,
        db_dir       = os.path.join(FCS_GX_DIR, "gxdb"),
        sif          = os.path.join(FCS_GX_DIR, "fcs-gx.sif"),
        fcs_py       = os.path.join(FCS_GX_DIR, "fcs.py"),
        outdir       = lambda w: os.path.join(
            config["OUT_FOLDER"], "GEP2_results", w.species, w.asm_id,
            "decontamination", "fcs-gx"
        ),
        asm_count    = get_fcs_gx_asm_count,
        species_name = lambda w: w.species
    threads: cpu_func("fcs_gx")
    resources:
        mem_mb = mem_func("fcs_gx"),
        runtime = time_func("fcs_gx")
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "F01_fcs_gx.log"
        )
    benchmark:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "F01_fcs_gx_benchmark.txt"
        )
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        exec > {log} 2>&1

        echo "[GEP2] =========================================="
        echo "[GEP2] FCS-GX contamination screening"
        echo "[GEP2] Species: {wildcards.species}"
        echo "[GEP2] Assembly ID: {wildcards.asm_id}"
        echo "[GEP2] Assembly count: {params.asm_count}"
        echo "[GEP2] =========================================="

        # Ensure singularity/apptainer is available
        if command -v apptainer &>/dev/null; then
            CONTAINER_CMD="apptainer"
            if ! command -v singularity &>/dev/null; then
                echo "[GEP2] Creating singularity -> apptainer wrapper"
                mkdir -p "$HOME/.local/bin"
                cat > "$HOME/.local/bin/singularity" <<'WRAPPER'
#!/bin/bash
exec apptainer "$@"
WRAPPER
                chmod +x "$HOME/.local/bin/singularity"
                export PATH="$HOME/.local/bin:$PATH"
            fi
        elif command -v singularity &>/dev/null; then
            CONTAINER_CMD="singularity"
        else
            echo "[GEP2] ERROR: Neither apptainer nor singularity found."
            exit 1
        fi
        echo "[GEP2] Container runtime: $CONTAINER_CMD"

        export FCS_DEFAULT_IMAGE="{params.sif}"
        mkdir -p {params.outdir}

        # Resolve tax-id from GoaT API
        echo "[GEP2] Resolving tax-id for: {params.species_name}"

    TAX_ID=$(python3 -c "
import json, sys
from urllib.request import urlopen, Request
from urllib.parse import urlencode, quote
from urllib.error import HTTPError

species = '{params.species_name}'.replace('_', ' ')
params = urlencode({{
    'query': 'tax_name(' + species + ')',
    'result': 'taxon'
}}, quote_via=quote)
url = 'https://goat.genomehubs.org/api/v2/search?' + params

print('DEBUG: url = ' + url, file=sys.stderr)

try:
    req = Request(url)
    req.add_header('Accept', 'application/json')
    resp = urlopen(req, timeout=60)
    data = json.loads(resp.read().decode())
    results = data.get('results', [])
    if not results:
        print('ERROR: Species not found in GoaT: ' + species, file=sys.stderr)
        sys.exit(1)
    taxon_id = results[0]['result']['taxon_id']
    print(taxon_id)
except HTTPError as e:
    body = e.read().decode() if e.fp else 'no body'
    print('ERROR: HTTP ' + str(e.code) + ': ' + body[:500], file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print('ERROR: GoaT API request failed: ' + str(e), file=sys.stderr)
    sys.exit(1)
")

        if [ -z "$TAX_ID" ]; then
            echo "[GEP2] ERROR: Could not resolve tax-id for {params.species_name}"
            exit 1
        fi
        echo "[GEP2] Resolved tax-id: $TAX_ID"

        # Copy database to fast storage (RAM) if available
        FAST_TMP="${{GEP2_FAST_TMP:-}}"
        USE_RAM_DB=false
        RAM_DB_DIR=""

        if [ -n "$FAST_TMP" ] && [ "$FAST_TMP" != "$GEP2_TMP" ]; then
            echo "[GEP2] FAST_TEMP_DIR is set: $FAST_TMP"
            echo "[GEP2] Checking available space for FCS-GX database..."

            # Ensure the fast tmp directory exists
            mkdir -p "$FAST_TMP" 2>/dev/null || true

            AVAIL_KB=$(df -k "$FAST_TMP" 2>/dev/null | awk 'NR==2 {{print $4}}' || echo "0")
            AVAIL_GB=$((AVAIL_KB / 1024 / 1024))

            if [ "$AVAIL_GB" -ge 480 ]; then
                echo "[GEP2] Sufficient space (${{AVAIL_GB}} GB). Copying database to RAM..."
                RAM_DB_DIR="${{FAST_TMP}}/GEP2_fcs_gx_db_{wildcards.species}_{wildcards.asm_id}"
                mkdir -p "$RAM_DB_DIR/gxdb"

                echo "[GEP2] Copying FCS-GX database to $RAM_DB_DIR/gxdb ..."
                cp {params.db_dir}/* "$RAM_DB_DIR/gxdb/"
                echo "[GEP2] Database copy to RAM complete"

                USE_RAM_DB=true
                GXDB_LOC="$RAM_DB_DIR/gxdb"
            else
                echo "[GEP2] Insufficient space in FAST_TEMP_DIR (${{AVAIL_GB}} GB < 480 GB)"
                echo "[GEP2] Running FCS-GX from disk (slower but functional)"
                GXDB_LOC="{params.db_dir}"
            fi
        else
            echo "[GEP2] FAST_TEMP_DIR not configured or same as TEMP_DIR"
            echo "[GEP2] Running FCS-GX from disk"
            GXDB_LOC="{params.db_dir}"
        fi

        # Cleanup trap (remove RAM copy on exit)
        cleanup_ram() {{
            if [ "$USE_RAM_DB" = true ] && [ -n "$RAM_DB_DIR" ] && [ -d "$RAM_DB_DIR" ]; then
                echo "[GEP2] Cleaning up RAM database copy: $RAM_DB_DIR"
                rm -rf "$RAM_DB_DIR"
            fi
        }}
        trap cleanup_ram EXIT

        echo "[GEP2] Using GX database: $GXDB_LOC"

        # Screen each assembly
        ASSEMBLIES="{input.assemblies}"
        ASM_COUNT={params.asm_count}
        ASM_IDX=0
        FAILED=0

        for ASM_FILE in $ASSEMBLIES; do
            ASM_IDX=$((ASM_IDX + 1))

            # Extract basename (strip extensions)
            ASM_BASENAME=$(basename "$ASM_FILE")
            for EXT in .fna.gz .fa.gz .fasta.gz .fna .fa .fasta .gz; do
                ASM_BASENAME="${{ASM_BASENAME%$EXT}}"
            done

            echo ""
            echo "[GEP2] Assembly $ASM_IDX/$ASM_COUNT: $ASM_BASENAME"
            echo "[GEP2] File: $ASM_FILE"

            ASM_OUTDIR="{params.outdir}/$ASM_BASENAME"
            mkdir -p "$ASM_OUTDIR"

            # Run FCS-GX screening
            echo "[GEP2] Running FCS-GX screen (tax-id: $TAX_ID)..."
            python3 "{params.fcs_py}" screen genome \
                --fasta "$ASM_FILE" \
                --out-dir "$ASM_OUTDIR" \
                --gx-db "$GXDB_LOC" \
                --tax-id "$TAX_ID" \
            || {{
                echo "[GEP2] FCS-GX failed for $ASM_BASENAME"
                FAILED=$((FAILED + 1))
                continue
            }}

            # Show outputs
            echo "[GEP2] FCS-GX outputs for $ASM_BASENAME:"
            ls -la "$ASM_OUTDIR/" 2>/dev/null || true

            # Print contamination summary
            REPORT=$(find "$ASM_OUTDIR" -name "*.fcs_gx_report.txt" -type f | head -1)
            if [ -n "$REPORT" ] && [ -f "$REPORT" ]; then
                echo ""
                echo "[GEP2] === Contamination Report: $ASM_BASENAME ==="
                head -50 "$REPORT"
                echo "[GEP2] === End Report ==="
            fi

            echo "[GEP2] FCS-GX completed for $ASM_BASENAME"
        done

        echo ""
        if [ "$FAILED" -gt 0 ]; then
            echo "[GEP2] ERROR: $FAILED out of $ASM_COUNT assemblies failed"
            exit 1
        fi

        echo "[GEP2] FCS-GX completed successfully for all $ASM_COUNT assemblies"
        touch {output.flag}
        """

rule F02_download_taxdump:
    """Download NCBI new_taxdump for BlobToolKit (runs once, stored in DB_FOLDER)"""
    output:
        flag = os.path.join(TAXDUMP_DIR, "taxdump.done")
    params:
        taxdump_dir = TAXDUMP_DIR
    threads: cpu_func("download_data")
    resources:
        mem_mb = mem_func("download_data"),
        runtime = time_func("download_data")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(TAXDUMP_DIR, "download_taxdump.log")
    benchmark:
        os.path.join(TAXDUMP_DIR, "download_taxdump_benchmark.txt")
    shell:
        """
        set -euo pipefail
        exec > "{log}" 2>&1

        echo "[GEP2] Downloading NCBI new_taxdump to {params.taxdump_dir}"
        mkdir -p {params.taxdump_dir}

        # Download new_taxdump archive
        TAXDUMP_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
        TARBALL="{params.taxdump_dir}/new_taxdump.tar.gz"

        curl -L --retry 3 --retry-delay 10 -o "$TARBALL" "$TAXDUMP_URL"

        echo "[GEP2] Extracting taxdump archive"
        tar -xzf "$TARBALL" -C {params.taxdump_dir}

        # Verify essential files exist
        for f in nodes.dmp names.dmp; do
            if [ ! -f "{params.taxdump_dir}/$f" ]; then
                echo "[GEP2] ERROR: $f not found after extraction" >&2
                exit 1
            fi
        done

        # Clean up tarball
        rm -f "$TARBALL"

        touch {output.flag}
        echo "[GEP2] NCBI taxdump download complete"
        """


rule F03_build_uniprot_diamond_db:
    """Download UniProt + build Diamond DB with taxonomy (runs once per DB type)"""
    input:
        taxdump_flag = os.path.join(TAXDUMP_DIR, "taxdump.done")
    output:
        flag = os.path.join(DIAMOND_DB_DIR, DIAMOND_DB_TYPE + ".done")
    params:
        db_dir       = DIAMOND_DB_DIR,
        db_type      = DIAMOND_DB_TYPE,
        db_name      = DIAMOND_DB_TYPE,
        taxdump_dir  = TAXDUMP_DIR,
        diamond_db   = DIAMOND_DB_PATH
    threads: cpu_func("build_diamond_db")
    resources:
        mem_mb = mem_func("build_diamond_db"),
        runtime = time_func("build_diamond_db")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(DIAMOND_DB_DIR, "build_diamond_db.log")
    benchmark:
        os.path.join(DIAMOND_DB_DIR, "build_diamond_db_benchmark.txt")
    shell:
        """
        set -euo pipefail
        exec > "{log}" 2>&1

        echo "[GEP2] Building Diamond DB: {params.db_type}"
        echo "[GEP2] Output directory: {params.db_dir}"
        mkdir -p {params.db_dir}

        # -------------------------------------------------------------------
        # SwissProt path  (~90 MB download, ~1 GB Diamond DB)
        # -------------------------------------------------------------------
        if [ "{params.db_type}" = "swissprot" ]; then

            SPROT_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
            FASTA="{params.db_dir}/uniprot_sprot.fasta.gz"

            echo "[GEP2] Downloading SwissProt FASTA"
            curl -L --retry 3 --retry-delay 10 -o "$FASTA" "$SPROT_URL"

            # Build accession -> taxid mapping from FASTA headers
            # Header format: >sp|ACCESSION|ENTRY_NAME Description OS=... OX=TAXID ...
            echo "[GEP2] Extracting accession-to-taxid mapping from headers"
            TAXMAP="{params.db_dir}/accession2taxid.tsv"
            echo -e "accession\\taccession.version\\ttaxid\\tgi" > "$TAXMAP"
            zcat "$FASTA" | grep "^>" | \
                sed -n 's/^>sp|\\([^|]*\\)|.* OX=\\([0-9]*\\) .*/\\1\\t\\1\\t\\2\\t0/p' \
                >> "$TAXMAP"

            NSEQS=$(zcat "$FASTA" | grep -c "^>")
            NMAPPED=$(tail -n +2 "$TAXMAP" | wc -l)
            echo "[GEP2] SwissProt: $NSEQS sequences, $NMAPPED taxid mappings"

            if [ "$NMAPPED" -lt 1000 ]; then
                echo "[GEP2] ERROR: Too few taxid mappings ($NMAPPED). Check FASTA format." >&2
                exit 1
            fi

            # Build Diamond DB with taxonomy
            # Note: diamond makedb auto-appends .dmnd to --db path
            echo "[GEP2] Running diamond makedb (SwissProt)"
            diamond makedb \
                --in "$FASTA" \
                --db "{params.db_dir}/{params.db_name}" \
                --taxonmap "$TAXMAP" \
                --taxonnodes "{params.taxdump_dir}/nodes.dmp" \
                --taxonnames "{params.taxdump_dir}/names.dmp" \
                --threads {threads}

            # Clean up intermediate files (keep .dmnd only)
            rm -f "$FASTA" "$TAXMAP"

        # -------------------------------------------------------------------
        # Reference Proteomes path  (~30-50 GB download, ~30-50 GB Diamond DB)
        # Following the official BlobToolKit recipe:
        #   https://blobtoolkit.genomehubs.org/install/
        # -------------------------------------------------------------------
        elif [ "{params.db_type}" = "reference_proteomes" ]; then

            cd {params.db_dir}

            # Discover the tarball filename (version-stamped, e.g. Reference_Proteomes_2025_01.tar.gz)
            REF_PROT_BASE="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes"
            echo "[GEP2] Discovering Reference Proteomes tarball filename"
            TARBALL_NAME=$(curl -sL "$REF_PROT_BASE/" | grep -oP 'Reference_Proteomes_[0-9_]+\\.tar\\.gz' | head -1)

            if [ -z "$TARBALL_NAME" ]; then
                echo "[GEP2] ERROR: Could not find Reference_Proteomes tarball at $REF_PROT_BASE" >&2
                exit 1
            fi
            echo "[GEP2] Found: $TARBALL_NAME"

            # Download the tarball
            echo "[GEP2] Downloading Reference Proteomes (this may take several hours)"
            curl -L --retry 3 --retry-delay 30 -o "$TARBALL_NAME" "$REF_PROT_BASE/$TARBALL_NAME"

            # Extract
            echo "[GEP2] Extracting tarball"
            tar xf "$TARBALL_NAME"

            # Concatenate all protein FASTA files (exclude DNA and additional/isoform sets)
            echo "[GEP2] Concatenating protein FASTA files"
            touch reference_proteomes.fasta.gz
            find . -mindepth 2 | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional' | \
                xargs cat >> reference_proteomes.fasta.gz

            # Build taxid mapping from idmapping files
            echo "[GEP2] Building taxid mapping"
            echo -e "accession\\taccession.version\\ttaxid\\tgi" > reference_proteomes.taxid_map
            find . -name "*.idmapping.gz" | \
                xargs zcat | grep "NCBI_TaxID" | \
                awk '!seen[$1]++ {{print $1 "\\t" $1 "\\t" $3 "\\t" 0}}' \
                >> reference_proteomes.taxid_map

            NSEQS=$(zcat reference_proteomes.fasta.gz | grep -c "^>")
            NMAPPED=$(tail -n +2 reference_proteomes.taxid_map | wc -l)
            echo "[GEP2] Reference Proteomes: $NSEQS sequences, $NMAPPED taxid mappings"

            # Build Diamond DB with taxonomy
            echo "[GEP2] Running diamond makedb (Reference Proteomes)"
            diamond makedb \
                --in reference_proteomes.fasta.gz \
                --db "{params.db_name}" \
                --taxonmap reference_proteomes.taxid_map \
                --taxonnodes "{params.taxdump_dir}/nodes.dmp" \
                --taxonnames "{params.taxdump_dir}/names.dmp" \
                --threads {threads}

            # Verify the .dmnd file was created
            if [ ! -f "{params.db_name}.dmnd" ]; then
                echo "[GEP2] ERROR: Diamond DB not created" >&2
                exit 1
            fi

            # Clean up large intermediate files (keep .dmnd only)
            echo "[GEP2] Cleaning up intermediate files"
            rm -f "$TARBALL_NAME" reference_proteomes.fasta.gz reference_proteomes.taxid_map
            rm -rf Archaea Bacteria Eukaryota Viruses

            cd -

        else
            echo "[GEP2] ERROR: Unknown DB type '{params.db_type}'" >&2
            exit 1
        fi

        # Verify final .dmnd exists
        if [ ! -f "{params.diamond_db}" ]; then
            echo "[GEP2] ERROR: Expected Diamond DB not found at {params.diamond_db}" >&2
            exit 1
        fi

        DB_SIZE=$(du -h "{params.diamond_db}" | cut -f1)
        echo "[GEP2] Diamond DB size: $DB_SIZE"

        touch {output.flag}
        echo "[GEP2] Diamond DB build complete: {params.db_type}"
        """

rule F04_diamond_blastx:
    """Diamond blastx search for blobtools taxonomic assignment (per assembly)"""
    input:
        asm = get_blob_asm_input,
        db_flag = os.path.join(DIAMOND_DB_DIR, DIAMOND_DB_TYPE + ".done")
    output:
        hits = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "decontamination", "blobtools", "{asm_basename}",
            "diamond_blastx.out"
        )
    params:
        diamond_db  = DIAMOND_DB_PATH.replace(".dmnd", ""),
        sensitivity = DIAMOND_SENSITIVITY,
        db_type     = DIAMOND_DB_TYPE,
        outdir      = lambda w: os.path.join(
            config["OUT_FOLDER"], "GEP2_results", w.species, w.asm_id,
            "decontamination", "blobtools", w.asm_basename
        )
    threads: cpu_func("diamond_blastx")
    resources:
        mem_mb = mem_func("diamond_blastx"),
        runtime = time_func("diamond_blastx")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "F04_diamond_blastx_{asm_basename}.log"
        )
    benchmark:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "F04_diamond_blastx_{asm_basename}_benchmark.txt"
        )
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        exec > "{log}" 2>&1

        echo "[GEP2] Running Diamond blastx for {wildcards.species}/{wildcards.asm_id}/{wildcards.asm_basename}"
        echo "[GEP2] Assembly:    {input.asm}"
        echo "[GEP2] Diamond DB:  {params.diamond_db} ({params.db_type})"
        echo "[GEP2] Sensitivity: {params.sensitivity}"
        echo "[GEP2] Threads:     {threads}"

        mkdir -p {params.outdir}

        # Get work directory for temp files (diamond uses temp space for large queries)
        WORK_DIR="$(gep2_get_workdir 50)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GEP2_diamond_{wildcards.species}_{wildcards.asm_basename}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT

        diamond blastx \
            --query {input.asm} \
            --db {params.diamond_db} \
            --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            {params.sensitivity} \
            --max-target-seqs 1 \
            --evalue 1e-25 \
            --threads {threads} \
            --tmpdir "$TEMP_DIR" \
            --block-size 4 \
            --out {output.hits}

        # Verify output
        if [ ! -f {output.hits} ]; then
            echo "[GEP2] ERROR: Diamond output not created" >&2
            exit 1
        fi

        NHITS=$(wc -l < {output.hits})
        echo "[GEP2] Diamond blastx completed: $NHITS hits"

        if [ "$NHITS" -eq 0 ]; then
            echo "[GEP2] WARNING: Zero hits found. This may indicate a problem with the assembly or database."
        fi

        echo "[GEP2] Diamond blastx complete: {output.hits}"
        """

# Maps the best available reads to the assembly to generate a sorted BAM file
# for blobtools coverage calculation.
rule F05_map_reads_for_blob:
    """Map reads to assembly for blobtools coverage (per assembly basename)"""
    input:
        asm   = get_blob_asm_input,
        reads = get_blob_reads_input
    output:
        bam = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "decontamination", "blobtools", "{asm_basename}",
            "reads_mapped.bam"
        )
    params:
        read_type = get_blob_read_type_param,
        preset    = get_minimap2_preset,
        outdir    = lambda w: os.path.join(
            config["OUT_FOLDER"], "GEP2_results", w.species, w.asm_id,
            "decontamination", "blobtools", w.asm_basename
        )
    threads: cpu_func("map_reads_blob")
    resources:
        mem_mb  = mem_func("map_reads_blob"),
        runtime = time_func("map_reads_blob")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "F05_map_reads_blob_{asm_basename}.log"
        )
    benchmark:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "F05_map_reads_blob_{asm_basename}_benchmark.txt"
        )
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        exec > "{log}" 2>&1

        echo "[GEP2] Mapping reads for blobtools coverage"
        echo "[GEP2] Species:    {wildcards.species}"
        echo "[GEP2] Assembly:   {wildcards.asm_basename}"
        echo "[GEP2] Read type:  {params.read_type}"
        echo "[GEP2] Preset:     {params.preset}"
        echo "[GEP2] Threads:    {threads}"

        mkdir -p {params.outdir}

        # Get work directory for temp files
        WORK_DIR="$(gep2_get_workdir 100)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GEP2_mapblob_{wildcards.species}_{wildcards.asm_basename}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT

        # Split threads: ~75% mapping, ~25% sorting (minimum 2 each)
        TOTAL={threads}
        MAP_THREADS=$(( TOTAL * 3 / 4 - 2 ))
        SORT_THREADS=$(( TOTAL - MAP_THREADS - 2 ))
        VIEW_THREADS=2
        [ "$MAP_THREADS" -lt 2 ] && MAP_THREADS=2
        [ "$SORT_THREADS" -lt 2 ] && SORT_THREADS=2

        echo "[GEP2] Thread split: minimap2=$MAP_THREADS, samtools=$VIEW_THREADS, sambamba=$SORT_THREADS"

        # Determine read inputs based on read type
        READ_TYPE="{params.read_type}"
        ALL_READS=({input.reads})
        NUM_READS=${{#ALL_READS[@]}}

        echo "[GEP2] Input files: $NUM_READS"

        if [ "$READ_TYPE" = "hifi" ] || [ "$READ_TYPE" = "ont" ]; then
            # Long reads: single-end
            # Concatenate all read files as input (minimap2 accepts multiple)
            echo "[GEP2] Long-read mapping: minimap2 -ax {params.preset}"

            minimap2 -ax {params.preset} \
                -t $MAP_THREADS \
                --secondary=no \
                -Q \
                {input.asm} \
                ${{ALL_READS[@]}} \
            | samtools view -bu -@ $VIEW_THREADS \
            | sambamba sort \
                -t $SORT_THREADS \
                -m $(( {resources.mem_mb} * 60 / 100 ))M \
                --tmpdir="$TEMP_DIR" \
                -o {output.bam} \
                /dev/stdin

        elif [ "$READ_TYPE" = "illumina" ] || [ "$READ_TYPE" = "10x" ] || [ "$READ_TYPE" = "hic" ]; then
            # Short/paired reads
            # Multiple read pairs: map each pair, then merge sorted BAMs
            echo "[GEP2] Paired-read mapping: minimap2 -ax {params.preset}"

            if [ "$NUM_READS" -eq 2 ]; then
                # Single pair - pipe directly
                minimap2 -ax {params.preset} \
                    -t $MAP_THREADS \
                    --secondary=no \
                    -Q \
                    {input.asm} \
                    "${{ALL_READS[0]}}" "${{ALL_READS[1]}}" \
                | samtools view -bu -@ $VIEW_THREADS \
                | sambamba sort \
                    -t $SORT_THREADS \
                    -m $(( {resources.mem_mb} * 60 / 100 ))M \
                    --tmpdir="$TEMP_DIR" \
                    -o {output.bam} \
                    /dev/stdin

            elif [ "$((NUM_READS % 2))" -eq 0 ] && [ "$NUM_READS" -gt 2 ]; then
                # Multiple pairs - map each, merge
                PARTIAL_BAMS=()
                PAIR_IDX=0

                for (( i=0; i<NUM_READS; i+=2 )); do
                    PAIR_IDX=$((PAIR_IDX + 1))
                    PARTIAL_BAM="$TEMP_DIR/partial_${{PAIR_IDX}}.bam"
                    echo "[GEP2] Mapping pair $PAIR_IDX: ${{ALL_READS[$i]}} + ${{ALL_READS[$((i+1))]}}"

                    minimap2 -ax {params.preset} \
                        -t $MAP_THREADS \
                        --secondary=no \
                        -Q \
                        {input.asm} \
                        "${{ALL_READS[$i]}}" "${{ALL_READS[$((i+1))]}}" \
                    | samtools view -bu -@ $VIEW_THREADS \
                    | sambamba sort \
                        -t $SORT_THREADS \
                        -m $(( {resources.mem_mb} * 60 / 100 ))M \
                        --tmpdir="$TEMP_DIR" \
                        -o "$PARTIAL_BAM" \
                        /dev/stdin

                    PARTIAL_BAMS+=("$PARTIAL_BAM")
                done

                echo "[GEP2] Merging $PAIR_IDX partial BAMs..."
                sambamba merge \
                    -t {threads} \
                    {output.bam} \
                    "${{PARTIAL_BAMS[@]}}"

            else
                echo "[GEP2] ERROR: Odd number of paired-end read files ($NUM_READS)"
                exit 1
            fi
        else
            echo "[GEP2] ERROR: Unknown read type: $READ_TYPE"
            exit 1
        fi

        # Verify output
        if [ ! -f {output.bam} ]; then
            echo "[GEP2] ERROR: BAM file not created" >&2
            exit 1
        fi

        BAM_SIZE=$(stat -c%s {output.bam})
        echo "[GEP2] BAM size: $((BAM_SIZE / 1024 / 1024)) MB"

        # Quick validation: check BAM header
        if ! sambamba view -H {output.bam} > /dev/null 2>&1; then
            echo "[GEP2] ERROR: Output BAM file is not valid" >&2
            exit 1
        fi

        # Count mapped reads
        MAPPED=$(sambamba flagstat -t {threads} {output.bam} 2>/dev/null | head -1 | awk '{{print $1}}')
        echo "[GEP2] Total alignments: $MAPPED"

        echo "[GEP2] Read mapping complete: {output.bam}"
        """

rule F06_create_blobdir:
    """Create BlobToolKit dataset directory (per assembly basename)"""
    input:
        asm  = get_blob_asm_input,
        hits = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "decontamination", "blobtools", "{asm_basename}",
            "diamond_blastx.out"
        ),
        bam  = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "decontamination", "blobtools", "{asm_basename}",
            "reads_mapped.bam"
        ),
        taxdump_flag = os.path.join(TAXDUMP_DIR, "taxdump.done")
    output:
        flag = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "decontamination", "blobtools", "{asm_basename}",
            "BlobDir", "meta.json"
        )
    params:
        blobdir     = lambda w: os.path.join(
            config["OUT_FOLDER"], "GEP2_results", w.species, w.asm_id,
            "decontamination", "blobtools", w.asm_basename, "BlobDir"
        ),
        taxdump_dir = TAXDUMP_DIR,
        outdir      = lambda w: os.path.join(
            config["OUT_FOLDER"], "GEP2_results", w.species, w.asm_id,
            "decontamination", "blobtools", w.asm_basename
        )
    threads: cpu_func("create_blobdir")
    resources:
        mem_mb  = mem_func("create_blobdir"),
        runtime = time_func("create_blobdir")
    container: CONTAINERS["blobtoolkit"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "F06_create_blobdir_{asm_basename}.log"
        )
    benchmark:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "F06_create_blobdir_{asm_basename}_benchmark.txt"
        )
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        exec > "{log}" 2>&1

        echo "[GEP2] Creating BlobDir dataset"
        echo "[GEP2] Species:     {wildcards.species}"
        echo "[GEP2] Assembly:    {wildcards.asm_basename}"
        echo "[GEP2] FASTA:       {input.asm}"
        echo "[GEP2] Hits file:   {input.hits}"
        echo "[GEP2] BAM file:    {input.bam}"
        echo "[GEP2] Taxdump:     {params.taxdump_dir}"
        echo "[GEP2] BlobDir:     {params.blobdir}"

        # Clean up any previous BlobDir (blobtools doesn't overwrite cleanly)
        if [ -d "{params.blobdir}" ]; then
            echo "[GEP2] Removing previous BlobDir..."
            rm -rf "{params.blobdir}"
        fi

        mkdir -p {params.outdir}

        # Step 1: Create BlobDir from assembly FASTA
        echo ""
        echo "[GEP2] Step 1/3: blobtools create (assembly -> GC + lengths)"

        blobtools create \
            --fasta {input.asm} \
            --taxdump {params.taxdump_dir} \
            {params.blobdir}

        if [ ! -f "{params.blobdir}/meta.json" ]; then
            echo "[GEP2] ERROR: blobtools create failed - no meta.json produced" >&2
            exit 1
        fi

        NSCAFF=$(python3 -c "
import json
with open('{params.blobdir}/meta.json') as f:
    meta = json.load(f)
records = meta.get('records', 0)
print(records)
")
        echo "[GEP2] BlobDir created: $NSCAFF scaffolds"

        # Step 2: Add diamond blastx hits (taxonomic assignment)
        echo ""
        echo "[GEP2] Step 2/3: blobtools add --hits (taxonomic assignment)"

        NHITS=$(wc -l < {input.hits})
        echo "[GEP2] Diamond hits: $NHITS"

        if [ "$NHITS" -gt 0 ]; then
            blobtools add \
                --hits {input.hits} \
                --taxdump {params.taxdump_dir} \
                --taxrule bestsumorder \
                {params.blobdir}

            echo "[GEP2] Taxonomic hits added successfully"
        else
            echo "[GEP2] WARNING: Zero diamond hits - skipping taxonomic assignment"
            echo "[GEP2]          BlobDir will show all scaffolds as 'no-hit'"
        fi

        # Step 3: Add read coverage from BAM
        echo ""
        echo "[GEP2] Step 3/3: blobtools add --cov (read coverage)"
        
        # Re-index BAM as CSI (blobtools prefers .csi over .bai)
        echo "[GEP2] Indexing BAM..."
        samtools index -c {input.bam}

        BAM_SIZE=$(stat -c%s {input.bam} 2>/dev/null || echo "0")
        echo "[GEP2] BAM file size: $((BAM_SIZE / 1024 / 1024)) MB"

        blobtools add \
            --cov {input.bam} \
            --threads {threads} \
            {params.blobdir}

        echo "[GEP2] Coverage data added successfully"

        # Verify final BlobDir
        echo ""
        echo "[GEP2] Verifying BlobDir contents..."

        # Check essential files exist
        if [ ! -f "{params.blobdir}/meta.json" ]; then
            echo "[GEP2] ERROR: BlobDir is missing meta.json" >&2
            exit 1
        fi

        # List BlobDir contents for the log
        echo "[GEP2] BlobDir contents:"
        ls -lh {params.blobdir}/ | head -30

        # Count JSON field files (each represents a data layer)
        NFIELDS=$(find {params.blobdir} -name "*.json" -not -name "meta.json" | wc -l)
        echo "[GEP2] Data fields: $NFIELDS"

        echo ""
        echo "[GEP2] BlobDir created: {params.blobdir}"
        echo "[GEP2]    Scaffolds: $NSCAFF"
        echo "[GEP2]    Diamond hits: $NHITS"
        echo "[GEP2]    Data fields: $NFIELDS"
        """

rule F07_blobplots:
    """Generate blob, snail, and cumulative plots from BlobDir (per assembly)"""
    input:
        meta = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "decontamination", "blobtools", "{asm_basename}",
            "BlobDir", "meta.json"
        )
    output:
        done = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "decontamination", "blobtools", "{asm_basename}",
            "blobplots.done"
        )
    params:
        blobdir = lambda w: os.path.join(
            config["OUT_FOLDER"], "GEP2_results", w.species, w.asm_id,
            "decontamination", "blobtools", w.asm_basename, "BlobDir"
        ),
        outdir  = lambda w: os.path.join(
            config["OUT_FOLDER"], "GEP2_results", w.species, w.asm_id,
            "decontamination", "blobtools", w.asm_basename
        )
    threads: cpu_func("blobplots")
    resources:
        mem_mb  = mem_func("blobplots"),
        runtime = time_func("blobplots")
    container: CONTAINERS["blobtoolkit"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "F07_blobplots_{asm_basename}.log"
        )
    benchmark:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "F07_blobplots_{asm_basename}_benchmark.txt"
        )
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        exec > "{log}" 2>&1

        echo "[GEP2] Generating BlobToolKit plots"
        echo "[GEP2] Species:   {wildcards.species}"
        echo "[GEP2] Assembly:  {wildcards.asm_basename}"
        echo "[GEP2] BlobDir:   {params.blobdir}"
        echo "[GEP2] Output:    {params.outdir}"

        # Snail plot (assembly contiguity summary)
        echo ""
        echo "[GEP2] Generating snail plot..."
        blobtools view \
            --plot \
            --view snail \
            --format png \
            --out {params.outdir} \
            {params.blobdir} \
        || echo "[GEP2] WARNING: Snail plot generation failed (non-fatal)"

        blobtools view \
            --plot \
            --view snail \
            --format svg \
            --out {params.outdir} \
            {params.blobdir} \
        || echo "[GEP2] WARNING: Snail SVG plot generation failed (non-fatal)"

        # Blob plot (GC vs coverage, colored by phylum)
        echo ""
        echo "[GEP2] Generating blob plot..."
        blobtools view \
            --plot \
            --view blob \
            --format png \
            --out {params.outdir} \
            {params.blobdir} \
        || echo "[GEP2] WARNING: Blob plot generation failed (non-fatal)"

        blobtools view \
            --plot \
            --view blob \
            --format svg \
            --out {params.outdir} \
            {params.blobdir} \
        || echo "[GEP2] WARNING: Blob SVG plot generation failed (non-fatal)"

        # Cumulative plot (scaffold length by phylum)
        echo ""
        echo "[GEP2] Generating cumulative plot..."
        blobtools view \
            --plot \
            --view cumulative \
            --format png \
            --out {params.outdir} \
            {params.blobdir} \
        || echo "[GEP2] WARNING: Cumulative plot generation failed (non-fatal)"

        blobtools view \
            --plot \
            --view cumulative \
            --format svg \
            --out {params.outdir} \
            {params.blobdir} \
        || echo "[GEP2] WARNING: Cumulative SVG plot generation failed (non-fatal)"

        # Rename outputs to predictable names
        # blobtools view --plot generates files named after the BlobDir
        # Rename to consistent names for downstream reporting
        echo ""
        echo "[GEP2] Organizing output files..."

        cd {params.outdir}

        # Find and rename generated files (blobtools uses BlobDir name as prefix)
        for PLOT_TYPE in snail blob cumulative; do
            # PNG
            PNG_FILE=$(find . -maxdepth 1 -name "*.$PLOT_TYPE.png" -type f | head -1)
            if [ -n "$PNG_FILE" ] && [ -f "$PNG_FILE" ]; then
                mv "$PNG_FILE" "${{PLOT_TYPE}}_plot.png"
                echo "[GEP2] Created: ${{PLOT_TYPE}}_plot.png"
            fi

            # SVG
            SVG_FILE=$(find . -maxdepth 1 -name "*.$PLOT_TYPE.svg" -type f | head -1)
            if [ -n "$SVG_FILE" ] && [ -f "$SVG_FILE" ]; then
                mv "$SVG_FILE" "${{PLOT_TYPE}}_plot.svg"
                echo "[GEP2] Created: ${{PLOT_TYPE}}_plot.svg"
            fi
        done

        cd -

        # Summary
        echo ""
        echo "[GEP2] Plot files:"
        ls -lh {params.outdir}/*_plot.* 2>/dev/null || echo "[GEP2] No plot files found"

        NPLOTS=$(find {params.outdir} -name "*_plot.*" -type f | wc -l)
        echo "[GEP2] Total plots generated: $NPLOTS"

        if [ "$NPLOTS" -eq 0 ]; then
            echo "[GEP2] WARNING: No plots were generated. Check BlobDir contents."
            echo "[GEP2]          Continuing anyway - BlobDir is still valid for interactive viewing."
        fi

        touch {output.done}
        echo "[GEP2] Blobplots complete: {params.outdir}"
        """
