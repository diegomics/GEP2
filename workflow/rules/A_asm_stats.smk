# ═══════════════════════════════════════════════════════════════════════════════
# GEP2 - Assembly Statistics Rules
# ═══════════════════════════════════════════════════════════════════════════════

# Note: The following are defined in the main Snakefile:
#   - BUSCO_DB_DIR: Path to shared BUSCO lineages
#   - get_assembly_input(): Helper to get assembly file path
#   - get_assembly_files(): Helper to get all assembly files for species/assembly
#   - get_assembly_basename(): Extract basename from filepath


# ═══════════════════════════════════════════════════════════════════════════════
# RULES
# ═══════════════════════════════════════════════════════════════════════════════

rule A00_download_compleasm_db:
    """Download compleasm/BUSCO database (runs only once, stored in pipeline folder)"""
    output:
        flag = os.path.join(BUSCO_DB_DIR, "eukaryota_odb12.done"),
        placement_flag = os.path.join(BUSCO_DB_DIR, "placement_files.done")
    params:
        db_dir = BUSCO_DB_DIR
    threads: cpu_func("download_db")
    resources:
        mem_mb = mem_func("download_db"),
        runtime = time_func("download_db")
    container: CONTAINERS["compleasm"]
    log:
        os.path.join(BUSCO_DB_DIR, "download_compleasm_db.log")
    benchmark:
        os.path.join(BUSCO_DB_DIR, "download_compleasm_db_benchmark.txt")
    shell:
        """
        echo "[GEP2] Downloading compleasm eukaryota database to {params.db_dir}"
        mkdir -p {params.db_dir}
        
        compleasm download eukaryota --odb odb12 --library_path {params.db_dir} 2>&1 | tee {log}
        
        touch {output.flag}
        touch {output.placement_flag}
        
        echo "[GEP2] compleasm database download complete"
        """


rule A01_gfastats:
    """Calculate assembly statistics with gfastats"""
    input:
        asm = get_assembly_input
    output:
        stats = "{outdir}/{species}/{assembly}/gfastats/{asm_basename}_stats.txt"
    threads: cpu_func("gfastats")
    resources:
        mem_mb = mem_func("gfastats"),
        runtime = time_func("gfastats")
    container: CONTAINERS["gep2_base"]
    shell:
        """
        echo "[GEP2] Running gfastats on {input.asm}"
        
        mkdir -p $(dirname {output.stats})
        
        gfastats {input.asm} \
            --threads {threads} \
            --nstar-report \
            --discover-paths \
            > {output.stats}
        
        echo "[GEP2] gfastats completed: {output.stats}"
        """


rule A02_compleasm:
    """Assess assembly completeness with compleasm"""
    input:
        asm = get_assembly_input,
        db_flag = os.path.join(BUSCO_DB_DIR, "eukaryota_odb12.done")
    output:
        summary = "{outdir}/{species}/{assembly}/compleasm/{asm_basename}/{asm_basename}_summary.txt",
        archive = "{outdir}/{species}/{assembly}/compleasm/{asm_basename}/{asm_basename}_results.tar.gz"
    params:
        outdir = lambda w: os.path.join(w.outdir, w.species, w.assembly, "compleasm", w.asm_basename),
        shared_db = BUSCO_DB_DIR
    threads: cpu_func("compleasm")
    resources:
        mem_mb = mem_func("compleasm"),
        runtime = time_func("compleasm")
    container: CONTAINERS["compleasm"]
    log:
        "{outdir}/{species}/{assembly}/logs/A02_compleasm_{asm_basename}.log"
    benchmark:
        "{outdir}/{species}/{assembly}/logs/A02_compleasm_{asm_basename}_benchmark.txt"
    shell:
        r'''
        set -euo pipefail
        mkdir -p "{params.outdir}" "$(dirname {log})"
        
        exec > "{log}" 2>&1
        
        echo "[GEP2] Running compleasm on {input.asm}"
        
        # Set up temp directory
        BASE="${{GEP2_TMP:-$TMPDIR}}"
        mkdir -p "$BASE" 2>/dev/null || BASE="$TMPDIR"
        
        WORKDIR="$(mktemp -d "$BASE/GEP2_compleasm_{wildcards.species}_{wildcards.asm_basename}_XXXXXX")"
        TMPDB="$WORKDIR/compleasm_db"
        mkdir -p "$TMPDB"
        
        echo "[GEP2] Preparing isolated lineage DB in $TMPDB"
        cp -rL "{params.shared_db}"/* "$TMPDB/"
        
        echo "[GEP2] Running compleasm --autolineage --retrocopy"
        compleasm run \
            -a "{input.asm}" \
            -o "{params.outdir}" \
            -t {threads} \
            --library_path "$TMPDB" \
            --autolineage \
            --retrocopy
        
        echo "[GEP2] Packaging outputs and cleaning directory"
        TB="{output.archive}"
        cd "{params.outdir}"
        rm -f "$TB"
        
        # Archive everything except summary.txt and the archive itself
        # Use tar -czf directly instead of piping to avoid "short read" errors
        tar -czf "$TB" --exclude='summary.txt' --exclude='*_results.tar.gz' .
        
        # Verify archive
        gzip -t "$TB"
        tar -tzf "$TB" >/dev/null
        
        # Handle summary.txt
        if [ -f "summary.txt" ]; then
            mv summary.txt "{output.summary}"
        else
            echo "[GEP2] Warning: summary.txt not found in output"
            FOUND_SUMMARY=$(find . -name "summary.txt" -type f | head -1)
            if [ -n "$FOUND_SUMMARY" ]; then
                echo "[GEP2] Found summary at: $FOUND_SUMMARY"
                cp "$FOUND_SUMMARY" "{output.summary}"
            else
                echo "[GEP2] Creating empty summary file"
                touch "{output.summary}"
            fi
        fi
        
        # Clean up
        find . -mindepth 1 -maxdepth 1 ! -name '*.txt' ! -name '*_results.tar.gz' -exec rm -rf {{}} \;
        
        echo "[GEP2] Done: archive created at $TB"
        
        cd /
        rm -rf "$WORKDIR"
        
        echo "[GEP2] compleasm completed successfully"
        '''
