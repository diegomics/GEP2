# -------------------------------------------------------------------------------
# GEP2 - Report Generation Rules
# -------------------------------------------------------------------------------

# Note: The following are defined in the main Snakefile:
#   - get_assembly_files(): Get assembly files for species/assembly
#   - get_assembly_basename(): Extract basename from filepath
#   - get_priority_read_type(): Get highest priority read type for species
#   - get_kmer_length(): Get k-mer length based on config/read type
#   - _as_bool(): Convert config value to boolean


# -------------------------------------------------------------------------------
# INPUT FUNCTIONS
# -------------------------------------------------------------------------------

def get_report_gfastats_inputs(wildcards):
    """Get all gfastats output files for this assembly."""
    asm_files = get_assembly_files(wildcards.species, wildcards.asm_id)
    
    gfastats_files = []
    for asm_key, asm_path in sorted(asm_files.items()):
        if asm_path and asm_path != "None":
            asm_basename = get_assembly_basename(asm_path)
            gfastats_files.append(os.path.join(
                config["OUT_FOLDER"], "GEP2_results", wildcards.species, 
                wildcards.asm_id, "gfastats", f"{asm_basename}_stats.txt"
            ))
    
    return gfastats_files


def get_report_compleasm_inputs(wildcards):
    """Get all compleasm tar.gz files for this assembly (if enabled).
    The tar.gz contains {lineage}_odb{version}/full_table.tsv which will
    be extracted at runtime for --compleasm-full."""
    if not _as_bool(config.get("RUN_COMPL", True)):
        return []
    
    asm_files = get_assembly_files(wildcards.species, wildcards.asm_id)
    
    compleasm_files = []
    for asm_key, asm_path in sorted(asm_files.items()):
        if asm_path and asm_path != "None":
            asm_basename = get_assembly_basename(asm_path)
            compleasm_files.append(os.path.join(
                config["OUT_FOLDER"], "GEP2_results", wildcards.species,
                wildcards.asm_id, "compleasm", asm_basename, 
                f"{asm_basename}_results.tar.gz"
            ))
    
    return compleasm_files


def get_report_merqury_inputs(wildcards):
    """Get Merqury output files if k-mer analysis was run."""
    # Global toggle
    if not _as_bool(config.get("KMER_STATS", True)):
        return {'qv': [], 'completeness': []}
    
    # Per-assembly skip
    if _should_skip_analysis(wildcards.species, wildcards.asm_id, "kmer"):
        return {'qv': [], 'completeness': []}
    
    # Check if this assembly has reads
    try:
        asm_data = samples_config["sp_name"][wildcards.species]["asm_id"][wildcards.asm_id]
        read_type_dict = asm_data.get("read_type", {})
        
        has_any_reads = False
        for read_type, rt_data in read_type_dict.items():
            if read_type and read_type != "None" and rt_data:
                read_files = rt_data.get("read_files", {})
                if any(v and v != "None" for v in read_files.values()):
                    has_any_reads = True
                    break
        
        if not has_any_reads:
            return {'qv': [], 'completeness': []}
            
    except (KeyError, TypeError, AttributeError):
        return {'qv': [], 'completeness': []}
    
    merqury_dir = os.path.join(
        config["OUT_FOLDER"], "GEP2_results", wildcards.species, 
        wildcards.asm_id, "merqury"
    )
    
    return {
        'qv': [os.path.join(merqury_dir, f"{wildcards.asm_id}.qv")],
        'completeness': [os.path.join(merqury_dir, f"{wildcards.asm_id}.completeness.stats")]
    }


def get_report_genomescope_input(wildcards):
    """Get GenomeScope2 plot if k-mer analysis was run."""
    # Global toggle
    if not _as_bool(config.get("KMER_STATS", True)):
        return []
    
    # Per-assembly skip
    if _should_skip_analysis(wildcards.species, wildcards.asm_id, "kmer"):
        return []
    
    read_type = get_priority_read_type_for_assembly(wildcards.species, wildcards.asm_id)
    
    if not read_type:
        return []
    
    kmer_len = get_kmer_length(read_type)
    
    return [os.path.join(
        config["OUT_FOLDER"], "GEP2_results", wildcards.species,
        wildcards.asm_id, f"k{kmer_len}", "genomescope2", f"{wildcards.asm_id}_linear_plot.png"
    )]


def get_report_inspector_inputs(wildcards):
    """Get Inspector output files if long read analysis was run."""
    # Global toggle
    if not _as_bool(config.get("RUN_INSP", True)):
        return []
    
    # Per-assembly skip
    if _should_skip_analysis(wildcards.species, wildcards.asm_id, "insp"):
        return []
    
    # Check for long reads
    long_rt = _get_long_read_type_for_assembly(wildcards.species, wildcards.asm_id)
    if not long_rt:
        return []
    
    # Get Inspector outputs for each assembly file
    asm_files = get_assembly_files(wildcards.species, wildcards.asm_id)
    results = []
    
    for asm_key, asm_path in asm_files.items():
        if not asm_path or asm_path == "None":
            continue
        
        asm_basename = get_assembly_basename(asm_path)
        inspector_dir = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", wildcards.species,
            wildcards.asm_id, "inspector", asm_basename
        )
        results.append(os.path.join(inspector_dir, "summary_statistics"))
    
    return results


def get_report_hic_inputs(wildcards):
    """Get Hi-C analysis output files if analysis was run."""
    # Global toggle
    if not _as_bool(config.get("RUN_HIC", True)):
        return []
    
    # Per-assembly skip
    if _should_skip_analysis(wildcards.species, wildcards.asm_id, "hic"):
        return []
    
    # Check for Hi-C reads
    if not _has_hic_reads_for_assembly(wildcards.species, wildcards.asm_id):
        return []
    
    # Get Hi-C outputs for each assembly file
    asm_files = get_assembly_files(wildcards.species, wildcards.asm_id)
    results = []
    
    for asm_key, asm_path in asm_files.items():
        if not asm_path or asm_path == "None":
            continue
        
        asm_basename = get_assembly_basename(asm_path)
        hic_dir = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", wildcards.species,
            wildcards.asm_id, "hic", asm_basename
        )
        
        results.extend([
            os.path.join(hic_dir, f"{asm_basename}.pretext"),
            os.path.join(hic_dir, f"{asm_basename}.pairtools_stats.txt")
        ])
    
    return results


def get_report_hic_snapshots(wildcards):
    """Get Hi-C snapshot PNGs if available (only when not in high-res mode).
    
    Note: These are NOT included as required inputs since PretextSnapshot may fail.
    The script will check if files exist at runtime.
    """
    # Global toggle
    if not _as_bool(config.get("RUN_HIC", True)):
        return []
    
    # Snapshots are not created in high-res mode
    if _as_bool(config.get("HIC_HIGH_RES", False)):
        return []
    
    # Per-assembly skip
    if _should_skip_analysis(wildcards.species, wildcards.asm_id, "hic"):
        return []
    
    # Check for Hi-C reads
    if not _has_hic_reads_for_assembly(wildcards.species, wildcards.asm_id):
        return []
    
    # Get snapshot PNGs for each assembly file
    asm_files = get_assembly_files(wildcards.species, wildcards.asm_id)
    results = []
    
    for asm_key, asm_path in sorted(asm_files.items()):
        if not asm_path or asm_path == "None":
            continue
        
        asm_basename = get_assembly_basename(asm_path)
        snapshot_path = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", wildcards.species,
            wildcards.asm_id, "hic", asm_basename,
            f"{asm_basename}_snapshots", f"{asm_basename}_FullMap.png"
        )
        results.append(snapshot_path)
    
    return results


def get_blobplot_inputs(wildcards):
    """Get all blobplot files for this assembly (if enabled)."""
    if not _as_bool(config.get("RUN_BLOB", False)):
        return []
    
    asm_files = get_assembly_files(wildcards.species, wildcards.asm_id)
    
    blob_files = []
    for asm_key, asm_path in sorted(asm_files.items()):
        if asm_path and asm_path != "None":
            asm_basename = get_assembly_basename(asm_path)
            blob_files.append(os.path.join(
                config["OUT_FOLDER"], "GEP2_results", wildcards.species,
                wildcards.asm_id, "decontamination", "blobtools", asm_basename, 
                "blobplots.done"
            ))
    
    return blob_files

def get_blobplot_dirs(wildcards):
    """Get blobtools output directories for the reporting script."""
    if not _as_bool(config.get("RUN_BLOB", False)):
        return []
    
    asm_files = get_assembly_files(wildcards.species, wildcards.asm_id)
    
    blob_dirs = []
    for asm_key, asm_path in sorted(asm_files.items()):
        if asm_path and asm_path != "None":
            asm_basename = get_assembly_basename(asm_path)
            blob_dirs.append(os.path.join(
                config["OUT_FOLDER"], "GEP2_results", wildcards.species,
                wildcards.asm_id, "decontamination", "blobtools", asm_basename
            ))
    
    return blob_dirs


def get_report_fcs_inputs(wildcards):
    """Get fcs-gx done flag for this assembly (if enabled).
    The actual report filenames contain an unpredictable taxid,
    so we depend on the sentinel flag and find files at runtime."""
    if not _as_bool(config.get("RUN_FCS", False)):
        return []
    
    return [os.path.join(
        config["OUT_FOLDER"], "GEP2_results", wildcards.species,
        wildcards.asm_id, "decontamination", "fcs-gx", "fcs_gx.done"
    )]


def get_report_fcs_dirs(wildcards):
    """Get per-assembly fcs-gx directories for runtime file discovery."""
    if not _as_bool(config.get("RUN_FCS", False)):
        return []
    
    asm_files = get_assembly_files(wildcards.species, wildcards.asm_id)
    
    fcs_dirs = []
    for asm_key, asm_path in sorted(asm_files.items()):
        if asm_path and asm_path != "None":
            asm_basename = get_assembly_basename(asm_path)
            fcs_dirs.append(os.path.join(
                config["OUT_FOLDER"], "GEP2_results", wildcards.species,
                wildcards.asm_id, "decontamination", "fcs-gx", asm_basename
            ))
    
    return fcs_dirs


def get_all_report_inputs(wildcards):
    """Collect all inputs for the report rule."""
    inputs = []
    
    # Always need gfastats
    inputs.extend(get_report_gfastats_inputs(wildcards))
    
    # Compleasm if enabled
    inputs.extend(get_report_compleasm_inputs(wildcards))
    
    # Merqury if enabled and has reads
    merqury = get_report_merqury_inputs(wildcards)
    inputs.extend(merqury['qv'])
    inputs.extend(merqury['completeness'])
    
    # GenomeScope2 if enabled
    inputs.extend(get_report_genomescope_input(wildcards))

    # Inspector if enabled and has long reads
    inputs.extend(get_report_inspector_inputs(wildcards))
    
    # Hi-C pretext files (but NOT snapshots - those are optional and may not exist)
    inputs.extend(get_report_hic_inputs(wildcards))
    
    # NOTE: Hi-C snapshots are NOT included here as required inputs
    # because PretextSnapshot can fail. They're passed as params and
    # checked at runtime in the shell command.
    
    # Blobplots if enabled (depends on sentinel .done file)
    inputs.extend(get_blobplot_inputs(wildcards))
    
    # FCS-GX is NOT included here as a required input
    # because it may not exist. It's passed as params and
    # checked at runtime in the shell command.
    
    return inputs


# -------------------------------------------------------------------------------
# RULES
# -------------------------------------------------------------------------------

rule Z00_generate_report:
    """Generate final markdown report aggregating all analysis results."""
    input:
        deps = get_all_report_inputs
    output:
        report = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "{asm_id}_report.md"
        )
    params:
        species = lambda w: w.species,
        asm_id = lambda w: w.asm_id,
        gfastats = lambda w: get_report_gfastats_inputs(w),
        compleasm = lambda w: get_report_compleasm_inputs(w),
        merqury_qv = lambda w: get_report_merqury_inputs(w)['qv'],
        merqury_completeness = lambda w: get_report_merqury_inputs(w)['completeness'],
        genomescope_plot = lambda w: get_report_genomescope_input(w),
        merqury_dir = lambda w: os.path.join(
            config["OUT_FOLDER"], "GEP2_results", w.species, w.asm_id, "merqury"
        ),
        inspector = lambda w: get_report_inspector_inputs(w),
        hic_snapshots = lambda w: get_report_hic_snapshots(w),
        blobplots = lambda w: get_blobplot_dirs(w),
        fcs_gx_dirs = lambda w: get_report_fcs_dirs(w),
        script_path = str(SCRIPTS_DIR / "make_gep2_report.py")
    container: CONTAINERS["gep2_base"]
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 30
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "Z00_generate_report.log"
        )
    shell:
        """
        exec > {log} 2>&1
        
        cmd="python {params.script_path} -s {params.species} -a {params.asm_id} -g {params.gfastats}"
        
        # Compleasm: extract full_table.tsv from tar.gz files
        COMPLEASM_FULLS=""
        COMPLEASM_CLEANUP=""
        for targz in {params.compleasm}; do
            if [ -n "$targz" ]; then
                extract_dir=$(dirname "$targz")
                tar -xzf "$targz" -C "$extract_dir"
                full_table=$(find "$extract_dir" -name "full_table.tsv" -path "*_odb*" | head -1)
                if [ -n "$full_table" ]; then
                    COMPLEASM_FULLS="$COMPLEASM_FULLS $full_table"
                    COMPLEASM_CLEANUP="$COMPLEASM_CLEANUP $(dirname "$full_table")"
                fi
            fi
        done
        
        if [ -n "$COMPLEASM_FULLS" ]; then
            cmd="$cmd --compleasm-full $COMPLEASM_FULLS"
        fi
        
        if [ -n "{params.merqury_qv}" ]; then
            cmd="$cmd -q {params.merqury_qv}"
        fi
        
        if [ -n "{params.merqury_completeness}" ]; then
            cmd="$cmd -m {params.merqury_completeness}"
        fi
        
        if [ -n "{params.genomescope_plot}" ]; then
            cmd="$cmd --genomescope-plot {params.genomescope_plot}"
        fi
        
        if [ -n "{params.merqury_qv}" ]; then
            cmd="$cmd --merqury-dir {params.merqury_dir}"
        fi
        
        if [ -n "{params.inspector}" ]; then
            cmd="$cmd --Inspector {params.inspector}"
        fi
        
        # Hi-C snapshots are optional - only add files that actually exist
        HIC_SNAPSHOTS=""
        for snapshot in {params.hic_snapshots}; do
            if [ -f "$snapshot" ]; then
                HIC_SNAPSHOTS="$HIC_SNAPSHOTS $snapshot"
            else
                echo "[GEP2] Hi-C snapshot not found (skipping): $snapshot"
            fi
        done
        
        if [ -n "$HIC_SNAPSHOTS" ]; then
            cmd="$cmd --hic $HIC_SNAPSHOTS"
        fi
        

        BLOB_PNGS=""
        for blobdir in {params.blobplots}; do
            if [ -d "$blobdir" ]; then
                # Blob plot: BlobDir.blob.circle.png
                BLOB_PNG=$(find "$blobdir" -maxdepth 1 -name "*.blob.circle.png" -type f | head -1)
                if [ -n "$BLOB_PNG" ]; then
                    BLOB_PNGS="$BLOB_PNGS $BLOB_PNG"
                fi
                # Could also grab snail/cumulative here if the script supports them
            fi
        done

        if [ -n "$BLOB_PNGS" ]; then
            cmd="$cmd --blob $BLOB_PNGS"
        fi


        # FCS-GX: find report files at runtime (filenames contain unpredictable taxid)
        FCS_GX_FILES=""
        for fcs_dir in {params.fcs_gx_dirs}; do
            report=$(find "$fcs_dir" -name "*fcs_gx_report.txt" | head -1)
            if [ -n "$report" ]; then
                FCS_GX_FILES="$FCS_GX_FILES $report"
            fi
        done
        
        if [ -n "$FCS_GX_FILES" ]; then
            cmd="$cmd --fcs-gx $FCS_GX_FILES"
        fi
        
        cmd="$cmd -o {output.report}"
        
        echo "[GEP2] Command: $cmd"
        $cmd
        
        # Cleanup extracted compleasm directories
        for cleanup_dir in $COMPLEASM_CLEANUP; do
            if [ -d "$cleanup_dir" ]; then
                rm -rf "$cleanup_dir"
                echo "[GEP2] Cleaned up extracted compleasm: $cleanup_dir"
            fi
        done
        """