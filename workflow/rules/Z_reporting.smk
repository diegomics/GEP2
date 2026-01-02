# ═══════════════════════════════════════════════════════════════════════════════
# GEP2 - Report Generation Rules
# ═══════════════════════════════════════════════════════════════════════════════

# Note: The following are defined in the main Snakefile:
#   - get_assembly_files(): Get assembly files for species/assembly
#   - get_assembly_basename(): Extract basename from filepath
#   - get_priority_read_type(): Get highest priority read type for species
#   - get_kmer_length(): Get k-mer length based on config/read type
#   - _as_bool(): Convert config value to boolean


# ═══════════════════════════════════════════════════════════════════════════════
# INPUT FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

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
    """Get all compleasm output files for this assembly (if enabled)."""
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
                f"{asm_basename}_summary.txt"
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
    
    inspector_dir = os.path.join(
        config["OUT_FOLDER"], "GEP2_results", wildcards.species,
        wildcards.asm_id, "inspector"
    )
    
    return [os.path.join(inspector_dir, "summary_statistics")]


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
    
    return inputs


# ═══════════════════════════════════════════════════════════════════════════════
# RULES
# ═══════════════════════════════════════════════════════════════════════════════

rule Z00_generate_report:
    """Generate final markdown report aggregating all analysis results."""
    input:
        deps = get_all_report_inputs
    output:
        report = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "{asm_id}_report.md"
        ),
        pdf = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "{asm_id}_report.pdf"
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
        script_path = os.path.join(BASEDIR, "scripts", "misc", "make_gep2_report.py")
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
        cmd="python {params.script_path} -s {params.species} -a {params.asm_id} -g {params.gfastats}"
        
        if [ -n "{params.compleasm}" ]; then
            cmd="$cmd -c {params.compleasm}"
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
        
        cmd="$cmd --also-pdf -o {output.report}"
        
        echo "Command: $cmd" > {log}
        $cmd >> {log} 2>&1
        """
