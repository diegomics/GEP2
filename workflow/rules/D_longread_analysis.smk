# ═══════════════════════════════════════════════════════════════════════════════
# GEP2 - Long Read Analysis Rules
# ═══════════════════════════════════════════════════════════════════════════════

# Note: The following are defined in the main Snakefile:
#   - samples_config: Parsed sample configuration
#   - get_assembly_files(): Get assembly files for species/assembly
#   - get_assembly_input(): Get primary assembly file
#   - normalize_read_type(): Normalize read type names
#   - _should_skip_analysis(): Check if analysis should be skipped
#   - _as_bool(): Convert config value to boolean


# ═══════════════════════════════════════════════════════════════════════════════
# INPUT FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def _get_long_read_type_for_assembly(species, asm_id):
    """Get the priority long read type available for this assembly.
    
    Returns 'hifi', 'ont', or None (prefers hifi over ont).
    """
    try:
        asm_data = samples_config["sp_name"][species]["asm_id"][asm_id]
        read_type_dict = asm_data.get("read_type", {})
        
        has_hifi = False
        has_ont = False
        
        for read_type, rt_data in read_type_dict.items():
            if not read_type or read_type == "None" or not rt_data:
                continue
            
            rt_normalized = normalize_read_type(read_type)
            read_files = rt_data.get("read_files", {})
            has_reads = any(v and v != "None" for v in read_files.values())
            
            if has_reads:
                if rt_normalized == "hifi":
                    has_hifi = True
                elif rt_normalized == "ont":
                    has_ont = True
        
        # Prefer hifi over ont
        if has_hifi:
            return "hifi"
        elif has_ont:
            return "ont"
        else:
            return None
            
    except (KeyError, TypeError, AttributeError):
        return None


def _get_long_reads_for_inspector(species, asm_id, read_type):
    """Get all processed long read files for Inspector.
    
    Returns list of processed read file paths.
    """
    reads = []
    
    try:
        asm_data = samples_config["sp_name"][species]["asm_id"][asm_id]
        read_type_dict = asm_data.get("read_type", {})
        
        for rt_key, rt_data in read_type_dict.items():
            if not rt_key or rt_key == "None" or not rt_data:
                continue
            
            rt_normalized = normalize_read_type(rt_key)
            if rt_normalized != read_type:
                continue
            
            read_files = rt_data.get("read_files", {})
            
            for path_key, path_value in sorted(read_files.items()):
                if not path_value or path_value == "None":
                    continue
                
                # Extract base name from path
                if isinstance(path_value, str) and "," in path_value:
                    paths = [p.strip() for p in path_value.split(",")]
                else:
                    paths = [str(path_value)]
                
                for p in paths:
                    basename = os.path.basename(p)
                    # Extract the base identifier
                    base = re.sub(r'^(hifi|ont)_Path\d+_', '', basename, flags=re.IGNORECASE)
                    base = base.replace(".fq.gz", "").replace(".fastq.gz", "")
                    base = base.replace("_filtered", "").replace("_corrected", "")
                    
                    # Find the path index
                    idx = path_key.replace("Path", "") if "Path" in str(path_key) else "1"
                    
                    # Construct processed read path
                    base_dir = os.path.join(
                        config["OUT_FOLDER"], "GEP2_results", "data", species,
                        "reads", read_type
                    )
                    
                    if read_type == "hifi":
                        if config.get("TRIM_HIFI", True):
                            read_path = os.path.join(
                                base_dir, "processed", 
                                f"hifi_Path{idx}_{base}_filtered.fq.gz"
                            )
                        else:
                            read_path = os.path.join(
                                base_dir, f"hifi_Path{idx}_{base}.fq.gz"
                            )
                    elif read_type == "ont":
                        if config.get("CORRECT_ONT", True):
                            read_path = os.path.join(
                                base_dir, "processed",
                                f"ont_Path{idx}_{base}_corrected.fq.gz"
                            )
                        else:
                            read_path = os.path.join(
                                base_dir, f"ont_Path{idx}_{base}.fq.gz"
                            )
                    else:
                        continue
                    
                    if read_path not in reads:
                        reads.append(read_path)
                        
    except (KeyError, TypeError, AttributeError):
        pass
    
    return reads


def get_inspector_asm_input(wildcards):
    """Get assembly file for Inspector."""
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
    
    return get_assembly_input(wildcards)


def get_inspector_reads_input(wildcards):
    """Get long read files for Inspector."""
    # Global toggle
    if not _as_bool(config.get("RUN_INSP", True)):
        return []
    
    # Per-assembly skip
    if _should_skip_analysis(wildcards.species, wildcards.asm_id, "insp"):
        return []
    
    # Get priority long read type
    long_rt = _get_long_read_type_for_assembly(wildcards.species, wildcards.asm_id)
    if not long_rt:
        return []
    
    return _get_long_reads_for_inspector(wildcards.species, wildcards.asm_id, long_rt)


def get_inspector_datatype(wildcards):
    """Get Inspector datatype parameter based on read type."""
    long_rt = _get_long_read_type_for_assembly(wildcards.species, wildcards.asm_id)
    if long_rt == "hifi":
        return "hifi"
    elif long_rt == "ont":
        return "nanopore"
    else:
        return "hifi"  # fallback


# ═══════════════════════════════════════════════════════════════════════════════
# RULES
# ═══════════════════════════════════════════════════════════════════════════════

rule D01_run_inspector:
    """Run Inspector for assembly evaluation using long reads."""
    input:
        asm = get_inspector_asm_input,
        reads = get_inspector_reads_input
    output:
        summary = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "inspector", "{asm_basename}", "summary_statistics"
        ),
        valid_contig = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "inspector", "{asm_basename}", "valid_contig.fa"
        )
    params:
        outdir = lambda w: os.path.join(
            config["OUT_FOLDER"], "GEP2_results", w.species, w.asm_id, 
            "inspector", w.asm_basename
        ),
        datatype = get_inspector_datatype,
        reads_str = lambda w: " ".join(get_inspector_reads_input(w))
    threads: cpu_func("inspector")
    resources:
        mem_mb = mem_func("inspector"),
        runtime = time_func("inspector")
    container: CONTAINERS["inspector"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "D00_inspector_{asm_basename}.log"
        )
    benchmark:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "D00_inspector_{asm_basename}_benchmark.txt"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Running Inspector for {wildcards.species}/{wildcards.asm_id}/{wildcards.asm_basename}"
        echo "[GEP2] Assembly: {input.asm}"
        echo "[GEP2] Reads: {params.reads_str}"
        echo "[GEP2] Datatype: {params.datatype}"
        
        mkdir -p {params.outdir}
        
        # Create temp directory
        TEMP_DIR="$(mktemp -d "$GEP2_TMP/GEP2_inspector_{wildcards.species}_{wildcards.asm_basename}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        
        cd "$TEMP_DIR"
        
        inspector.py \
            -c {input.asm} \
            -r {params.reads_str} \
            -o {params.outdir} \
            --datatype {params.datatype} \
            -t {threads}
        
        echo "[GEP2] ✅ Inspector completed successfully"
        echo ""
        echo "=== Summary Statistics ==="
        cat {output.summary}
        """
