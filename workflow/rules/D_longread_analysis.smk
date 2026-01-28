# -------------------------------------------------------------------------------
# GEP2 - Long Read Analysis Rules
# -------------------------------------------------------------------------------

# Note: The following are defined in the main Snakefile:
#   - samples_config: Parsed sample configuration
#   - get_assembly_files(): Get assembly files for species/assembly
#   - get_assembly_input(): Get primary assembly file
#   - normalize_read_type(): Normalize read type names
#   - _should_skip_analysis(): Check if analysis should be skipped
#   - _as_bool(): Convert config value to boolean


# -------------------------------------------------------------------------------
# INPUT FUNCTIONS
# -------------------------------------------------------------------------------

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
    """Get specific assembly file for Inspector based on asm_basename."""
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
    
    # Find the assembly file matching asm_basename
    asm_files = get_assembly_files(wildcards.species, wildcards.asm_id)
    
    for asm_key, asm_path in asm_files.items():
        if not asm_path or asm_path == "None":
            continue
        
        if get_assembly_basename(asm_path) == wildcards.asm_basename:
            # Check if it needs to be downloaded
            if is_ncbi_assembly_accession(asm_path) or is_url(asm_path):
                # Downloaded assembly path
                return os.path.join(
                    config["OUT_FOLDER"], "GEP2_results", "downloaded_data",
                    wildcards.species, "assemblies", 
                    f"{wildcards.asm_basename}.fna.gz"
                )
            else:
                # Local path
                return asm_path
    
    return []


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


# -------------------------------------------------------------------------------
# RULES
# -------------------------------------------------------------------------------

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
            "inspector", "{asm_basename}", "valid_contig.fa.gz"
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
        
        # Create temp directory
        WORK_DIR="$(gep2_get_workdir 100)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GEP2_inspector_{wildcards.species}_{wildcards.asm_basename}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        
        cd "$TEMP_DIR"
        
        # Run Inspector, output to temp directory
        inspector.py \
            -c {input.asm} \
            -r {params.reads_str} \
            -o "$TEMP_DIR/output" \
            --datatype {params.datatype} \
            -t {threads}
        
        cd output

        # Compress and archive
        pigz -c -p {threads} valid_contig.fa > valid_contig.fa.gz

        tar -cf ae_merge_workspace.tar ae_merge_workspace/
        pigz -p {threads} ae_merge_workspace.tar

        tar -cf base_error_workspace.tar base_error_workspace/
        pigz -p {threads} base_error_workspace.tar

        tar -cf debreak_workspace.tar debreak_workspace/
        pigz -p {threads} debreak_workspace.tar

        tar -cf map_depth.tar map_depth/
        pigz -p {threads} map_depth.tar


        # Verify
        pigz -t valid_contig.fa.gz
        tar -tzf ae_merge_workspace.tar.gz >/dev/null
        tar -tzf base_error_workspace.tar.gz >/dev/null
        tar -tzf debreak_workspace.tar.gz >/dev/null
        tar -tzf map_depth.tar.gz >/dev/null

        # Clean
        rm valid_contig.fa
        rm -r ae_merge_workspace
        rm -r base_error_workspace
        rm -r debreak_workspace
        rm -r map_depth

        cd "$TEMP_DIR"

        # Copy results to final location
        echo "[GEP2] Copying results to {params.outdir}"
        mkdir -p {params.outdir}
        cp -r "$TEMP_DIR/output/"* {params.outdir}/
        
        echo "[GEP2] ✅ Inspector completed successfully"
        """
