# -------------------------------------------------------------------------------
# GEP2 - Hi-C Analysis Rules
# -------------------------------------------------------------------------------
#
# Workflow:
#   1. chromap index      → index assembly
#   2. chromap map        → .pairs file
#   3. pairtools stats    → stats from .pairs
#   4. genome file        → scaffold sizes for cooler
#   5. cooler cload pairs → .cool file
#   6. cooler zoomify     → .mcool file
#   7. pairtools split    → .bam from .pairs
#   8. PretextMap         → .pretext from .bam
#
# Note: The following are defined in the main Snakefile:
#   - samples_config: Parsed sample configuration
#   - get_assembly_files(): Get assembly files for species/assembly
#   - get_assembly_input(): Get primary assembly file
#   - get_assembly_basename(): Extract basename from filepath
#   - normalize_read_type(): Normalize read type names
#   - _should_skip_analysis(): Check if analysis should be skipped
#   - _as_bool(): Convert config value to boolean


# -------------------------------------------------------------------------------
# INPUT FUNCTIONS
# -------------------------------------------------------------------------------

def _has_hic_reads_for_assembly(species, asm_id):
    """Check if assembly has Hi-C reads available."""
    try:
        asm_data = samples_config["sp_name"][species]["asm_id"][asm_id]
        read_type_dict = asm_data.get("read_type", {})
        
        for read_type, rt_data in read_type_dict.items():
            if not read_type or read_type == "None" or not rt_data:
                continue
            
            rt_normalized = normalize_read_type(read_type)
            if rt_normalized != "hic":
                continue
            
            read_files = rt_data.get("read_files", {})
            if any(v and v != "None" for v in read_files.values()):
                return True
        
        return False
        
    except (KeyError, TypeError, AttributeError):
        return False


def _get_hic_reads_for_assembly(species, asm_id):
    """Get all Hi-C read files for this assembly.
    
    Returns dict with 'r1' and 'r2' lists of file paths.
    """
    r1_files = []
    r2_files = []
    
    try:
        asm_data = samples_config["sp_name"][species]["asm_id"][asm_id]
        read_type_dict = asm_data.get("read_type", {})
        
        for rt_key, rt_data in read_type_dict.items():
            if not rt_key or rt_key == "None" or not rt_data:
                continue
            
            rt_normalized = normalize_read_type(rt_key)
            if rt_normalized != "hic":
                continue
            
            read_files = rt_data.get("read_files", {})
            
            # Group paths by Path index
            path_groups = {}
            for path_key, path_value in read_files.items():
                if not path_value or path_value == "None":
                    continue
                
                # Extract path index
                idx = "1"
                if "Path" in str(path_key):
                    idx = str(path_key).replace("Path", "")
                
                if idx not in path_groups:
                    path_groups[idx] = []
                
                # Handle comma-separated paths
                if isinstance(path_value, str) and "," in path_value:
                    paths = [p.strip() for p in path_value.split(",")]
                else:
                    paths = [str(path_value)]
                
                path_groups[idx].extend(paths)
            
            # Process each path group
            for idx, paths in sorted(path_groups.items()):
                # Extract base names and construct processed paths
                seen_bases = set()
                for p in paths:
                    basename = os.path.basename(p)
                    base = re.sub(r'^(hic)_Path\d+_', '', basename, flags=re.IGNORECASE)
                    base = base.replace(".fq.gz", "").replace(".fastq.gz", "")
                    base = base.replace("_1", "").replace("_2", "")
                    
                    if base in seen_bases:
                        continue
                    seen_bases.add(base)
                    
                    base_dir = os.path.join(
                        config["OUT_FOLDER"], "GEP2_results", "data", species,
                        "reads", "hic"
                    )
                    
                    # Hi-C reads are paired-end
                    r1_path = os.path.join(base_dir, f"hic_Path{idx}_{base}_1.fq.gz")
                    r2_path = os.path.join(base_dir, f"hic_Path{idx}_{base}_2.fq.gz")
                    
                    if r1_path not in r1_files:
                        r1_files.append(r1_path)
                    if r2_path not in r2_files:
                        r2_files.append(r2_path)
                        
    except (KeyError, TypeError, AttributeError):
        pass
    
    return {"r1": r1_files, "r2": r2_files}


def _should_run_hic(species, asm_id):
    """Check if Hi-C analysis should run for this assembly."""
    # Global toggle
    if not _as_bool(config.get("RUN_HIC", True)):
        return False
    
    # Per-assembly skip
    if _should_skip_analysis(species, asm_id, "hic"):
        return False
    
    # Check for Hi-C reads
    if not _has_hic_reads_for_assembly(species, asm_id):
        return False
    
    return True


def get_hic_asm_input(wildcards):
    """Get specific assembly file for Hi-C analysis based on asm_basename."""
    if not _should_run_hic(wildcards.species, wildcards.asm_id):
        return []
    
    asm_files = get_assembly_files(wildcards.species, wildcards.asm_id)
    
    for asm_key, asm_path in asm_files.items():
        if not asm_path or asm_path == "None":
            continue
        if get_assembly_basename(asm_path) == wildcards.asm_basename:
            # Return the actual input path (could be downloaded or local)
            return get_assembly_input_for_basename(wildcards, wildcards.asm_basename)
    
    return []


def get_assembly_input_for_basename(wildcards, asm_basename):
    """Get assembly input path for a specific basename."""
    asm_files = get_assembly_files(wildcards.species, wildcards.asm_id)
    
    for asm_key, asm_path in asm_files.items():
        if not asm_path or asm_path == "None":
            continue
        if get_assembly_basename(asm_path) == asm_basename:
            # Check if it needs to be downloaded
            if asm_path.startswith("GCA_") or asm_path.startswith("GCF_") or asm_path.startswith("http"):
                # Downloaded assembly
                return os.path.join(
                    config["OUT_FOLDER"], "GEP2_results", "downloaded_data",
                    wildcards.species, "assemblies", f"{asm_basename}.fna.gz"
                )
            else:
                return asm_path
    
    return []


def get_hic_reads_r1(wildcards):
    """Get Hi-C R1 read files."""
    if not _should_run_hic(wildcards.species, wildcards.asm_id):
        return []
    
    hic_reads = _get_hic_reads_for_assembly(wildcards.species, wildcards.asm_id)
    return hic_reads["r1"]


def get_hic_reads_r2(wildcards):
    """Get Hi-C R2 read files."""
    if not _should_run_hic(wildcards.species, wildcards.asm_id):
        return []
    
    hic_reads = _get_hic_reads_for_assembly(wildcards.species, wildcards.asm_id)
    return hic_reads["r2"]

def _get_track_bedgraphs(wildcards):
    """Get list of bedgraph files to add as tracks based on config."""
    base_dir = os.path.join(
        config["OUT_FOLDER"], "GEP2_results", wildcards.species, wildcards.asm_id,
        "hic", wildcards.asm_basename, "tracks"
    )
    
    tracks = []
    
    if _as_bool(config.get("GAP_TRACK", False)):
        tracks.append(os.path.join(base_dir, "gap_density.bedgraph"))
    
    # TELO_TRACK can be: auto, off/false, or a custom motif string
    telo_setting = str(config.get("TELO_TRACK", "off")).strip().lower()
    if telo_setting not in ["off", "false", "no", "0", ""]:
        tracks.append(os.path.join(base_dir, "telo_density.bedgraph"))
    
    if _as_bool(config.get("COVER_TRACK", False)):
        # Only add if long reads are available
        if _has_long_reads_for_assembly(wildcards.species, wildcards.asm_id):
            tracks.append(os.path.join(base_dir, "coverage.bedgraph"))
    
    if _as_bool(config.get("REP_TRACK", False)):
        tracks.append(os.path.join(base_dir, "sdust_density.bedgraph"))
        tracks.append(os.path.join(base_dir, "ldust_density.bedgraph"))
    
    return tracks


def _has_long_reads_for_assembly(species, asm_id):
    """Check if assembly has long reads (HiFi or ONT) available for coverage track."""
    try:
        asm_data = samples_config["sp_name"][species]["asm_id"][asm_id]
        read_type_dict = asm_data.get("read_type", {})
        
        for read_type, rt_data in read_type_dict.items():
            if not read_type or read_type == "None" or not rt_data:
                continue
            
            rt_normalized = normalize_read_type(read_type)
            if rt_normalized in ["hifi", "ont"]:
                read_files = rt_data.get("read_files", {})
                if any(v and v != "None" for v in read_files.values()):
                    return True
        
        return False
        
    except (KeyError, TypeError, AttributeError):
        return False


def _get_long_reads_for_coverage(wildcards):
    """Get long read files for coverage track (prefers HiFi over ONT).
    
    Uses processed reads if available/enabled, otherwise raw reads.
    """
    try:
        asm_data = samples_config["sp_name"][wildcards.species]["asm_id"][wildcards.asm_id]
        read_type_dict = asm_data.get("read_type", {})
        
        # Master switch for read processing
        reads_proc_enabled = _as_bool(config.get("READS_PROC", False))
        
        # Prefer HiFi, fall back to ONT
        for preferred_type in ["hifi", "ont"]:
            for read_type, rt_data in read_type_dict.items():
                if not read_type or read_type == "None" or not rt_data:
                    continue
                
                rt_normalized = normalize_read_type(read_type)
                if rt_normalized != preferred_type:
                    continue
                
                read_files = rt_data.get("read_files", {})
                files = []
                
                # Determine if we should use processed or raw reads
                if rt_normalized == "hifi":
                    use_processed = reads_proc_enabled and _as_bool(config.get("TRIM_HIFI", False))
                    proc_suffix = "_filtered.fq.gz"
                elif rt_normalized == "ont":
                    use_processed = reads_proc_enabled and _as_bool(config.get("CORRECT_ONT", False))
                    proc_suffix = "_corrected.fq.gz"
                else:
                    use_processed = False
                    proc_suffix = "_proc.fq.gz"
                
                for path_key, path_value in read_files.items():
                    if not path_value or path_value == "None":
                        continue
                    
                    # Get path index
                    idx = "1"
                    if "Path" in str(path_key):
                        idx = str(path_key).replace("Path", "")
                    
                    # Handle comma-separated paths (multiple files per path)
                    if isinstance(path_value, str):
                        paths = [p.strip() for p in path_value.split(",") if p.strip()]
                    else:
                        paths = [str(path_value)]
                    
                    for p in paths:
                        basename = os.path.basename(str(p))
                        # Remove extension and path prefix to get accession
                        base = re.sub(r'^(' + rt_normalized + r')_Path\d+_', '', basename, flags=re.IGNORECASE)
                        base = base.replace(".fq.gz", "").replace(".fastq.gz", "")
                        base = re.sub(r'_(1|2)$', '', base)  # Remove _1 or _2 suffix if present
                        
                        base_dir = os.path.join(
                            config["OUT_FOLDER"], "GEP2_results", "data", wildcards.species,
                            "reads", rt_normalized
                        )
                        
                        if use_processed:
                            # Processed reads path
                            proc_path = os.path.join(
                                base_dir, "processed",
                                f"{rt_normalized}_Path{idx}_{base}{proc_suffix}"
                            )
                            files.append(proc_path)
                        else:
                            # Raw reads path (symlink in data folder)
                            raw_path = os.path.join(
                                base_dir,
                                f"{rt_normalized}_Path{idx}_{base}.fq.gz"
                            )
                            files.append(raw_path)
                
                if files:
                    return files
        
        return []
        
    except (KeyError, TypeError, AttributeError):
        return []


def _get_long_read_type_for_coverage(wildcards):
    """Get the type of long reads being used for coverage (hifi or ont)."""
    try:
        asm_data = samples_config["sp_name"][wildcards.species]["asm_id"][wildcards.asm_id]
        read_type_dict = asm_data.get("read_type", {})
        
        # Prefer HiFi, fall back to ONT (same order as _get_long_reads_for_coverage)
        for preferred_type in ["hifi", "ont"]:
            for read_type, rt_data in read_type_dict.items():
                if not read_type or read_type == "None" or not rt_data:
                    continue
                
                rt_normalized = normalize_read_type(read_type)
                if rt_normalized != preferred_type:
                    continue
                
                read_files = rt_data.get("read_files", {})
                if any(v and v != "None" for v in read_files.values()):
                    return rt_normalized
        
        return "hifi"  # Default fallback
        
    except (KeyError, TypeError, AttributeError):
        return "hifi"


# -------------------------------------------------------------------------------
# RULES
# -------------------------------------------------------------------------------

rule E00_chromap_index:
    """Create chromap index for the assembly."""
    input:
        asm = get_hic_asm_input
    output:
        index = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "chromap_index.index"
        )
    params:
        outdir = lambda w: os.path.join(
            config["OUT_FOLDER"], "GEP2_results", w.species, w.asm_id, 
            "hic", w.asm_basename
        )
    threads: cpu_func("chromap_index")
    resources:
        mem_mb = mem_func("chromap_index"),
        runtime = time_func("chromap_index")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E00_chromap_index_{asm_basename}.log"
        )
    benchmark:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E00_chromap_index_{asm_basename}_benchmark.txt"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Creating chromap index for {wildcards.species}/{wildcards.asm_id}/{wildcards.asm_basename}"
        echo "[GEP2] Assembly: {input.asm}"
        
        mkdir -p {params.outdir}
        
        chromap -i -r {input.asm} -o {output.index}
        
        echo "[GEP2] ✅ Chromap index created successfully"
        """


rule E01_chromap_map:
    """Map Hi-C reads to assembly with chromap, output .pairs file."""
    input:
        asm = get_hic_asm_input,
        index = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "chromap_index.index"
        ),
        r1 = get_hic_reads_r1,
        r2 = get_hic_reads_r2
    output:
        pairs = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}.pairs.gz"
        )
    params:
        outdir = lambda w: os.path.join(
            config["OUT_FOLDER"], "GEP2_results", w.species, w.asm_id, 
            "hic", w.asm_basename
        ),
        r1_str = lambda w: ",".join(get_hic_reads_r1(w)),
        r2_str = lambda w: ",".join(get_hic_reads_r2(w)),
        mapq = config.get("HIC_MAPQ", 0)
    threads: cpu_func("chromap_map")
    resources:
        mem_mb = mem_func("chromap_map"),
        runtime = time_func("chromap_map")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E01_chromap_map_{asm_basename}.log"
        )
    benchmark:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E01_chromap_map_{asm_basename}_benchmark.txt"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Mapping Hi-C reads for {wildcards.species}/{wildcards.asm_id}/{wildcards.asm_basename}"
        echo "[GEP2] Assembly: {input.asm}"
        echo "[GEP2] R1 files: {params.r1_str}"
        echo "[GEP2] R2 files: {params.r2_str}"
        echo "[GEP2] MAPQ threshold: {params.mapq}"
        
        mkdir -p {params.outdir}
        
        # Create temp directory
        WORK_DIR="$(gep2_get_workdir 100)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GEP2_chromap_{wildcards.species}_{wildcards.asm_basename}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        
        # Map with chromap, output pairs format
        chromap --preset hic \
            -q {params.mapq} \
            -x {input.index} \
            -r {input.asm} \
            -1 {params.r1_str} \
            -2 {params.r2_str} \
            -t {threads} \
            --pairs \
            -o "$TEMP_DIR/{wildcards.asm_basename}.pairs"
        
        echo "[GEP2] Compressing pairs file..."
        bgzip -@ {threads} -c "$TEMP_DIR/{wildcards.asm_basename}.pairs" > {output.pairs}
        
        echo "[GEP2] ✅ Chromap mapping completed successfully"
        """


rule E02_pairtools_stats:
    """Generate statistics from .pairs file."""
    input:
        pairs = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}.pairs.gz"
        )
    output:
        stats = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}.pairtools_stats.txt"
        )
    threads: cpu_func("light_task")
    resources:
        mem_mb = mem_func("light_task"),
        runtime = time_func("light_task")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E02_pairtools_stats_{asm_basename}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Running pairtools stats for {wildcards.asm_basename}"
        
        pairtools stats {input.pairs} -o {output.stats}
        
        echo "[GEP2] ✅ Pairtools stats completed"
        """


rule E03_create_genome_file:
    """Create genome file (chrom sizes) for cooler."""
    input:
        asm = get_hic_asm_input
    output:
        genome = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}.genome"
        )
    threads: cpu_func("light_task")
    resources:
        mem_mb = mem_func("light_task"),
        runtime = time_func("light_task")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E03_create_genome_file_{asm_basename}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Creating genome file for {wildcards.asm_basename}"
        
        mkdir -p $(dirname {output.genome})
        
        # Use seqkit fx2tab with -l for length, then extract columns 1 and 4
        # Output format: name \t seq \t ... \t length (when using -l)
        seqkit fx2tab -l {input.asm} | \
            awk 'BEGIN{{OFS="\\t"}} {{print $1, $NF}}' > {output.genome}
        
        echo "[GEP2] ✅ Genome file created"
        echo ""
        echo "=== Scaffold sizes (first 20) ==="
        head -20 {output.genome}
        
        # Sanity check - verify no zero-length scaffolds
        ZEROS=$(awk '$2 == 0' {output.genome} | wc -l)
        if [ "$ZEROS" -gt 0 ]; then
            echo "[GEP2] ⚠️  WARNING: Found $ZEROS scaffolds with zero length!"
        fi
        """


rule E04_cooler_cload:
    """Create .cool file from .pairs using cooler cload pairs."""
    input:
        pairs = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}.pairs.gz"
        ),
        genome = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}.genome"
        )
    output:
        cool = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}.cool"
        )
    params:
        binsize = config.get("HIC_BINSIZE", 10000)
    threads: cpu_func("light_task")
    resources:
        mem_mb = mem_func("pretext_map"),
        runtime = time_func("pretext_map")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E04_cooler_cload_{asm_basename}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Creating .cool file for {wildcards.asm_basename}"
        echo "[GEP2] Bin size: {params.binsize}"
        
        cooler cload pairs \
            -c1 2 -p1 3 -c2 4 -p2 5 \
            {input.genome}:{params.binsize} \
            {input.pairs} \
            {output.cool}
        
        echo "[GEP2] ✅ Cool file created"
        """


rule E05_cooler_zoomify:
    """Create multi-resolution .mcool file from .cool."""
    input:
        cool = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}.cool"
        )
    output:
        mcool = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}.mcool"
        )
    threads: cpu_func("pretext_map")
    resources:
        mem_mb = mem_func("pretext_map"),
        runtime = time_func("pretext_map")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E05_cooler_zoomify_{asm_basename}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Creating multi-resolution .mcool for {wildcards.asm_basename}"
        
        cooler zoomify \
            -p {threads} \
            {input.cool} \
            -o {output.mcool}
        
        echo "[GEP2] ✅ Mcool file created"
        """

rule E06_pretext_map:
    """Create PretextMap using the 'text hack' method, optionally create snapshot."""
    input:
        pairs = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}.pairs.gz"
        )
    output:
        pretext = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}.pretext"
        )
    params:
        highres_flag = "--highRes" if _as_bool(config.get("HIC_HIGH_RES", False)) else "",
        run_snapshot = not _as_bool(config.get("HIC_HIGH_RES", False)),
        outdir = lambda w: os.path.join(
            config["OUT_FOLDER"], "GEP2_results", w.species, w.asm_id,
            "hic", w.asm_basename
        )
    threads: cpu_func("pretext_map")
    resources:
        mem_mb = mem_func("pretext_map"),
        runtime = time_func("pretext_map")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E06_pretext_map_{asm_basename}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Patching pairs header and creating PretextMap..."
        echo "[GEP2] High resolution mode: {params.highres_flag}"
        
        zcat {input.pairs} \
        | sed -e 's/format v1.0.0/format v1.0/' \
              -e '/^#shape/d' \
              -e 's/ pair_type//' \
        | awk 'BEGIN{{OFS=sprintf("%c", 9)}} /^#/ {{print; next}} {{print $1,$2,$3,$4,$5,$6,$7}}' \
        | PretextMap \
            --sortby length \
            --sortorder descend \
            {params.highres_flag} \
            -o {output.pretext}
        
        echo "[GEP2] ✅ PretextMap created successfully"
        ls -lh {output.pretext}
        
        # Create snapshot (only if not in high-res mode)
        if [ "{params.run_snapshot}" = "True" ]; then
            echo ""
            echo "[GEP2] Creating PretextSnapshot..."
            cd {params.outdir}
            
            if PretextSnapshot -m {output.pretext} -r 3840 --sequences "=full" 2>&1; then
                echo "[GEP2] ✅ PretextSnapshot created"
                ls -la {wildcards.asm_basename}_snapshots/ 2>/dev/null || true
            else
                echo "[GEP2] ⚠️  PretextSnapshot failed (non-fatal)"
            fi
        else
            echo ""
            echo "[GEP2] ⚠️  Skipping PretextSnapshot (not compatible with high-res mode)"
        fi
        """

rule E07_windows_bed:
    """Generate 1kb windows BED file for track generation."""
    input:
        asm = get_hic_asm_input,
        genome = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}.genome"
        )
    output:
        windows = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "tracks", "windows_1kb.bed"
        )
    threads: cpu_func("light_task")
    resources:
        mem_mb = mem_func("light_task"),
        runtime = time_func("light_task")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E07_windows_bed_{asm_basename}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Generating 1kb windows for {wildcards.asm_basename}"
        
        mkdir -p $(dirname {output.windows})
        
        # Create windows using pre-generated genome file
        bedtools makewindows -g {input.genome} -w 1000 > {output.windows}
        
        echo "[GEP2] ✅ Windows BED created"
        wc -l {output.windows}
        """


rule E08_gap_track:
    """Generate gap/N density track in 1kb windows."""
    input:
        asm = get_hic_asm_input
    output:
        bedgraph = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "tracks", "gap_density.bedgraph"
        )
    threads: cpu_func("track_light")
    resources:
        mem_mb = mem_func("track_light"),
        runtime = time_func("track_light")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E08_gap_track_{asm_basename}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Generating gap density track for {wildcards.asm_basename}"
        
        mkdir -p $(dirname {output.bedgraph})
        
        # Use seqkit sliding windows and count Ns directly
        seqkit sliding -W 1000 -s 1000 {input.asm} | \
        seqkit fx2tab | \
        awk 'BEGIN{{OFS="\\t"}} 
             {{
                 # Parse header like "scaffold_1_sliding:1-1000"
                 split($1, a, "_sliding:")
                 split(a[2], b, "-")
                 # Count N/n characters in sequence
                 seq = $2
                 gsub(/[^Nn]/, "", seq)
                 n_count = length(seq)
                 printf "%s\\t%d\\t%d\\t%d\\n", a[1], b[1]-1, b[2], n_count
             }}' > {output.bedgraph}
        
        echo "[GEP2] ✅ Gap density track created"
        head -5 {output.bedgraph}
        """


rule E09_sdust_track:
    """Generate low complexity track using sdust."""
    input:
        asm = get_hic_asm_input,
        windows = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "tracks", "windows_1kb.bed"
        )
    output:
        bedgraph = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "tracks", "sdust_density.bedgraph"
        )
    threads: cpu_func("track_heavy")
    resources:
        mem_mb = mem_func("track_heavy"),
        runtime = time_func("track_heavy")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E09_sdust_track_{asm_basename}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Generating sdust low-complexity track for {wildcards.asm_basename}"
        
        # Run sdust and calculate coverage per window
        sdust {input.asm} | \
        bedtools coverage -a {input.windows} -b stdin | \
        awk 'BEGIN{{OFS="\\t"}} 
             {{print $1, $2, $3, $5}}' > {output.bedgraph}
        
        echo "[GEP2] ✅ Sdust track created"
        head -5 {output.bedgraph}
        """


rule E10_ldust_track:
    """Generate long dust (tandem repeats) track."""
    input:
        asm = get_hic_asm_input,
        windows = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "tracks", "windows_1kb.bed"
        )
    output:
        bedgraph = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "tracks", "ldust_density.bedgraph"
        )
    threads: cpu_func("track_heavy")
    resources:
        mem_mb = mem_func("track_heavy"),
        runtime = time_func("track_heavy")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E10_ldust_track_{asm_basename}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Generating longdust track for {wildcards.asm_basename}"
        
        WORK_DIR="$(gep2_get_workdir 10)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GEP2_ldust_{wildcards.asm_basename}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        
        # Run longdust (outputs BED-like format)
        longdust {input.asm} > "$TEMP_DIR/raw_longdust.bed"
        
        # Bin into windows
        bedtools coverage -a {input.windows} -b "$TEMP_DIR/raw_longdust.bed" | \
        awk 'BEGIN{{OFS="\\t"}} 
             {{print $1, $2, $3, $5}}' > {output.bedgraph}
        
        echo "[GEP2] ✅ Longdust track created"
        head -5 {output.bedgraph}
        """


rule E11_coverage_track:
    """Generate long read coverage track."""
    input:
        asm = get_hic_asm_input,
        reads = _get_long_reads_for_coverage
    output:
        bedgraph = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "tracks", "coverage.bedgraph"
        )
    params:
        mm2_preset = lambda w: "map-hifi" if _get_long_read_type_for_coverage(w) == "hifi" else "map-ont"
    threads: cpu_func("minimap2")
    resources:
        mem_mb = mem_func("minimap2"),
        runtime = time_func("minimap2")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E11_coverage_track_{asm_basename}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Generating coverage track for {wildcards.asm_basename}"
        echo "[GEP2] Input reads: {input.reads}"
        echo "[GEP2] Minimap2 preset: {params.mm2_preset}"
        
        WORK_DIR="$(gep2_get_workdir 150)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GEP2_coverage_{wildcards.asm_basename}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        
        mkdir -p $(dirname {output.bedgraph})
        
        # Concatenate all read files if multiple
        READS_INPUT="{input.reads}"
        
        echo "[GEP2] Mapping reads with minimap2 (-ax {params.mm2_preset})..."
        minimap2 -ax {params.mm2_preset} -t {threads} {input.asm} $READS_INPUT | \
            samtools sort -@ {threads} -o "$TEMP_DIR/sorted.bam" -
        
        samtools index -@ {threads} "$TEMP_DIR/sorted.bam"
        
        echo "[GEP2] Calculating coverage with sambamba..."
        sambamba depth window -w 1000 -t {threads} "$TEMP_DIR/sorted.bam" | \
        awk 'BEGIN{{OFS="\\t"}}
            !/^#/ && NR>1 {{printf "%s\\t%s\\t%s\\t%d\\n", $1, $2, $3, $5}}' > {output.bedgraph}
        
        echo "[GEP2] ✅ Coverage track created"
        head -5 {output.bedgraph}
        """


rule E12_telo_track:
    """Generate telomere density track using tidk."""
    input:
        asm = get_hic_asm_input
    output:
        bedgraph = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "tracks", "telo_density.bedgraph"
        )
    params:
        species_name = lambda w: w.species.replace("_", " "),
        telo_setting = config.get("TELO_TRACK", "auto")
    threads: cpu_func("track_light")
    resources:
        mem_mb = mem_func("track_light"),
        runtime = time_func("track_light")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E12_telo_track_{asm_basename}.log"
        )
    shell:
        """
        exec > {log} 2>&1
        
        echo "[GEP2] Generating telomere track for {wildcards.asm_basename}"
        echo "[GEP2] Species: {params.species_name}"
        echo "[GEP2] TELO_TRACK setting: {params.telo_setting}"
        
        WORK_DIR="$(gep2_get_workdir 10)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GEP2_telo_{wildcards.asm_basename}_XXXXXX")"
        
        cleanup() {{
            rm -rf "$TEMP_DIR"
        }}
        trap cleanup EXIT
        
        mkdir -p $(dirname {output.bedgraph})
        
        # Decompress assembly if needed (tidk can't read gzip)
        if [[ "{input.asm}" == *.gz ]]; then
            echo "[GEP2] Decompressing assembly for tidk..."
            gunzip -c {input.asm} > "$TEMP_DIR/assembly.fa"
            ASM_FILE="$TEMP_DIR/assembly.fa"
        else
            ASM_FILE="{input.asm}"
        fi
        
        TRACK_CREATED=false
        TELO_SETTING="{params.telo_setting}"
        TELO_LOWER=$(echo "$TELO_SETTING" | tr '[:upper:]' '[:lower:]')
        
        # Check if user provided a custom motif (not "auto")
        if [ "$TELO_LOWER" != "auto" ]; then
            # User provided a custom telomere motif
            CUSTOM_MOTIF=$(echo "$TELO_SETTING" | tr '[:lower:]' '[:upper:]')
            echo "[GEP2] Using custom telomere motif: $CUSTOM_MOTIF"
            
            # tidk may return non-zero even on success, so use || true
            tidk search \
                --string "$CUSTOM_MOTIF" \
                --window 1000 \
                --output telomere \
                --dir "$TEMP_DIR" \
                "$ASM_FILE" 2>&1 || true
            
            echo "[GEP2] tidk search completed, looking for output files..."
            ls -la "$TEMP_DIR"/ 2>/dev/null || true
            
            RAW_TSV=$(ls "$TEMP_DIR"/*windows*.tsv 2>/dev/null | head -1) || true
            
            if [ -n "$RAW_TSV" ] && [ -s "$RAW_TSV" ]; then
                echo "[GEP2] Found TSV file: $RAW_TSV"
                head -5 "$RAW_TSV"
                
                awk 'BEGIN{{OFS="\\t"}}
                     NR>1 && NF>=4 {{
                         chrom=$1
                         window=$2
                         start=(window-1)*1000
                         end=window*1000
                         total=int($3)+int($4)
                         printf "%s\\t%d\\t%d\\t%d\\n", chrom, start, end, total
                     }}' "$RAW_TSV" > {output.bedgraph}
                TRACK_CREATED=true
                echo "[GEP2] ✅ Telomere track created with custom motif $CUSTOM_MOTIF"
            fi
        else
            # Auto mode: try GoaT first, then tidk explore
            echo "[GEP2] Auto mode: detecting telomere motif..."
            
            # Initialize tidk database (may fail if no network, that's ok)
            tidk build 2>/dev/null || true
            
            # Query GoaT for taxonomic lineage and find matching TIDK clade
            echo "[GEP2] Querying GoaT for taxonomic lineage..."
            CLADE=""
            CLADE=$(python3 << 'PYEOF' || echo ""
import requests
import sys

TIDK_CLADES = {{
    "Crassiclitellata", "Hirudinida", "Phyllodocida", "Eucoccidiorida", 
    "Coleoptera", "Hemiptera", "Hymenoptera", "Lepidoptera", "Odonata", 
    "Orthoptera", "Plecoptera", "Symphypleona", "Trichoptera", 
    "Cheilostomatida", "Chlamydomonadales", "Accipitriformes", "Anura", 
    "Aplousobranchia", "Caprimulgiformes", "Carangiformes", "Carcharhiniformes", 
    "Carnivora", "Chiroptera", "Cypriniformes", "Labriformes", "Perciformes", 
    "Phlebobranchia", "Pleuronectiformes", "Rodentia", "Salmoniformes", 
    "Syngnathiformes", "Actiniaria", "Forcipulatida", "Cardiida", "Pectinida", 
    "Trochida", "Venerida", "Heteronemertea", "Apiales", "Asterales", "Buxales", 
    "Caryophyllales", "Fabales", "Fagales", "Hypnales", "Lamiales", "Malpighiales", 
    "Myrtales", "Poales", "Rosales", "Sapindales", "Solanales"
}}

species = "{params.species_name}"
encoded = requests.utils.quote(species)
url = f'https://goat.genomehubs.org/api/v2/search?query=tax_name%28{{encoded}}%29&result=taxon'

try:
    response = requests.get(url, timeout=30)
    data = response.json()
    
    if data.get('results'):
        lineage = {{node['scientific_name'] for node in data['results'][0]['result']['lineage']}}
        matches = lineage.intersection(TIDK_CLADES)
        if matches:
            print(list(matches)[0])
            sys.exit(0)
except Exception as e:
    print(f"GoaT query failed: {{e}}", file=sys.stderr)

sys.exit(1)
PYEOF
)
            
            if [ -n "$CLADE" ]; then
                echo "[GEP2] Found matching TIDK clade: $CLADE"
                
                # tidk may return non-zero even on success, so use || true
                tidk find \
                    --clade "$CLADE" \
                    --window 1000 \
                    --output telomere \
                    --dir "$TEMP_DIR" \
                    "$ASM_FILE" 2>&1 || true
                
                echo "[GEP2] tidk find completed, checking for output files..."
                ls -la "$TEMP_DIR"/ 2>/dev/null || true
                
                RAW_BG=$(ls "$TEMP_DIR"/*.bedgraph 2>/dev/null | head -1) || true
                
                if [ -n "$RAW_BG" ] && [ -s "$RAW_BG" ]; then
                    echo "[GEP2] Found bedgraph file: $RAW_BG"
                    head -5 "$RAW_BG"
                    
                    awk 'BEGIN{{OFS="\\t"}}
                         !/^#/ && NF>=4 {{printf "%s\\t%s\\t%s\\t%d\\n", $1, $2, $3, int($4)}}' \
                         "$RAW_BG" > {output.bedgraph}
                    TRACK_CREATED=true
                    echo "[GEP2] ✅ Telomere track created with clade $CLADE"
                else
                    echo "[GEP2] ⚠️  No bedgraph file found after tidk find"
                fi
            else
                echo "[GEP2] ⚠️  No matching TIDK clade found in lineage"
            fi
            
            # Fallback: use tidk explore if no clade match or tidk find failed
            if [ "$TRACK_CREATED" = "false" ]; then
                echo "[GEP2] Trying tidk explore to find telomere motif..."
                
                tidk explore \
                    --minimum 5 \
                    --maximum 12 \
                    "$ASM_FILE" > "$TEMP_DIR/explore_output.txt" 2>&1 || true
                
                echo "[GEP2] tidk explore output:"
                head -30 "$TEMP_DIR/explore_output.txt" 2>/dev/null || echo "(no output)"
                
                TOP_MOTIF=$(awk 'NR>1 && /^[ACGT]+/ {{print $1; exit}}' "$TEMP_DIR/explore_output.txt" 2>/dev/null) || true
                
                if [ -n "$TOP_MOTIF" ]; then
                    echo "[GEP2] Found telomere motif via explore: $TOP_MOTIF"
                    
                    tidk search \
                        --string "$TOP_MOTIF" \
                        --window 1000 \
                        --output telomere \
                        --dir "$TEMP_DIR" \
                        "$ASM_FILE" 2>&1 || true
                    
                    echo "[GEP2] tidk search completed, looking for output files..."
                    ls -la "$TEMP_DIR"/ 2>/dev/null || true
                    
                    RAW_TSV=$(ls "$TEMP_DIR"/*windows*.tsv 2>/dev/null | head -1) || true
                    
                    if [ -n "$RAW_TSV" ] && [ -s "$RAW_TSV" ]; then
                        echo "[GEP2] Found TSV file: $RAW_TSV"
                        head -5 "$RAW_TSV"
                        
                        awk 'BEGIN{{OFS="\\t"}}
                             NR>1 && NF>=4 {{
                                 chrom=$1
                                 window=$2
                                 start=(window-1)*1000
                                 end=window*1000
                                 total=int($3)+int($4)
                                 printf "%s\\t%d\\t%d\\t%d\\n", chrom, start, end, total
                             }}' "$RAW_TSV" > {output.bedgraph}
                        TRACK_CREATED=true
                        echo "[GEP2] ✅ Telomere track created with motif $TOP_MOTIF"
                    else
                        echo "[GEP2] ⚠️  No TSV file found after tidk search"
                    fi
                else
                    echo "[GEP2] ⚠️  Could not extract telomere motif from explore output"
                fi
            fi
        fi
        
        # If still no track, create empty file
        if [ "$TRACK_CREATED" = "false" ]; then
            echo "[GEP2] ⚠️  Could not generate telomere track, creating empty file"
            touch {output.bedgraph}
        fi
        
        echo ""
        echo "[GEP2] Output file:"
        ls -la {output.bedgraph}
        echo ""
        echo "[GEP2] First 10 lines:"
        head -10 {output.bedgraph} || echo "(empty)"
        
        # Always exit successfully - empty track is acceptable
        exit 0
        """


rule E13_add_pretext_tracks:
    """Add all enabled tracks to the pretext file."""
    input:
        pretext = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}.pretext"
        ),
        tracks = _get_track_bedgraphs
    output:
        pretext = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}_tracks.pretext"
        )
    threads: cpu_func("track_light")
    resources:
        mem_mb = mem_func("track_light"),
        runtime = time_func("track_light")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E13_add_pretext_tracks_{asm_basename}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Adding tracks to PretextMap for {wildcards.asm_basename}"
        
        # Start with a copy of the base pretext
        cp {input.pretext} {output.pretext}
        
        # Track name mapping (basename -> display name)
        declare -A TRACK_NAMES=(
            ["gap_density"]="gaps"
            ["telo_density"]="telo"
            ["coverage"]="cov"
            ["sdust_density"]="sdust"
            ["ldust_density"]="ldust"
        )
        
        # Add each track sequentially
        TRACKS=({input.tracks})
        
        if [ ${{#TRACKS[@]}} -eq 0 ]; then
            echo "[GEP2] No tracks to add"
            exit 0
        fi
        
        for TRACK_FILE in "${{TRACKS[@]}}"; do
            if [ ! -f "$TRACK_FILE" ]; then
                echo "[GEP2] ⚠️  Track file not found: $TRACK_FILE"
                continue
            fi
            
            # Get track name from filename
            BASENAME=$(basename "$TRACK_FILE" .bedgraph)
            TRACK_NAME="${{TRACK_NAMES[$BASENAME]:-$BASENAME}}"
            
            # Check if bedgraph has data (more than just header)
            DATA_LINES=$(grep -v '^#' "$TRACK_FILE" | wc -l)
            
            if [ "$DATA_LINES" -eq 0 ]; then
                echo "[GEP2] ⚠️  Track $TRACK_NAME has no data, skipping"
                continue
            fi
            
            echo "[GEP2] Adding track: $TRACK_NAME ($DATA_LINES data lines)"
            
            # Create temp file for intermediate output
            TEMP_PRETEXT=$(mktemp)
            
            # Add track using PretextGraph (skip header lines, ensure integers)
            grep -v '^#' "$TRACK_FILE" | \
            awk 'BEGIN{{OFS="\\t"}} NF>=4 {{printf "%s\\t%d\\t%d\\t%d\\n", $1, $2, $3, int($4)}}' | \
            PretextGraph \
                -i {output.pretext} \
                -n "$TRACK_NAME" \
                -o "$TEMP_PRETEXT"
            
            # Replace output with updated version
            mv "$TEMP_PRETEXT" {output.pretext}
            
            echo "[GEP2] ✅ Added $TRACK_NAME"
        done
        
        echo ""
        echo "[GEP2] ✅ All tracks added to pretext"
        ls -lh {output.pretext}
        """