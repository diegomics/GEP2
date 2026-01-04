# ═══════════════════════════════════════════════════════════════════════════════
# GEP2 - Hi-C Analysis Rules
# ═══════════════════════════════════════════════════════════════════════════════
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


# ═══════════════════════════════════════════════════════════════════════════════
# INPUT FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

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


# ═══════════════════════════════════════════════════════════════════════════════
# RULES
# ═══════════════════════════════════════════════════════════════════════════════

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
        r2_str = lambda w: ",".join(get_hic_reads_r2(w))
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
        
        mkdir -p {params.outdir}
        
        # Create temp directory
        TEMP_DIR="$(mktemp -d "$GEP2_TMP/GEP2_chromap_{wildcards.species}_{wildcards.asm_basename}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        
        # Map with chromap, output pairs format
        chromap --preset hic \
            -x {input.index} \
            -r {input.asm} \
            -1 {params.r1_str} \
            -2 {params.r2_str} \
            -t {threads} \
            --pairs \
            -o "$TEMP_DIR/{wildcards.asm_basename}.pairs"
        
        echo "[GEP2] Compressing pairs file..."
        pigz -p {threads} -c "$TEMP_DIR/{wildcards.asm_basename}.pairs" > {output.pairs}
        
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
        echo ""
        echo "=== Pairs Statistics ==="
        head -50 {output.stats}
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
        
        # Index if needed and create genome file
        samtools faidx {input.asm}
        awk '{{print $1"\t"$2}}' {input.asm}.fai > {output.genome}
        
        echo "[GEP2] ✅ Genome file created"
        echo ""
        echo "=== Scaffold sizes (first 20) ==="
        head -20 {output.genome}
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


rule E06_pairtools_split:
    """Convert .pairs to sorted BAM using pairtools."""
    input:
        pairs = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}.pairs.gz"
        ),
        asm = get_hic_asm_input
    output:
        bam = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}_hic.bam"
        )
    threads: cpu_func("pretext_map")
    resources:
        mem_mb = mem_func("pretext_map"),
        runtime = time_func("pretext_map")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E06_pairtools_split_{asm_basename}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Converting pairs to BAM for {wildcards.asm_basename}"
        
        # Create temp directory
        TEMP_DIR="$(mktemp -d "$GEP2_TMP/GEP2_pairtools_{wildcards.species}_{wildcards.asm_basename}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        
        # Split pairs to SAM
        pairtools split \
            --output-sam "$TEMP_DIR/{wildcards.asm_basename}.sam" \
            {input.pairs}
        
        echo "[GEP2] Converting SAM to sorted BAM..."
        samtools sort \
            -@ {threads} \
            -o {output.bam} \
            -T "$TEMP_DIR/sort_tmp" \
            "$TEMP_DIR/{wildcards.asm_basename}.sam"
        
        echo "[GEP2] Indexing BAM..."
        samtools index -@ {threads} {output.bam}
        
        echo "[GEP2] ✅ BAM file created"
        """


rule E07_pretext_map:
    """Create PretextMap from Hi-C alignments."""
    input:
        bam = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}_hic.bam"
        )
    output:
        pretext = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "hic", "{asm_basename}", "{asm_basename}.pretext"
        )
    params:
        mapq = config.get("HIC_MAPQ", 10)
    threads: cpu_func("pretext_map")
    resources:
        mem_mb = mem_func("pretext_map"),
        runtime = time_func("pretext_map")
    container: CONTAINERS["hic_analysis"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E07_pretext_map_{asm_basename}.log"
        )
    benchmark:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "E07_pretext_map_{asm_basename}_benchmark.txt"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Creating PretextMap for {wildcards.asm_basename}"
        echo "[GEP2] BAM file: {input.bam}"
        echo "[GEP2] MAPQ threshold: {params.mapq}"
        
        samtools view -h {input.bam} | PretextMap \
            -o {output.pretext} \
            --sortby length \
            --sortorder descend \
            --mapq {params.mapq}
        
        echo "[GEP2] ✅ PretextMap created successfully"
        
        # Show file size as a basic sanity check
        echo ""
        echo "=== Output File ==="
        ls -lh {output.pretext}
        """
