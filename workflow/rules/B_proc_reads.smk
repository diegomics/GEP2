# ═══════════════════════════════════════════════════════════════════════════════
# GEP2 - Read Processing Rules
# ═══════════════════════════════════════════════════════════════════════════════

# Note: The following are defined in the main Snakefile:
#   - _MISSING_IN: Placeholder for missing input files
#   - _mk_dirs(): Create log directories
#   - _is_10x(): Check if read type is 10x
#   - _get_centralized_reads_dict(): Get centralized read files
#   - _enumerate_centralized_groups(): Generate read groups
#   - _pick_centralized_group(): Pick read group by wildcards
#   - normalize_read_type(): Normalize read type names
#   - DOWNLOAD_MANIFEST, DOWNLOAD_MANIFEST_DICT: Download manifest data


# ═══════════════════════════════════════════════════════════════════════════════
# RULE ORDER
# ═══════════════════════════════════════════════════════════════════════════════

ruleorder: B04_check_uli_primers > B00_centralize_reads
ruleorder: B05_filter_hifi_adapters > B00_centralize_reads
ruleorder: B03_ont_correction > B00_centralize_reads
ruleorder: B02_trim_pe_fastp > B00_centralize_reads
ruleorder: B00_centralize_reads > B00_link_long_reads > B01_compress_long_reads
ruleorder: B00_centralize_reads > B00_link_pe_reads > B01_compress_pe_reads


# ═══════════════════════════════════════════════════════════════════════════════
# INPUT FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def _get_source_for_centralized_read(wildcards):
    """Find the source file that needs to be symlinked to the centralized location."""
    filename = wildcards.filename
    
    # Extract identifier from filename (e.g., "illumina_Path1_ERR5709830_1.fq.gz" -> "ERR5709830_1")
    identifier = re.sub(r'^(hifi|ont|illumina|10x|hic)_Path\d+_', '', filename, flags=re.IGNORECASE)
    identifier = identifier.replace('.fq.gz', '').replace('.fastq.gz', '')
    
    base_path = os.path.join(
        config["OUT_FOLDER"], "GEP2_results", "downloaded_data", 
        wildcards.species, "reads"  # NO assembly folder
    )
    
    # Search in downloaded_data (species-level, no assembly folder)
    for rt_variant in [wildcards.read_type, wildcards.read_type.upper(), wildcards.read_type.lower()]:
        for ext in [".fastq.gz", ".fq.gz"]:
            pattern = os.path.join(base_path, rt_variant, f"{identifier}{ext}")
            if os.path.exists(pattern):
                return pattern
            # Also try glob for flexibility
            matches = glob.glob(pattern)
            if matches:
                return matches[0]
    
    # Check original paths in samples_config
    try:
        sp_data = samples_config["sp_name"][wildcards.species]
        for asm_id, asm_data in sp_data["asm_id"].items():
            read_type_dict = asm_data.get("read_type", {})
            
            for read_type, rt_data in read_type_dict.items():
                if normalize_read_type(read_type) != wildcards.read_type.lower():
                    continue
                
                read_files = rt_data.get("read_files", {})
                for path_key, path_value in read_files.items():
                    if not path_value or path_value == "None":
                        continue
                    
                    if isinstance(path_value, str) and "," in path_value:
                        paths = [p.strip() for p in path_value.split(",")]
                    elif isinstance(path_value, list):
                        paths = path_value
                    else:
                        paths = [str(path_value)]
                    
                    for path in paths:
                        path_base = os.path.basename(path)
                        path_id = path_base.replace('.fq.gz', '').replace('.fastq.gz', '').replace('.fq', '').replace('.fastq', '')
                        
                        if identifier in path_id or path_id in identifier:
                            if os.path.exists(path):
                                return path
    except (KeyError, TypeError, AttributeError):
        pass
    
    # Fallback: construct expected download path (NO assembly folder)
    expected_path = os.path.join(
        base_path, wildcards.read_type.lower(), f"{identifier}.fastq.gz"
    )
    return expected_path


def _linkable_long_src(w):
    """Get source file for linking long reads."""
    g = _pick_centralized_group(w)
    yaml_path = str(g["r"])
    
    if not yaml_path.endswith(".gz"):
        return _MISSING_IN
    
    # Convert centralized path to downloaded path
    # From: .../data/{species}/reads/{read_type}/hifi_Path1_ERR12205285.fq.gz
    # To:   .../downloaded_data/{species}/reads/{read_type}/ERR12205285.fastq.gz
    
    if "/data/" in yaml_path and "_Path" in os.path.basename(yaml_path):
        # Extract accession from filename (e.g., hifi_Path1_ERR12205285.fq.gz -> ERR12205285)
        basename = os.path.basename(yaml_path)
        # Remove read_type prefix and Path index, and change extension
        parts = basename.replace(".fq.gz", "").split("_")
        # Find the accession part (starts with ERR, SRR, DRR, or is the last non-empty part)
        accession = None
        for part in parts:
            if part.startswith(("ERR", "SRR", "DRR")) or (not accession and part and part not in ["hifi", "ont", "hic", "illumina", "10x"]):
                accession = part
        
        if accession:
            downloaded_path = yaml_path.replace("/data/", "/downloaded_data/")
            downloaded_path = os.path.join(
                os.path.dirname(downloaded_path),
                f"{accession}.fastq.gz"
            )
            return downloaded_path
    
    return yaml_path


def _compressible_long_src(w):
    """Get source file for compressing long reads."""
    g = _pick_centralized_group(w)
    return g["r"] if not str(g["r"]).endswith(".gz") else _MISSING_IN


def _linkable_pe_r1(w):
    """Get R1 source for linking PE reads."""
    if config.get("TRIM_PE", True):
        return _MISSING_IN
    g = _pick_centralized_group(w)
    yaml_path = str(g["r1"])
    
    if not yaml_path.endswith(".gz"):
        return _MISSING_IN
    
    # Convert centralized path to downloaded path
    return _convert_centralized_to_downloaded(yaml_path, w.read_type)


def _linkable_pe_r2(w):
    """Get R2 source for linking PE reads."""
    if config.get("TRIM_PE", True):
        return _MISSING_IN
    g = _pick_centralized_group(w)
    yaml_path = str(g["r2"])
    
    if not yaml_path.endswith(".gz"):
        return _MISSING_IN
    
    # Convert centralized path to downloaded path
    return _convert_centralized_to_downloaded(yaml_path, w.read_type)


def _convert_centralized_to_downloaded(yaml_path, read_type):
    """Convert a centralized data path to the corresponding downloaded_data path."""
    if "/data/" not in yaml_path or "_Path" not in os.path.basename(yaml_path):
        return yaml_path
    
    basename = os.path.basename(yaml_path)
    # Remove read_type prefix and Path index: illumina_Path1_SRR123_1.fq.gz -> SRR123_1
    parts = basename.replace(".fq.gz", "").replace(".fastq.gz", "").split("_")
    
    # Reconstruct accession with pair suffix
    # e.g., ["illumina", "Path1", "SRR123", "1"] -> "SRR123_1"
    accession_parts = []
    skip_next = False
    for i, part in enumerate(parts):
        if skip_next:
            skip_next = False
            continue
        if part.lower() in ["hifi", "ont", "hic", "illumina", "10x"]:
            continue
        if part.startswith("Path"):
            continue
        accession_parts.append(part)
    
    accession = "_".join(accession_parts)
    
    downloaded_path = yaml_path.replace("/data/", "/downloaded_data/")
    downloaded_path = os.path.join(
        os.path.dirname(downloaded_path),
        f"{accession}.fastq.gz"
    )
    return downloaded_path


def _compressible_pe_r1(w):
    """Get R1 source for compressing PE reads."""
    if config.get("TRIM_PE", True):
        return _MISSING_IN
    g = _pick_centralized_group(w)
    needs = (not str(g["r1"]).endswith(".gz")) or (not str(g["r2"]).endswith(".gz"))
    return g["r1"] if needs else _MISSING_IN


def _compressible_pe_r2(w):
    """Get R2 source for compressing PE reads."""
    if config.get("TRIM_PE", True):
        return _MISSING_IN
    g = _pick_centralized_group(w)
    needs = (not str(g["r1"]).endswith(".gz")) or (not str(g["r2"]).endswith(".gz"))
    return g["r2"] if needs else _MISSING_IN


def _get_ont_reads_for_correction(w):
    """Get ONT reads for correction if enabled."""
    if not config.get("CORRECT_ONT", True):
        return _MISSING_IN
    
    if w.read_type.lower() != "ont":
        return _MISSING_IN
    
    g = _pick_centralized_group(w)
    return g["r"]


def _get_uli_input(w):
    """Get HiFi reads for ULI primer check."""
    g = _pick_centralized_group(w)
    return g["r"]


def _get_hifiasm_input(w):
    """Get HiFi reads for adapter filtering."""
    g = _pick_centralized_group(w)
    return g["r"]


def _get_fastqc_pe_input(w):
    """Get input for FastQC (PE reads)."""
    return os.path.join(
        config["OUT_FOLDER"], "GEP2_results", "data", w.species,
        "reads", w.read_type, f"{w.read_type}_Path{w.idx}_{w.base}_{w.read}.fq.gz"
    )


def _get_fastqc_long_input(w):
    """Get input for FastQC (long reads)."""
    return os.path.join(
        config["OUT_FOLDER"], "GEP2_results", "data", w.species,
        "reads", w.read_type, f"{w.read_type}_Path{w.idx}_{w.base}.fq.gz"
    )


def _get_nanoplot_input(w):
    """Get input for NanoPlot."""
    base_dir = os.path.join(
        config["OUT_FOLDER"], "GEP2_results", "data", w.species,
        "reads", w.read_type
    )
    
    if "_filtered" in w.filename or "_corrected" in w.filename:
        return os.path.join(base_dir, "processed", f"{w.filename}.fq.gz")
    else:
        return os.path.join(base_dir, f"{w.filename}.fq.gz")


def _get_qc_reports_for_multiqc(w):
    """Collect all QC reports for a species/read_type combination."""
    species = w.species
    read_type = w.read_type
    read_type_lower = read_type.lower()
    
    reports = []
    report_dir = os.path.join(
        config["OUT_FOLDER"], "GEP2_results", "data", species,
        "reads", read_type_lower, "processed", "reports"
    )

    for grp in _enumerate_centralized_groups(species, read_type):
        idx = grp["idx"]
        base = grp["base"]
        
        if grp["kind"] == "long":
            rt_lower = read_type_lower
            
            if read_type_lower == "hifi":
                reports.extend([
                    os.path.join(report_dir, f"{rt_lower}_Path{idx}_{base}_uli_mqc.yaml"),
                    os.path.join(report_dir, f"{rt_lower}_Path{idx}_{base}_fastqc.zip")
                ])
                
                if config.get("TRIM_HIFI", True):
                    reports.extend([
                        os.path.join(report_dir, f"{rt_lower}_Path{idx}_{base}_haf_mqc.yaml"),
                        os.path.join(report_dir, f"{rt_lower}_Path{idx}_{base}_filtered_nanoplot")
                    ])
                else:
                    reports.append(
                        os.path.join(report_dir, f"{rt_lower}_Path{idx}_{base}_nanoplot")
                    )
                    
            elif read_type_lower == "ont":
                reports.append(
                    os.path.join(report_dir, f"{rt_lower}_Path{idx}_{base}_fastqc.zip")
                )
                
                if config.get("CORRECT_ONT", True):
                    reports.append(
                        os.path.join(report_dir, f"{rt_lower}_Path{idx}_{base}_corrected_nanoplot")
                    )
                else:
                    reports.append(
                        os.path.join(report_dir, f"{rt_lower}_Path{idx}_{base}_nanoplot")
                    )

        else:
            reports.extend([
                os.path.join(report_dir, f"{read_type_lower}_Path{idx}_{base}_1_fastqc.zip"),
                os.path.join(report_dir, f"{read_type_lower}_Path{idx}_{base}_2_fastqc.zip")
            ])                    

            if config.get("TRIM_PE", True):
                reports.append(os.path.join(report_dir, f"{read_type_lower}_Path{idx}_{base}_fastp.json"))
    
    return reports


# ═══════════════════════════════════════════════════════════════════════════════
# RULES - Read Centralization
# ═══════════════════════════════════════════════════════════════════════════════

rule B00_centralize_reads:
    """Create symlinks from downloaded_data to centralized data/reads location."""
    input:
        source = _get_source_for_centralized_read
    output:
        centralized = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "{filename}"
        )
    wildcard_constraints:
        read_type = r"[^/]+",
        filename = r"[^/]+\.fq\.gz"
    threads: cpu_func("light_task")
    resources:
        mem_mb = mem_func("light_task"),
        runtime = time_func("light_task")
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs", "B00_centralize_{filename}.log"
        )
    run:
        _mk_dirs(log[0])
        
        src = input.source
        dest = output.centralized
        
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        
        if os.path.exists(src):
            real_src = os.path.realpath(src)
            cmd = f'ln -sf "{real_src}" "{dest}" && echo "Centralized {src} to {dest}" > {log} 2>&1'
        else:
            cmd = f'echo "Warning: Source file not found: {src}, creating placeholder symlink" > {log} 2>&1 && ln -sf "{src}" "{dest}"'
        
        shell(cmd)


# ═══════════════════════════════════════════════════════════════════════════════
# RULES - Long Reads: Link/Compress
# ═══════════════════════════════════════════════════════════════════════════════

rule B00_link_long_reads:
    """Create symlink for already-compressed long reads."""
    input:
        src = _linkable_long_src
    output:
        linked = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "{read_type}_Path{idx}_{base}.fq.gz"
        )
    wildcard_constraints:
        read_type = r"(hifi|ont|hic)",
        idx = r"\d+",
        base = r"[^/]+(?<!_corrected)(?<!_filtered)"
    threads: cpu_func("light_task")
    resources:
        mem_mb = mem_func("light_task"),
        runtime = time_func("light_task")
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs",
            "B00_link_long_reads.{read_type}.Path{idx}.{base}.log"
        )
    run:
        grp = _pick_centralized_group(wildcards)
        
        if grp["kind"] != "long":
            raise ValueError("This target belongs to a paired-end job, not long reads.")
        if input.src == _MISSING_IN:
            raise ValueError("link_long_reads should only run when the source is .gz")
        
        _mk_dirs(log[0])
        
        os.makedirs(os.path.dirname(output.linked), exist_ok=True)
        
        cmd = f'ln -sf "$(realpath "{input.src}")" "{output.linked}" && echo "Linked {input.src} to {output.linked}" > {log} 2>&1'
        shell(cmd)


rule B01_compress_long_reads:
    """Compress uncompressed long reads."""
    input:
        src = _compressible_long_src
    output:
        compressed = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "{read_type}_Path{idx}_{base}.fq.gz"
        )
    wildcard_constraints:
        read_type = r"(hifi|ont|hic)",
        idx = r"\d+",
        base = r"[^/]+(?<!_corrected)(?<!_filtered)"
    threads: cpu_func("compress_reads")
    resources:
        mem_mb = mem_func("compress_reads"),
        runtime = time_func("compress_reads")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs",
            "B01_compress_long_reads.{read_type}.Path{idx}.{base}.log"
        )
    run:
        grp = _pick_centralized_group(wildcards)
        
        if grp["kind"] != "long":
            raise ValueError("This target belongs to a paired-end job, not long reads.")
        if input.src == _MISSING_IN:
            raise ValueError("compress_long_reads should only run when the source is not .gz")
        
        _mk_dirs(log[0])
        os.makedirs(os.path.dirname(output.compressed), exist_ok=True)
        
        cmd = f'pigz -p {threads} -c "{input.src}" > "{output.compressed}" 2> {log}'
        shell(cmd)


# ═══════════════════════════════════════════════════════════════════════════════
# RULES - Paired-End Reads: Link/Compress
# ═══════════════════════════════════════════════════════════════════════════════

rule B00_link_pe_reads:
    """Create symlinks for already-compressed PE reads (when TRIM_PE is off)."""
    input:
        r1 = _linkable_pe_r1,
        r2 = _linkable_pe_r2
    output:
        r1 = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "{read_type}_Path{idx}_{base}_1.fq.gz"
        ),
        r2 = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "{read_type}_Path{idx}_{base}_2.fq.gz"
        )
    wildcard_constraints:
        read_type = r"(illumina|10x)",
        idx = r"\d+",
        base = r"[^/]+(?<!_trimmed)"
    threads: cpu_func("light_task")
    resources:
        mem_mb = mem_func("light_task"),
        runtime = time_func("light_task")
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs",
            "B00_link_pe_reads.{read_type}.Path{idx}.{base}.log"
        )
    run:
        grp = _pick_centralized_group(wildcards)
        
        if grp["kind"] != "pe":
            raise ValueError("This target belongs to a long-read job, not paired-end.")
        if input.r1 == _MISSING_IN or input.r2 == _MISSING_IN:
            raise ValueError("link_pe_reads should only run when both sources are .gz")
        
        _mk_dirs(log[0])
        os.makedirs(os.path.dirname(output.r1), exist_ok=True)
        
        cmd = f'''
        ln -sf "$(realpath "{input.r1}")" "{output.r1}"
        ln -sf "$(realpath "{input.r2}")" "{output.r2}"
        echo "Linked PE reads" > {log} 2>&1
        '''
        shell(cmd)


rule B01_compress_pe_reads:
    """Compress uncompressed PE reads (when TRIM_PE is off)."""
    input:
        r1 = _compressible_pe_r1,
        r2 = _compressible_pe_r2
    output:
        r1 = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "{read_type}_Path{idx}_{base}_1.fq.gz"
        ),
        r2 = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "{read_type}_Path{idx}_{base}_2.fq.gz"
        )
    wildcard_constraints:
        read_type = r"(illumina|10x)",
        idx = r"\d+",
        base = r"[^/]+(?<!_trimmed)"
    threads: cpu_func("compress_reads")
    resources:
        mem_mb = mem_func("compress_reads"),
        runtime = time_func("compress_reads")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs",
            "B01_compress_pe_reads.{read_type}.Path{idx}.{base}.log"
        )
    run:
        grp = _pick_centralized_group(wildcards)
        
        if grp["kind"] != "pe":
            raise ValueError("This target belongs to a long-read job, not paired-end.")
        if input.r1 == _MISSING_IN or input.r2 == _MISSING_IN:
            raise ValueError("compress_pe_reads should only run when sources are not .gz")
        
        _mk_dirs(log[0])
        os.makedirs(os.path.dirname(output.r1), exist_ok=True)
        
        cmd = f'''
        pigz -p {threads} -c "{input.r1}" > "{output.r1}" 2> {log}
        pigz -p {threads} -c "{input.r2}" > "{output.r2}" 2>> {log}
        '''
        shell(cmd)


# ═══════════════════════════════════════════════════════════════════════════════
# RULES - PE Trimming
# ═══════════════════════════════════════════════════════════════════════════════

rule B02_trim_pe_fastp:
    """Trim and QC paired-end reads with fastp."""
    input:
        r1 = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "{read_type}_Path{idx}_{base}_1.fq.gz"
        ),
        r2 = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "{read_type}_Path{idx}_{base}_2.fq.gz"
        )
    output:
        r1 = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "{read_type}_Path{idx}_{base}_1_trimmed.fq.gz"
        ),
        r2 = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "{read_type}_Path{idx}_{base}_2_trimmed.fq.gz"
        ),
        json = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "reports", "{read_type}_Path{idx}_{base}_fastp.json"
        ),
        html = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "reports", "{read_type}_Path{idx}_{base}_fastp.html"
        )
    wildcard_constraints:
        read_type = r"(illumina|10x)",
        idx = r"\d+",
        base = r"[^/]+"
    threads: cpu_func("trim_pe")
    resources:
        mem_mb = mem_func("trim_pe"),
        runtime = time_func("trim_pe")
    container: CONTAINERS["gep2_base"]
    params:
        extra_opts = lambda w: "--trim_front1 23" if _is_10x(w.read_type) else ""
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs",
            "B02_fastp.{read_type}_Path{idx}_{base}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        mkdir -p $(dirname {output.r1})
        mkdir -p $(dirname {output.json})
        
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              --thread {threads} \
              --json {output.json} \
              --html {output.html} \
              --detect_adapter_for_pe \
              --length_required 50 \
              --overrepresentation_analysis \
              {params.extra_opts}
        
        echo "[GEP2] fastp completed for {wildcards.species}/{wildcards.read_type}"
        """


# ═══════════════════════════════════════════════════════════════════════════════
# RULES - Long Read Processing
# ═══════════════════════════════════════════════════════════════════════════════

rule B03_ont_correction:
    """Correct ONT reads using Hifiasm (if enabled)."""
    input:
        reads = _get_ont_reads_for_correction
    output:
        corrected = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "{read_type}_Path{idx}_{base}_corrected.fq.gz"
        )
    wildcard_constraints:
        read_type = r"(ont)",
        idx = r"\d+",
        base = r"[^/]+"
    threads: cpu_func("correct_ont")
    resources:
        mem_mb = mem_func("correct_ont"),
        runtime = time_func("correct_ont")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs",
            "B03_hifiasm_correct.{read_type}_Path{idx}_{base}.log"
        )
    shell:
        r'''
        set -euo pipefail
        
        mkdir -p $(dirname {output.corrected}) $(dirname {log})
        
        BASE="${{GEP2_TMP:-$TMPDIR}}"
        mkdir -p "$BASE"
        TEMP_DIR=$(mktemp -d "$BASE/GEP2_ont_correct_{wildcards.species}_{wildcards.read_type}_{wildcards.base}_XXXXXX")
        cd $TEMP_DIR
        
        echo "Starting ONT correction with Hifiasm..." > {log}
        
        hifiasm -o tmp --ont -t {threads} -l0 --write-ec {input.reads} >> hifi.log 2>&1 &
        HIFI_PID=$!
        
        LAST_SIZE=0
        STAGNANT_COUNT=0
        LOOP_COUNT=0
        
        while kill -0 $HIFI_PID 2>/dev/null; do
            LOOP_COUNT=$((LOOP_COUNT + 1))
            
            if grep -q "Reads has been written" hifi.log 2>/dev/null; then
                sleep 30
                kill -TERM $HIFI_PID 2>/dev/null || true
                break
            fi
            
            if [[ -f "tmp.ec.fq" ]]; then
                CURRENT_SIZE=$(stat -c%s "tmp.ec.fq" 2>/dev/null || echo 0)
                
                if [[ $CURRENT_SIZE -eq $LAST_SIZE ]] && [[ $CURRENT_SIZE -gt 1000000 ]]; then
                    STAGNANT_COUNT=$((STAGNANT_COUNT + 1))
                    
                    if [[ $STAGNANT_COUNT -ge 4 ]]; then
                        kill -TERM $HIFI_PID 2>/dev/null || true
                        break
                    fi
                else
                    STAGNANT_COUNT=0
                    LAST_SIZE=$CURRENT_SIZE
                fi
            fi
            
            if [[ $LOOP_COUNT -gt 2880 ]]; then
                kill -TERM $HIFI_PID 2>/dev/null || true
                break
            fi
            
            sleep 30
        done
        
        wait $HIFI_PID || EXIT_CODE=$?
        
        if [[ -f "tmp.ec.fq" ]] && [[ $(stat -c%s "tmp.ec.fq") -gt 1000000 ]]; then
            pigz -p {threads} -c tmp.ec.fq > {output.corrected}
        else
            echo "ERROR: No valid EC file produced" >> {log}
            exit 1
        fi
        
        cd /
        rm -rf $TEMP_DIR
        '''


rule B04_check_uli_primers:
    """Check for ULI primers in HiFi reads (100k subsample)."""
    input:
        reads = _get_uli_input
    output:
        report = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "reports",
            "{read_type}_Path{idx}_{base}_uli_mqc.yaml"
        )
    wildcard_constraints:
        read_type = r"(hifi)",
        idx = r"\d+",
        base = r"[^/]+(?<!_filtered)(?<!_corrected)"
    threads: cpu_func("reads_qc")
    resources:
        mem_mb = mem_func("reads_qc"),
        runtime = time_func("reads_qc")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs",
            "B04_uli.{read_type}_Path{idx}_{base}.log"
        )
    shell:
        r'''
        set -euo pipefail
        exec > {log} 2>&1
        
        mkdir -p $(dirname {output.report})
        
        BASE="${{GEP2_TMP:-$TMPDIR}}"
        mkdir -p "$BASE"
        TEMP_DIR="$(mktemp -d "$BASE/GEP2_uli_check_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        
        seqtk sample -s789 {input.reads} 100000 > subset.fq
        
        hits=$(awk 'NR%4==2' subset.fq | grep -E -c 'AAGCAGTGGTATCAACGCAGAGTACT|AGTACTCTGCGTTGATACCACTGCTT') || hits=0
        reads=$(awk 'NR%4==2' subset.fq | wc -l) || reads=0
        ulipct=$(awk -v h="$hits" -v r="$reads" 'BEGIN{{printf("%.2f", r>0 ? 100*h/r : 0)}}')
        
        cat > {output.report} << EOF
id: "uli_check"
section_name: "ULI Primer Check"
description: "ULI adapter presence in HiFi reads (100k subsample)"
plot_type: "generalstats"
headers:
  - uli_primer_pct:
      title: "ULI %"
      max: 100
      min: 0
      suffix: "%"
      scale: "RdYlGn-rev"
data:
  {wildcards.read_type}_Path{wildcards.idx}_{wildcards.base}:
    uli_primer_pct: $ulipct
EOF
        '''


rule B05_filter_hifi_adapters:
    """Filter HiFi reads with HiFiAdapterFilt."""
    input:
        reads = _get_hifiasm_input
    output:
        filtered = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "{read_type}_Path{idx}_{base}_filtered.fq.gz"
        ),
        report = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "reports",
            "{read_type}_Path{idx}_{base}_haf_mqc.yaml"
        ),
        stats = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "reports",
            "{read_type}_Path{idx}_{base}_haf.stats"
        )
    wildcard_constraints:
        read_type = r"(hifi)",
        idx = r"\d+",
        base = r"[^/]+(?<!_filtered)(?<!_corrected)"
    threads: cpu_func("filter_hifi")
    resources:
        mem_mb = mem_func("filter_hifi"),
        runtime = time_func("filter_hifi")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs",
            "B05_hifiadapterfilt.{read_type}_Path{idx}_{base}.log"
        )
    shell:
        r'''
        set -euo pipefail
        exec > {log} 2>&1
        
        mkdir -p $(dirname {output.filtered})
        mkdir -p $(dirname {output.report})
        
        BASE="${{GEP2_TMP:-$TMPDIR}}"
        mkdir -p "$BASE"
        TEMP_DIR="$(mktemp -d "$BASE/GEP2_hifi_filter_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        
        tot=$(zcat {input.reads} | awk 'NR%4==1' | wc -l)
        
        ln -s {input.reads} input.fq.gz
        
        /opt/HiFiAdapterFilt/hifiadapterfilt.sh -p input -t {threads} -o .
        
        mv input.filt.fastq.gz {output.filtered}
        
        if [[ -f "input.stats" ]]; then
            cp input.stats {output.stats}
        else
            echo "No stats file produced" > {output.stats}
        fi
        
        kept=$(zcat {output.filtered} | awk 'NR%4==1' | wc -l)
        rem=$((tot - kept))
        pct=$(awk -v r="$rem" -v t="$tot" 'BEGIN{{printf("%.2f", t>0 ? 100*r/t : 0)}}')
        
        cat > {output.report} << EOF
id: "hifiadapterfilt"
section_name: "HiFiAdapterFilt"
description: "HiFi adapter filtering statistics"
plot_type: "generalstats"
headers:
  - haf_removed_pct:
      title: "Removed %"
      max: 100
      min: 0
      suffix: "%"
      scale: "RdYlGn-rev"
data:
  {wildcards.read_type}_Path{wildcards.idx}_{wildcards.base}:
    haf_removed_pct: $pct
EOF
        '''


# ═══════════════════════════════════════════════════════════════════════════════
# RULES - QC
# ═══════════════════════════════════════════════════════════════════════════════

rule B06_fastqc_pe_reads:
    """Run FastQC on paired-end reads (250k subsample)."""
    input:
        reads = _get_fastqc_pe_input
    output:
        html = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "reports",
            "{read_type}_Path{idx}_{base}_{read}_fastqc.html"
        ),
        zip = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "reports",
            "{read_type}_Path{idx}_{base}_{read}_fastqc.zip"
        )
    wildcard_constraints:
        read_type = r"(illumina|10x)",
        read = r"(1|2)",
        idx = r"\d+",
        base = r"[^/]+"
    threads: cpu_func("reads_qc")
    resources:
        mem_mb = mem_func("reads_qc"),
        runtime = time_func("reads_qc")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs",
            "B06_fastqc.{read_type}_Path{idx}_{base}_{read}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        mkdir -p $(dirname {output.html})
        
        BASE="${{GEP2_TMP:-$TMPDIR}}"
        mkdir -p "$BASE"
        TEMP_DIR="$(mktemp -d "$BASE/GEP2_fastqc_XXXXXX")"
        cd "$TEMP_DIR"
        
        TEMP_SUB="{wildcards.read_type}_Path{wildcards.idx}_{wildcards.base}_{wildcards.read}.fq"
        seqtk sample -s123 {input.reads} 250000 > $TEMP_SUB
        
        fastqc -t {threads} -o . $TEMP_SUB
        
        mv "{wildcards.read_type}_Path{wildcards.idx}_{wildcards.base}_{wildcards.read}_fastqc.html" "{output.html}"
        mv "{wildcards.read_type}_Path{wildcards.idx}_{wildcards.base}_{wildcards.read}_fastqc.zip" "{output.zip}"
        
        cd /
        rm -rf "$TEMP_DIR"
        """


rule B06_fastqc_long_reads:
    """Run FastQC on long reads (250k subsample)."""
    input:
        reads = _get_fastqc_long_input
    output:
        html = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "reports",
            "{read_type}_Path{idx}_{base}_fastqc.html"
        ),
        zip = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "reports",
            "{read_type}_Path{idx}_{base}_fastqc.zip"
        )
    wildcard_constraints:
        read_type = r"(hifi|ont)",
        idx = r"\d+",
        base = r"[^/]+"
    threads: cpu_func("reads_qc")
    resources:
        mem_mb = mem_func("reads_qc"),
        runtime = time_func("reads_qc")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs",
            "B06_fastqc.{read_type}_Path{idx}_{base}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        mkdir -p $(dirname {output.html})
        
        BASE="${{GEP2_TMP:-$TMPDIR}}"
        mkdir -p "$BASE"
        TEMP_DIR="$(mktemp -d "$BASE/GEP2_fastqc_XXXXXX")"
        cd "$TEMP_DIR"
        
        TEMP_SUB="{wildcards.read_type}_Path{wildcards.idx}_{wildcards.base}.fq"
        seqtk sample -s123 {input.reads} 250000 > $TEMP_SUB
        
        fastqc -t {threads} -o . $TEMP_SUB
        
        mv "{wildcards.read_type}_Path{wildcards.idx}_{wildcards.base}_fastqc.html" "{output.html}"
        mv "{wildcards.read_type}_Path{wildcards.idx}_{wildcards.base}_fastqc.zip" "{output.zip}"
        
        cd /
        rm -rf "$TEMP_DIR"
        """


rule B07_nanoplot_long_reads:
    """Run NanoPlot QC on long reads."""
    input:
        reads = _get_nanoplot_input
    output:
        dir = directory(os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "reports",
            "{filename}_nanoplot"
        ))
    wildcard_constraints:
        read_type = r"(hifi|ont)",
        filename = r"[^/]+"
    threads: cpu_func("reads_qc")
    resources:
        mem_mb = mem_func("reads_qc"),
        runtime = time_func("reads_qc")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs",
            "B07_nanoplot.{filename}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        mkdir -p {output.dir}
        
        TEMP_DIR=$(mktemp -d "$TMPDIR/nanoplot_XXXXXX")
        cd $TEMP_DIR
        
        seqtk sample -s456 {input.reads} 250000 | gzip -c > subsample.fq.gz
        
        if NanoPlot -t {threads} --fastq subsample.fq.gz -o {output.dir} --tsv_stats; then
            echo "NanoPlot completed successfully"
        else
            cat > {output.dir}/NanoStats.txt << 'EOF'
General summary:
Number of reads:	0
Number of bases:	0
Median read length:	0.0
Mean read length:	0.0
Read length N50:	0
EOF
        fi
        
        cd /
        rm -rf $TEMP_DIR
        """


# ═══════════════════════════════════════════════════════════════════════════════
# RULES - MultiQC Aggregation
# ═══════════════════════════════════════════════════════════════════════════════

rule B08_aggregate_read_qc:
    """Aggregate all QC reports with MultiQC."""
    input:
        _get_qc_reports_for_multiqc
    output:
        html = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "reports",
            "multiqc_report.html"
        ),
        data = directory(os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "reports",
            "multiqc_report_data"
        ))
    threads: cpu_func("reads_qc")
    resources:
        mem_mb = mem_func("reads_qc"),
        runtime = time_func("reads_qc")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs",
            "B08_multiqc.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        mkdir -p $(dirname {output.html})
        
        TEMP_LIST=$(mktemp)
        for item in {input}; do
            echo "$item" >> $TEMP_LIST
        done
        
        multiqc -f -n multiqc_report \
            -o $(dirname {output.html}) \
            --file-list $TEMP_LIST
        
        rm $TEMP_LIST
        """
