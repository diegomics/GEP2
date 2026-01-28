# -------------------------------------------------------------------------------
# GEP2 - K-mer Analysis Rules
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# INPUT FUNCTIONS
# -------------------------------------------------------------------------------

def _get_processed_read_path(species, read_type, idx, base):
    """Get the processed read path for k-mer counting."""
    read_type_lower = read_type.lower()
    base_dir = os.path.join(
        config["OUT_FOLDER"], "GEP2_results", "data", species,
        "reads", read_type_lower
    )
    
    if read_type_lower in ["illumina", "10x"]:
        if config.get("TRIM_PE", True):
            return os.path.join(base_dir, "processed", f"{read_type_lower}_Path{idx}_{base}_1_trimmed.fq.gz")
        else:
            return os.path.join(base_dir, f"{read_type_lower}_Path{idx}_{base}_1.fq.gz")
    elif read_type_lower == "hifi":
        if config.get("TRIM_HIFI", True):
            return os.path.join(base_dir, "processed", f"hifi_Path{idx}_{base}_filtered.fq.gz")
        else:
            return os.path.join(base_dir, f"hifi_Path{idx}_{base}.fq.gz")
    elif read_type_lower == "ont":
        if config.get("CORRECT_ONT", True):
            return os.path.join(base_dir, "processed", f"ont_Path{idx}_{base}_corrected.fq.gz")
        else:
            return os.path.join(base_dir, f"ont_Path{idx}_{base}.fq.gz")
    else:
        return os.path.join(base_dir, f"{read_type_lower}_Path{idx}_{base}.fq.gz")


def get_per_read_kmer_input(wildcards):
    """Input function for per-read k-mer database construction."""
    species = wildcards.species
    read_type = wildcards.read_type.lower()
    base = wildcards.base
    
    # Find the idx for this base by looking up in centralized groups
    for grp in _enumerate_centralized_groups(species, read_type):
        if grp["base"] == base:
            idx = grp["idx"]
            return _get_processed_read_path(species, read_type, idx, base)
    
    # Fallback - construct path directly (try common patterns)
    base_dir = os.path.join(
        config["OUT_FOLDER"], "GEP2_results", "data", species,
        "reads", read_type
    )
    
    # Try to find any matching file
    if read_type == "hifi":
        if config.get("TRIM_HIFI", True):
            # Look for filtered files with any Path index
            import glob
            pattern = os.path.join(base_dir, "processed", f"hifi_Path*_{base}_filtered.fq.gz")
            matches = glob.glob(pattern)
            if matches:
                return matches[0]
        else:
            pattern = os.path.join(base_dir, f"hifi_Path*_{base}.fq.gz")
            matches = glob.glob(pattern)
            if matches:
                return matches[0]
    
    elif read_type == "ont":
        if config.get("CORRECT_ONT", True):
            pattern = os.path.join(base_dir, "processed", f"ont_Path*_{base}_corrected.fq.gz")
            matches = glob.glob(pattern)
            if matches:
                return matches[0]
        else:
            pattern = os.path.join(base_dir, f"ont_Path*_{base}.fq.gz")
            matches = glob.glob(pattern)
            if matches:
                return matches[0]
    
    elif read_type in ["illumina", "10x"]:
        if config.get("TRIM_PE", True):
            pattern = os.path.join(base_dir, "processed", f"{read_type}_Path*_{base}_1_trimmed.fq.gz")
            matches = glob.glob(pattern)
            if matches:
                return matches[0]
        else:
            pattern = os.path.join(base_dir, f"{read_type}_Path*_{base}_1.fq.gz")
            matches = glob.glob(pattern)
            if matches:
                return matches[0]
    
    raise ValueError(f"Could not find processed read file for {species}/{read_type}/{base}")


def get_assembly_kmer_db_inputs(wildcards):
    """Get the per-read k-mer DBs needed for an assembly."""
    # Check if k-mer analysis should be skipped for this assembly
    if _should_skip_analysis(wildcards.species, wildcards.asm_id, "kmer"):
        return []
    
    reads = _get_reads_for_assembly(wildcards.species, wildcards.asm_id)
    
    # Determine k-mer length (use priority read type for k-mer length decision)
    priority_rt = get_priority_read_type_for_assembly(wildcards.species, wildcards.asm_id)
    kmer_len = get_kmer_length(priority_rt) if priority_rt else 31
    
    inputs = []
    for r in reads:
        db_path = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", wildcards.species,
            "reads", r["read_type"], f"kmer_db_k{kmer_len}", f"{r['base']}.meryl"
        )
        inputs.append(db_path)
    
    return inputs


def get_merqury_db_input(wildcards):
    """Get the merged k-mer database path for this assembly."""
    # Check if k-mer analysis should be skipped for this assembly
    if _should_skip_analysis(wildcards.species, wildcards.asm_id, "kmer"):
        return []
    
    read_type = get_priority_read_type(wildcards.species)
    if not read_type:
        raise ValueError(f"No reads available for {wildcards.species}")
    kmer_len = get_kmer_length(read_type)
    return os.path.join(
        config["OUT_FOLDER"], "GEP2_results", wildcards.species,
        wildcards.asm_id, f"k{kmer_len}", f"{wildcards.asm_id}.meryl"
    )

def get_merqury_asm_inputs(wildcards):
    """Get assembly files for Merqury, in sorted order."""
    # Check if k-mer analysis should be skipped for this assembly
    if _should_skip_analysis(wildcards.species, wildcards.asm_id, "kmer"):
        return []
    
    asm_files = get_assembly_files(wildcards.species, wildcards.asm_id)
    return [v for k, v in sorted(asm_files.items()) if v and v != "None"]


def get_asm_count(wildcards):
    """Get number of assembly files for determining haploid/diploid mode."""
    asm_files = get_assembly_files(wildcards.species, wildcards.asm_id)
    return len([v for v in asm_files.values() if v and v != "None"])


# -------------------------------------------------------------------------------
# RULES - Per-Read K-mer Database Construction
# -------------------------------------------------------------------------------

rule C00_build_per_read_kmer_db:
    """Build Meryl k-mer database for a single read file."""
    input:
        reads = get_per_read_kmer_input
    output:
        meryl_db = directory(os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "kmer_db_k{kmer_len}", "{base}.meryl"
        ))
    wildcard_constraints:
        kmer_len = r"\d+",
        base = r"[^/]+"
    threads: cpu_func("kmer_count")
    resources:
        mem_mb = mem_func("kmer_count"),
        runtime = time_func("kmer_count")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs", "C00_build_kmer_db_k{kmer_len}_{base}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Building k-mer database for {wildcards.base}"
        echo "[GEP2] K-mer length: {wildcards.kmer_len}"
        echo "[GEP2] Input: {input.reads}"
        
        mkdir -p $(dirname {output.meryl_db})
        
        WORK_DIR="$(gep2_get_workdir 100)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GEP2_meryl_{wildcards.species}_{wildcards.base}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT

        cd $TEMP_DIR
        
        meryl k={wildcards.kmer_len} \\
              threads={threads} \\
              count \\
              {input.reads} \\
              output temp.meryl
        
        mv temp.meryl {output.meryl_db}
        
        echo "[GEP2] ✅ K-mer database complete: {output.meryl_db}"
        """


# -------------------------------------------------------------------------------
# RULES - Assembly-Specific K-mer Database (merge or symlink)
# -------------------------------------------------------------------------------

rule C00_merge_assembly_kmer_db:
    """Create assembly-specific k-mer database (symlink if 1 read, union-sum if multiple)."""
    input:
        dbs = get_assembly_kmer_db_inputs
    output:
        meryl_db = directory(os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "k{kmer_len}", "{asm_id}.meryl"
        )),
        hist = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "k{kmer_len}", "{asm_id}.hist"
        )
    wildcard_constraints:
        kmer_len = r"\d+"
    params:
        db_count = lambda w, input: len(input.dbs),
        db_list = lambda w, input: " ".join(input.dbs)
    threads: cpu_func("kmer_count")
    resources:
        mem_mb = mem_func("kmer_count"),
        runtime = time_func("kmer_count")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "C00_merge_kmer_db_k{kmer_len}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Creating assembly-specific k-mer database for {wildcards.asm_id}"
        echo "[GEP2] Input databases: {params.db_count}"
        
        WORK_DIR="$(gep2_get_workdir 100)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GEP2_merge_kmer_{wildcards.species}_{wildcards.asm_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        
        cd "$TEMP_DIR"
        
        if [ {params.db_count} -eq 1 ]; then
            echo "[GEP2] Single read - copying k-mer database"
            cp -r {input.dbs} merged.meryl
        else
            echo "[GEP2] Multiple reads - running union-sum"
            meryl k={wildcards.kmer_len} \
                threads={threads} \
                union-sum \
                output merged.meryl \
                {params.db_list}
        fi
        
        echo "[GEP2] Generating histogram"
        meryl histogram merged.meryl | sed 's/\\t/ /g' > merged.hist
        
        # Copy results to final location
        echo "[GEP2] Copying results to final location"
        mkdir -p $(dirname {output.meryl_db})
        rm -rf {output.meryl_db}
        cp -r merged.meryl {output.meryl_db}
        cp merged.hist {output.hist}
        
        echo "[GEP2] ✅ Assembly k-mer database complete"
        """



# -------------------------------------------------------------------------------
# RULES - GenomeScope2
# -------------------------------------------------------------------------------

rule C01_run_genomescope2:
    """Run GenomeScope2 analysis on assembly-specific k-mer histogram."""
    input:
        hist = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "k{kmer_len}", "{asm_id}.hist"
        )
    output:
        summary = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "k{kmer_len}", "genomescope2", "{asm_id}_summary.txt"
        ),
        model = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "k{kmer_len}", "genomescope2", "{asm_id}_model.txt"
        ),
        linear_plot = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "k{kmer_len}", "genomescope2", "{asm_id}_linear_plot.png"
        )
    wildcard_constraints:
        kmer_len = r"\d+"
    params:
        outdir = lambda w: os.path.join(
            config["OUT_FOLDER"], "GEP2_results", w.species, w.asm_id,
            f"k{w.kmer_len}", "genomescope2"
        ),
        ploidy = config.get("PLOIDY", 2)
    threads: cpu_func("genomescope")
    resources:
        mem_mb = mem_func("genomescope"),
        runtime = time_func("genomescope")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "C01_genomescope2_k{kmer_len}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Running GenomeScope2 for {wildcards.species}/{wildcards.asm_id}"
        echo "[GEP2] K-mer length: {wildcards.kmer_len}"
        echo "[GEP2] Ploidy: {params.ploidy}"
        
        mkdir -p {params.outdir}
        
        genomescope2 -i {input.hist} \\
                     -o {params.outdir} \\
                     -k {wildcards.kmer_len} \\
                     -p {params.ploidy} \\
                     -n {wildcards.asm_id}
        
        echo "[GEP2] ✅ GenomeScope2 complete"
        """


# -------------------------------------------------------------------------------
# RULES - Merqury
# -------------------------------------------------------------------------------

rule C02_run_merqury:
    """Run Merqury for assembly QV and completeness analysis."""
    input:
        meryl_db = get_merqury_db_input,
        assemblies = get_merqury_asm_inputs
    output:
        qv = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "merqury", "{asm_id}.qv"
        ),
        completeness = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "merqury", "{asm_id}.completeness.stats"
        )
    params:
        outdir = lambda w: os.path.join(
            config["OUT_FOLDER"], "GEP2_results", w.species, w.asm_id, "merqury"
        ),
        asm_count = get_asm_count,
        prefix = lambda w: w.asm_id
    threads: cpu_func("merqury")
    resources:
        mem_mb = mem_func("merqury"),
        runtime = time_func("merqury")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "{species}", "{asm_id}",
            "logs", "C02_merqury.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Running Merqury for {wildcards.species}/{wildcards.asm_id}"
        echo "[GEP2] K-mer database: {input.meryl_db}"
        echo "[GEP2] Assembly count: {params.asm_count}"
        
        export OMP_NUM_THREADS={threads}
        export MERQURY=/opt/conda/share/merqury 
        
        mkdir -p {params.outdir}
        
        WORK_DIR="$(gep2_get_workdir 50)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GEP2_merqury_{wildcards.species}_{wildcards.asm_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT

        cd $TEMP_DIR
        
        ln -sf {input.meryl_db} read_db.meryl
        
        link_assembly() {{
            local src="$1"
            local linkname="$2"
            local ext=""
            
            case "$src" in
                *.fasta.gz) ext=".fasta.gz" ;;
                *.fa.gz)    ext=".fasta.gz" ;;
                *.fna.gz)   ext=".fasta.gz" ;;
                *.fasta)    ext=".fasta" ;;
                *.fa)       ext=".fasta" ;;
                *.fna)      ext=".fasta" ;;
                *)          ext=".fasta" ;;
            esac
            
            ln -sf "$src" "${{linkname}}${{ext}}"
            echo "${{linkname}}${{ext}}"
        }}
        
        ASM_COUNT={params.asm_count}
        ASSEMBLIES="{input.assemblies}"
        
        if [ $ASM_COUNT -eq 1 ]; then
            echo "[GEP2] Running Merqury in HAPLOID mode"
            
            ASM1=$(echo "$ASSEMBLIES" | awk '{{print $1}}')
            ASM1_LINK=$(link_assembly "$ASM1" "asm1")
            
            merqury.sh read_db.meryl "$ASM1_LINK" {params.prefix}
            
        elif [ $ASM_COUNT -eq 2 ]; then
            echo "[GEP2] Running Merqury in DIPLOID mode"
            
            ASM1=$(echo "$ASSEMBLIES" | awk '{{print $1}}')
            ASM2=$(echo "$ASSEMBLIES" | awk '{{print $2}}')
            ASM1_LINK=$(link_assembly "$ASM1" "asm1")
            ASM2_LINK=$(link_assembly "$ASM2" "asm2")
            
            merqury.sh read_db.meryl "$ASM1_LINK" "$ASM2_LINK" {params.prefix}
            
        else
            echo "[GEP2] ERROR: Expected 1 or 2 assembly files, got $ASM_COUNT"
            exit 1
        fi
        
        # DEBUG: List all files created
        echo "[GEP2] Files created in temp directory:"
        ls -la
        echo ""
        echo "[GEP2] Looking for completeness files:"
        ls -la *completeness* 2>/dev/null || echo "No completeness files found"
        echo ""

        # Move results
        mv {params.prefix}.* {params.outdir}/ 2>/dev/null || true
        mv *.png {params.outdir}/ 2>/dev/null || true
        mv *.pdf {params.outdir}/ 2>/dev/null || true
        mv *.hist {params.outdir}/ 2>/dev/null || true
        mv *.wig {params.outdir}/ 2>/dev/null || true
        mv *.bed {params.outdir}/ 2>/dev/null || true
        mv asm*.meryl {params.outdir}/ 2>/dev/null || true
        mv completeness.stats {params.outdir}/{params.prefix}.completeness.stats 2>/dev/null || true
        
        echo "[GEP2] Files in output directory after move:"
        ls -la {params.outdir}/
        
        if [ ! -f {output.qv} ]; then
            echo "[GEP2] ERROR: QV file not created"
            exit 1
        fi
        
        echo "[GEP2] ✅ Merqury completed"
        echo "=== QV Summary ==="
        cat {output.qv}
        echo "=== Completeness Summary ==="
        cat {output.completeness}
        """

# -------------------------------------------------------------------------------
# RULES - Reads-Only Genome Profiling
# -------------------------------------------------------------------------------

def get_reads_only_kmer_db_inputs(wildcards):
    """Get per-read k-mer DBs for reads-only genome profiling."""
    kmer_len = wildcards.kmer_len
    read_type = wildcards.read_type
    species = wildcards.species
    
    inputs = []
    
    try:
        for asm_id, asm_data in samples_config["sp_name"][species]["asm_id"].items():
            if not _is_reads_only_entry(species, asm_id):
                continue
            
            read_type_dict = asm_data.get("read_type", {})
            
            for rt_key, rt_data in read_type_dict.items():
                if normalize_read_type(rt_key) != read_type:
                    continue
                
                read_files = rt_data.get("read_files", {})
                
                for path_key, path_value in sorted(read_files.items()):
                    if not path_value or path_value == "None":
                        continue
                    
                    # Extract base name
                    basename = os.path.basename(str(path_value))
                    base = re.sub(r'^(hifi|ont|illumina|10x|hic)_Path\d+_', '', basename, flags=re.IGNORECASE)
                    base = base.replace(".fq.gz", "").replace(".fastq.gz", "")
                    base = base.replace("_1", "").replace("_2", "")
                    base = base.replace("_filtered", "").replace("_corrected", "").replace("_trimmed", "")
                    
                    db_path = os.path.join(
                        config["OUT_FOLDER"], "GEP2_results", "data", species,
                        "reads", read_type, f"kmer_db_k{kmer_len}", f"{base}.meryl"
                    )
                    
                    if db_path not in inputs:
                        inputs.append(db_path)
                        
    except (KeyError, TypeError, AttributeError):
        pass
    
    return inputs


rule C10_merge_reads_only_kmer_db:
    """Merge k-mer databases for reads-only genome profiling."""
    input:
        dbs = get_reads_only_kmer_db_inputs
    output:
        meryl_db = directory(os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "k{kmer_len}", "{species}.meryl"
        )),
        hist = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "k{kmer_len}", "{species}.hist"
        )
    params:
        outdir = lambda w: os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", w.species,
            "reads", w.read_type, f"k{w.kmer_len}"
        ),
        db_list = lambda w, input: " ".join(input.dbs)
    threads: cpu_func("kmer_count")
    resources:
        mem_mb = mem_func("kmer_count"),
        runtime = time_func("kmer_count")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs",
            "C10_merge_reads_only_kmer_db_k{kmer_len}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Merging k-mer databases for reads-only profiling: {wildcards.species}"
        echo "[GEP2] Read type: {wildcards.read_type}"
        echo "[GEP2] K-mer length: {wildcards.kmer_len}"
        echo "[GEP2] Input databases: {params.db_list}"
        
        WORK_DIR="$(gep2_get_workdir 100)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GEP2_reads_only_kmer_{wildcards.species}_{wildcards.read_type}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        
        cd "$TEMP_DIR"
        
        DB_COUNT=$(echo "{params.db_list}" | wc -w)
        
        if [ "$DB_COUNT" -eq 1 ]; then
            echo "[GEP2] Single database - copying"
            cp -r {params.db_list} merged.meryl
        else
            echo "[GEP2] Multiple databases - merging with union-sum"
            meryl k={wildcards.kmer_len} \
                threads={threads} \
                union-sum \
                output merged.meryl \
                {params.db_list}
        fi
        
        echo "[GEP2] Generating histogram"
        meryl histogram merged.meryl | sed 's/\\t/ /g' > merged.hist
        
        # Copy results to final location
        echo "[GEP2] Copying results to final location"
        mkdir -p {params.outdir}
        rm -rf {output.meryl_db}
        cp -r merged.meryl {output.meryl_db}
        cp merged.hist {output.hist}
        
        echo "[GEP2] ✅ Reads-only k-mer database complete"
        """


rule C11_reads_only_genomescope2:
    """Run GenomeScope2 for reads-only genome profiling."""
    input:
        hist = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "k{kmer_len}", "{species}.hist"
        )
    output:
        summary = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "k{kmer_len}", "genomescope2",
            "{species}_summary.txt"
        ),
        model = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "k{kmer_len}", "genomescope2",
            "{species}_model.txt"
        ),
        linear_plot = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "k{kmer_len}", "genomescope2",
            "{species}_linear_plot.png"
        )
    params:
        outdir = lambda w: os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", w.species,
            "reads", w.read_type, f"k{w.kmer_len}", "genomescope2"
        ),
        ploidy = config.get("PLOIDY", 2)
    threads: cpu_func("genomescope")
    resources:
        mem_mb = mem_func("genomescope"),
        runtime = time_func("genomescope")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs",
            "C11_reads_only_genomescope2_k{kmer_len}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        
        echo "[GEP2] Running GenomeScope2 for reads-only profiling: {wildcards.species}"
        echo "[GEP2] K-mer length: {wildcards.kmer_len}"
        echo "[GEP2] Ploidy: {params.ploidy}"
        
        mkdir -p {params.outdir}
        
        genomescope2 -i {input.hist} \
                     -o {params.outdir} \
                     -k {wildcards.kmer_len} \
                     -p {params.ploidy} \
                     -n {wildcards.species}
        
        echo "[GEP2] ✅ GenomeScope2 complete"
        """