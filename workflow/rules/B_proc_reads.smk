# -------------------------------------------------------------------------------
# GEP2 - Read Processing Rules
# -------------------------------------------------------------------------------
#
# Routes per read type, with READS_PROC on (off: no QC, processing or Kraken2;
# every analysis uses the raw reads, and the base counts still run, on raw)
#
#   Illumina/10x/Hi-C  FastQC on raw reads (250k per mate)
#                      TRIM_PE: fastp -> processed/*_trimmed (10x: R1 barcode cut);
#                      fastp reports before/after stats
#   HiFi               FastQC (250k) + ULI primer check (100k) on raw reads
#                      FILTER_HIFI: bbduk drops reads with SMRTbell adapters
#                      -> processed/*_filtered; NanoPlot (250k) on the result
#   ONT                FastQC on raw reads (250k)
#                      CORRECT_ONT: hifiasm correction -> processed/*_corrected;
#                      NanoPlot (250k) on the result
#   (a flag off: no processing, NanoPlot on the raw reads)
#
# Later steps (k-mers, Inspector, Hi-C, mapping, base counts, Kraken2) use the
# processed file when its flag is on, the raw one otherwise (_get_processed_read_path).
# RUN_KRAKEN: Kraken2 on 50k reads/pairs of those files (G_read_screen.smk).
# Reports: reads/<type>/reports/, one MultiQC report per read type (B08).

# Note: The following are defined in the main Snakefile:
#   - _MISSING_IN: Placeholder for missing input files
#   - _mk_dirs(): Create log directories
#   - _is_10x(): Check if read type is 10x
#   - _get_centralized_reads_dict(): Get centralized read files
#   - _enumerate_centralized_groups(): Generate read groups
#   - _pick_centralized_group(): Pick read group by wildcards
#   - normalize_read_type(): Normalize read type names
#   - DOWNLOAD_MANIFEST, DOWNLOAD_MANIFEST_DICT: Download manifest data


# -------------------------------------------------------------------------------
# RULE ORDER
# -------------------------------------------------------------------------------

ruleorder: B04_check_uli_primers > B00_centralize_reads
ruleorder: B05_filter_hifi_adapters > B00_centralize_reads
ruleorder: B03_ont_correction > B00_centralize_reads
ruleorder: B02_trim_pe_fastp > B00_centralize_reads
ruleorder: B00_centralize_reads > B00_link_long_reads > B01_compress_long_reads
ruleorder: B00_centralize_reads > B00_link_pe_reads > B01_compress_pe_reads


# -------------------------------------------------------------------------------
# INPUT FUNCTIONS
# -------------------------------------------------------------------------------

# Load centralize map (at module level in B_proc_reads.smk, next to DOWNLOAD_MANIFEST)
_centralize_map_path = os.path.join(config["OUT_FOLDER"], "GEP2_results", "centralize_map.json")
if os.path.exists(_centralize_map_path):
    with open(_centralize_map_path) as f:
        CENTRALIZE_MAP = json.load(f)
else:
    CENTRALIZE_MAP = {}


def _get_source_for_centralized_read(wildcards):
    """Find the source file that needs to be symlinked to the centralized location."""
    centralized_path = os.path.join(
        config["OUT_FOLDER"], "GEP2_results", "data",
        wildcards.species, "reads", wildcards.read_type, wildcards.filename
    )

    if centralized_path in CENTRALIZE_MAP:
        return CENTRALIZE_MAP[centralized_path]

    raise ValueError(
        f"[GEP2] No source for centralized read:\n"
        f"  {centralized_path}\n"
        f"centralize_map.json has {len(CENTRALIZE_MAP)} entries, e.g.:\n  "
        + "\n  ".join(list(CENTRALIZE_MAP)[:3])
        + "\nA rule requested a read filename that centralization never produced."
    )


def _linkable_long_src(w):
    """Get source file for linking long reads."""
    g = _pick_centralized_group(w)
    yaml_path = str(g["r"])

    if not yaml_path.endswith(".gz"):
        return _MISSING_IN

    return _resolve_centralized_source(yaml_path, w.read_type, w.species)


def _compressible_long_src(w):
    """Get source file for compressing long reads."""
    g = _pick_centralized_group(w)
    return g["r"] if not str(g["r"]).endswith(".gz") else _MISSING_IN


def _linkable_pe_r1(w):
    """Get R1 source for linking PE reads."""
    if _as_bool(config.get("TRIM_PE", True)):
        return _MISSING_IN
    g = _pick_centralized_group(w)
    yaml_path = str(g["r1"])

    if not yaml_path.endswith(".gz"):
        return _MISSING_IN

    return _resolve_centralized_source(yaml_path, w.read_type, w.species)


def _linkable_pe_r2(w):
    """Get R2 source for linking PE reads."""
    if _as_bool(config.get("TRIM_PE", True)):
        return _MISSING_IN
    g = _pick_centralized_group(w)
    yaml_path = str(g["r2"])

    if not yaml_path.endswith(".gz"):
        return _MISSING_IN

    return _resolve_centralized_source(yaml_path, w.read_type, w.species)


def _resolve_centralized_source(yaml_path, read_type, species):
    """Given a centralized path under /data/, find the actual source file."""
    # If it's not a centralized path, return as-is
    if "/data/" not in yaml_path or "_Path" not in os.path.basename(yaml_path):
        return yaml_path

    # Look up in centralize map
    if yaml_path in CENTRALIZE_MAP:
        return CENTRALIZE_MAP[yaml_path]

    # Fallback: assume it was downloaded
    basename = os.path.basename(yaml_path)
    identifier = re.sub(r'^(hifi|ont|illumina|10x|hic)_Path\d+_', '', basename, flags=re.IGNORECASE)
    identifier = identifier.replace('.fq.gz', '').replace('.fastq.gz', '')

    parts = identifier.split("_")
    accession_parts = [p for p in parts if p.lower() not in ["hifi", "ont", "hic", "illumina", "10x"] and not p.startswith("Path")]
    accession = "_".join(accession_parts)

    downloaded_path = yaml_path.replace("/data/", "/downloaded_data/")
    downloaded_path = os.path.join(
        os.path.dirname(downloaded_path),
        f"{accession}.fastq.gz"
    )
    return downloaded_path


def _compressible_pe_r1(w):
    """Get R1 source for compressing PE reads."""
    if _as_bool(config.get("TRIM_PE", True)):
        return _MISSING_IN
    g = _pick_centralized_group(w)
    needs = (not str(g["r1"]).endswith(".gz")) or (not str(g["r2"]).endswith(".gz"))
    return g["r1"] if needs else _MISSING_IN


def _compressible_pe_r2(w):
    """Get R2 source for compressing PE reads."""
    if _as_bool(config.get("TRIM_PE", True)):
        return _MISSING_IN
    g = _pick_centralized_group(w)
    needs = (not str(g["r1"]).endswith(".gz")) or (not str(g["r2"]).endswith(".gz"))
    return g["r2"] if needs else _MISSING_IN


def _get_ont_reads_for_correction(w):
    """Get ONT reads for correction if enabled."""
    if not _as_bool(config.get("CORRECT_ONT", False)):
        return _MISSING_IN

    if w.read_type.lower() != "ont":
        return _MISSING_IN

    g = _pick_centralized_group(w)
    return g["r"]


def _get_uli_input(w):
    """Get HiFi reads for ULI primer check."""
    g = _pick_centralized_group(w)
    return g["r"]


def _get_hifi_input(w):
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

    # Anchored, not a substring search: a library legitimately named
    # e.g. "run_filtered_2019" contains "_filtered" but is a raw read set,
    # and would otherwise be looked for under processed/ and never found.
    # Processed long reads are always named "<base>_filtered" / "<base>_corrected".
    if w.filename.endswith(("_filtered", "_corrected")):
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
        "reads", read_type_lower, "reports"
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

                if _as_bool(config.get("FILTER_HIFI", True)):
                    reports.extend([
                        os.path.join(report_dir, f"{rt_lower}_Path{idx}_{base}_bbduk.stats"),
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

                if _as_bool(config.get("CORRECT_ONT", False)):
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

            if _as_bool(config.get("TRIM_PE", True)):
                reports.append(os.path.join(report_dir, f"{read_type_lower}_Path{idx}_{base}_fastp.json"))

    # Kraken2 read screen (READS_PROC + RUN_KRAKEN; G_read_screen.smk)
    reports.extend(get_kraken_multiqc_files(species, read_type_lower))

    # Read bases per file and the total (B10)
    reports.append(_read_bases_mqc_file(species, read_type_lower))

    return reports


# -------------------------------------------------------------------------------
# INPUT FUNCTIONS - Read base counts (B09/B10, used by the EAR yaml)
# -------------------------------------------------------------------------------
# Bases are counted on the reads the analyses use (_get_processed_read_path():
# processed when READS_PROC and the technology's flag are on, raw otherwise)

def _read_counts_dir(species, read_type):
    return os.path.join(
        config["OUT_FOLDER"], "GEP2_results", "data", species,
        "reads", read_type, "counts"
    )


def _read_count_file(species, read_type, read_file):
    """Per-file count table for one read file: counts/{file}.{raw|processed}.bases.tsv.
    The stage is part of the name, not guessed from a suffix like _trimmed, which a raw
    file could also carry."""
    stage = "processed" if os.path.basename(os.path.dirname(read_file)) == "processed" else "raw"
    name = os.path.basename(read_file)
    if name.endswith(".fq.gz"):
        name = name[:-len(".fq.gz")]
    return os.path.join(_read_counts_dir(species, read_type), f"{name}.{stage}.bases.tsv")


def _read_bases_summary_file(species, read_type):
    """Per-technology table: every counted file of the species plus the total."""
    return os.path.join(_read_counts_dir(species, read_type), f"{read_type}_bases_summary.tsv")


def _read_bases_mqc_file(species, read_type):
    """The same table as MultiQC custom content (B08), with the total in its description."""
    return os.path.join(_read_counts_dir(species, read_type), f"{read_type}_bases_mqc.yaml")


def _effective_read_entries_for_assembly(species, asm_id):
    """Read files used for one assembly, per normalized read type, in PathN order:
    {read_type: [(file, source), ...]}. Paired-end entries give both mates.
    'source' is the original file behind the centralized one (centralize_map.json),
    so one input file centralized under two names is recognised as the same file."""
    try:
        asm_data = samples_config["sp_name"][species]["asm_id"][asm_id] or {}
    except (KeyError, TypeError):
        return {}

    entries_by_type = {}
    for rt_key, rt_data in (asm_data.get("read_type") or {}).items():
        read_type = normalize_read_type(rt_key)
        if not read_type or not rt_data:
            continue
        paired = read_type in PE_READ_TYPES
        read_files = rt_data.get("read_files") or {}

        for path_key in sorted(read_files, key=_path_num):
            raw_paths = _split_entry(read_files[path_key])
            if not raw_paths:
                continue
            first = raw_paths[0]

            m = re.search(r'_Path(\d+)_', os.path.basename(first))
            if not m:
                raise ValueError(
                    f"[GEP2] Read base counts: {species}/{asm_id}/{read_type}/{path_key} is not "
                    f"a centralized read path ({first}); expected {read_type}_Path<N>_<name>.fq.gz"
                )
            base = read_base_from_path(first, paired)
            effective = _get_processed_read_path(species, read_type, m.group(1), base)
            effective = effective if isinstance(effective, list) else [effective]

            bucket = entries_by_type.setdefault(read_type, [])
            for i, f in enumerate(effective):
                raw = raw_paths[i] if i < len(raw_paths) else f
                source = os.path.normpath(str(CENTRALIZE_MAP.get(raw, raw)))
                if all(source != s for _, s in bucket):
                    bucket.append((f, source))

    return entries_by_type


def _species_read_files_by_source(species, read_type):
    """{source: file} over every entry of the species (reads-only entries included),
    in table order. An input file centralized under two names (e.g. assembly A lists
    reads a+b and assembly B only a, so 'a' becomes Path1 and Path3) keeps the first
    name, the same one the per-read k-mer databases are built from."""
    by_source = {}
    try:
        asm_ids = list(samples_config["sp_name"][species]["asm_id"])
    except (KeyError, TypeError):
        return by_source
    for asm_id in asm_ids:
        for f, source in _effective_read_entries_for_assembly(species, asm_id).get(read_type, []):
            by_source.setdefault(source, f)
    return by_source


def _effective_read_files_for_species(species, read_type):
    """Files counted for the species table: one per original input file."""
    return list(_species_read_files_by_source(species, read_type).values())


def _effective_read_files_for_assembly(species, asm_id):
    """{read_type: [file, ...]} counted for one assembly (its EAR DATA lines), using
    the species-level name of each input file so no file is counted twice."""
    files_by_type = {}
    for read_type, entries in _effective_read_entries_for_assembly(species, asm_id).items():
        canonical = _species_read_files_by_source(species, read_type)
        files_by_type[read_type] = [canonical.get(source, f) for f, source in entries]
    return files_by_type


def _get_count_bases_input(w):
    """B09: raw centralized file or its processed counterpart."""
    base_dir = os.path.join(
        config["OUT_FOLDER"], "GEP2_results", "data", w.species, "reads", w.read_type
    )
    if w.stage == "processed":
        base_dir = os.path.join(base_dir, "processed")
    return os.path.join(base_dir, f"{w.filename}.fq.gz")


def _get_bases_summary_inputs(w):
    """B10: per-file counts of every read file of this species/read type."""
    return [
        _read_count_file(w.species, w.read_type, f)
        for f in _effective_read_files_for_species(w.species, w.read_type)
    ]


# -------------------------------------------------------------------------------
# HELPERS - fastp adapters (B02, ADAPTERS_PE)
# -------------------------------------------------------------------------------

def _fastp_adapter_args():
    """fastp options for the configured ADAPTERS_PE (validated at load in Snakefile)."""
    if ADAPTERS_PE_MODE == "auto":
        return "--detect_adapter_for_pe"
    if ADAPTERS_PE_MODE == "overlap":
        return ""
    return f'--adapter_fasta "{ADAPTERS_PE_FASTA}"'


def _fastp_adapter_input(w):
    """The adapter FASTA as an input, so editing it re-trims (none for auto/overlap)."""
    return [ADAPTERS_PE_FASTA] if ADAPTERS_PE_MODE == "fasta" else []


# fastp prints a known adapter after its name line and a guessed one bare: the bare line is the signal
_FASTP_DETECT_RE = re.compile(r"^Detecting adapter sequence for (read[12])\.\.\.\s*$")


def _tandem_repeat_unit(seq, max_period=12, max_mismatch=0.1):
    """Repeat unit if seq is a short tandem repeat (telomere, satellite, ...), else None."""
    seq = seq.upper()
    for k in range(1, max_period + 1):
        n = len(seq) - k
        if n >= 2 * k and sum(seq[i] != seq[i + k] for i in range(n)) <= max_mismatch * n:
            return seq[:k]
    return None


def _fastp_unknown_adapter_guesses(text):
    """[(read, sequence)] for each adapter fastp guessed without matching its known list."""
    lines = text.splitlines()
    out = []
    for i, line in enumerate(lines):
        m = _FASTP_DETECT_RE.match(line.strip())
        if not m:
            continue
        nxt = lines[i + 1].strip() if i + 1 < len(lines) else ""
        if not nxt or nxt.startswith(">") or nxt.startswith("No adapter detected"):
            continue
        out.append((m.group(1), nxt))
    return out


# -------------------------------------------------------------------------------
# RULES - Read Centralization
# -------------------------------------------------------------------------------

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
            raise ValueError(
                f"[GEP2] Source file for centralization does not exist: {src}"
            )

        shell(cmd)


# -------------------------------------------------------------------------------
# RULES - Long Reads: Link/Compress
# -------------------------------------------------------------------------------

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
        read_type = r"(hifi|ont)",
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
        read_type = r"(hifi|ont)",
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


# -------------------------------------------------------------------------------
# RULES - Paired-End Reads: Link/Compress
# -------------------------------------------------------------------------------

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
        read_type = r"(illumina|10x|hic)",
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
        read_type = r"(illumina|10x|hic)",
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


# -------------------------------------------------------------------------------
# RULES - PE Trimming
# -------------------------------------------------------------------------------

rule B02_trim_pe_fastp:
    """Trim and QC paired-end reads with fastp (adapters: ADAPTERS_PE)."""
    input:
        r1 = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "{read_type}_Path{idx}_{base}_1.fq.gz"
        ),
        r2 = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "{read_type}_Path{idx}_{base}_2.fq.gz"
        ),
        adapters = _fastp_adapter_input
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
            "reads", "{read_type}", "reports", "{read_type}_Path{idx}_{base}_fastp.json"
        ),
        html = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "reports", "{read_type}_Path{idx}_{base}_fastp.html"
        )
    wildcard_constraints:
        read_type = r"(illumina|10x|hic)",
        idx = r"\d+",
        base = r"[^/]+"
    threads: cpu_func("trim_pe")
    resources:
        mem_mb = mem_func("trim_pe"),
        runtime = time_func("trim_pe")
    container: CONTAINERS["gep2_base"]
    params:
        extra_opts = lambda w: "--trim_front1 23 --trim_front2 0" if _is_10x(w.read_type) else "",
        adapter_mode = str(config.get("ADAPTERS_PE", "illumina")),
        adapter_args = _fastp_adapter_args()
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs",
            "B02_fastp.{read_type}_Path{idx}_{base}.log"
        )
    run:
        _mk_dirs(log[0])
        os.makedirs(os.path.dirname(output.r1), exist_ok=True)
        os.makedirs(os.path.dirname(output.json), exist_ok=True)

        with open(log[0], "w") as fh:
            fh.write(f"[GEP2] ADAPTERS_PE: {params.adapter_mode}\n")

        shell(
            'fastp -i "{input.r1}" -I "{input.r2}" '
            '-o "{output.r1}" -O "{output.r2}" '
            '--thread {threads} '
            '--json "{output.json}" --html "{output.html}" '
            '{params.adapter_args} '
            '--length_required 50 '
            '--overrepresentation_analysis '
            '{params.extra_opts} '
            '>> "{log}" 2>&1'
        )

        # auto: fastp may take an over-represented genome sequence (telomere, satellite)
        # for an adapter and trim it from every read that contains it
        if ADAPTERS_PE_MODE == "auto":
            with open(log[0]) as fh:
                guesses = _fastp_unknown_adapter_guesses(fh.read())
            for read, seq in guesses:
                unit = _tandem_repeat_unit(seq)
                hint = (f" It is a short tandem repeat (unit {unit}), e.g. a telomere or "
                        f"satellite repeat." if unit else "")
                msg = (f"[GEP2] WARNING: {wildcards.species}/{wildcards.read_type} "
                       f"Path{wildcards.idx}_{wildcards.base}: fastp guessed the {read} adapter "
                       f"as {seq[:30]}..., which is not a known adapter. It may be genome "
                       f"sequence, now trimmed from every read that contains it.{hint} Check "
                       f"the fastp report, or set ADAPTERS_PE: illumina and re-trim these reads.")
                with open(log[0], "a") as fh:
                    fh.write(msg + "\n")
                print(msg)

        with open(log[0], "a") as fh:
            fh.write(f"[GEP2] fastp completed for {wildcards.species}/{wildcards.read_type}\n")


# -------------------------------------------------------------------------------
# RULES - Long Read Processing
# -------------------------------------------------------------------------------

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

        WORK_DIR="$(gep2_get_workdir 250)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GEP2_ont_correct_{wildcards.species}_{wildcards.read_type}_{wildcards.base}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT

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

        '''


rule B04_check_uli_primers:
    """Check for ULI primers in HiFi reads (100k subsample)."""
    input:
        reads = _get_uli_input
    output:
        report = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "reports",
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

        WORK_DIR="$(gep2_get_workdir 5)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GEP2_uli_check_XXXXXX")"
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
    """Filter HiFi reads containing SMRTbell adapter sequences using bbduk."""
    input:
        reads = _get_hifi_input,
        adapters = workflow.source_path("../scripts/hifi_blocklist.fa")
    output:
        filtered = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "processed", "{read_type}_Path{idx}_{base}_filtered.fq.gz"
        ),
        bbduk_stats = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "reports",
            "{read_type}_Path{idx}_{base}_bbduk.stats"
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
            "B05_bbduk.{read_type}_Path{idx}_{base}.log"
        )
    shell:
        r'''
        set -euo pipefail
        exec > {log} 2>&1

        mkdir -p $(dirname {output.filtered})
        mkdir -p $(dirname {output.bbduk_stats})

        bbduk.sh -Xmx$(({resources.mem_mb} * 85 / 100))m \
            in={input.reads} \
            out={output.filtered} \
            ref={input.adapters} \
            stats={output.bbduk_stats} \
            k=27 hdist=1 rcomp=t \
            threads={threads} \
            overwrite=t
        '''


# -------------------------------------------------------------------------------
# RULES - QC
# -------------------------------------------------------------------------------

rule B06_fastqc_pe_reads:
    """Run FastQC on paired-end reads (250k subsample)."""
    input:
        reads = _get_fastqc_pe_input
    output:
        html = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "reports",
            "{read_type}_Path{idx}_{base}_{read}_fastqc.html"
        ),
        zip = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "reports",
            "{read_type}_Path{idx}_{base}_{read}_fastqc.zip"
        )
    wildcard_constraints:
        read_type = r"(illumina|10x|hic)",
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

        WORK_DIR="$(gep2_get_workdir 10)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GEP2_fastqc_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT

        cd "$TEMP_DIR"

        export _JAVA_OPTIONS="-Xmx$(({resources.mem_mb} * 90 / 100))m"

        TEMP_SUB="{wildcards.read_type}_Path{wildcards.idx}_{wildcards.base}_{wildcards.read}.fq"
        seqtk sample -s123 {input.reads} 250000 > $TEMP_SUB

        fastqc -t {threads} --outdir . $TEMP_SUB

        mv "{wildcards.read_type}_Path{wildcards.idx}_{wildcards.base}_{wildcards.read}_fastqc.html" "{output.html}"
        mv "{wildcards.read_type}_Path{wildcards.idx}_{wildcards.base}_{wildcards.read}_fastqc.zip" "{output.zip}"

        """


rule B06_fastqc_long_reads:
    """Run FastQC on long reads (250k subsample)."""
    input:
        reads = _get_fastqc_long_input
    output:
        html = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "reports",
            "{read_type}_Path{idx}_{base}_fastqc.html"
        ),
        zip = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "reports",
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

        WORK_DIR="$(gep2_get_workdir 10)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GEP2_fastqc_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT

        cd "$TEMP_DIR"

        export _JAVA_OPTIONS="-Xmx$(({resources.mem_mb} * 90 / 100))m"

        TEMP_SUB="{wildcards.read_type}_Path{wildcards.idx}_{wildcards.base}.fq"
        seqtk sample -s123 {input.reads} 250000 > $TEMP_SUB

        fastqc -t {threads} --outdir . $TEMP_SUB

        mv "{wildcards.read_type}_Path{wildcards.idx}_{wildcards.base}_fastqc.html" "{output.html}"
        mv "{wildcards.read_type}_Path{wildcards.idx}_{wildcards.base}_fastqc.zip" "{output.zip}"

        """


rule B07_nanoplot_long_reads:
    """Run NanoPlot QC on long reads."""
    input:
        reads = _get_nanoplot_input
    output:
        dir = directory(os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "reports",
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

        WORK_DIR="$(gep2_get_workdir 10)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GEP2_nanoplot_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT

        cd $TEMP_DIR

        echo "[GEP2] NanoPlot runs on a random subsample of up to 250000 reads:"
        echo "[GEP2] number of reads/bases in NanoStats.txt (and MultiQC) refer to the subsample"
        seqtk sample -s456 {input.reads} 250000 | gzip -c > subsample.fq.gz

        if NanoPlot -t {threads} --fastq subsample.fq.gz -o {output.dir} --tsv_stats; then
            echo "[GEP2] NanoPlot completed successfully"
        elif [ -s {output.dir}/NanoStats.txt ]; then
            echo "[GEP2] WARNING: NanoPlot failed after writing NanoStats.txt (kept); some plots may be missing"
        else
            echo "[GEP2] WARNING: NanoPlot failed before writing stats; writing a placeholder NanoStats.txt with zeros"
            cat > {output.dir}/NanoStats.txt << 'EOF'
General summary:
Number of reads:	0
Number of bases:	0
Median read length:	0.0
Mean read length:	0.0
Read length N50:	0
EOF
        fi

        """


# -------------------------------------------------------------------------------
# RULES - MultiQC Aggregation
# -------------------------------------------------------------------------------

rule B08_aggregate_read_qc:
    """Aggregate all QC reports with MultiQC."""
    input:
        _get_qc_reports_for_multiqc
    output:
        html = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "reports",
            "multiqc_report.html"
        ),
        data = directory(os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "reports",
            "multiqc_report_data"
        ))
    threads: cpu_func("reads_qc")
    resources:
        mem_mb = mem_func("reads_qc"),
        runtime = time_func("reads_qc")
    params:
        # section order and notes (+ R1/R2 rows under their pair for paired-end types)
        mqc_config = lambda w: " ".join(
            f'-c "{SCRIPTS_DIR / f}"'
            for f in ["multiqc_config.yaml"]
                     + (["multiqc_config_pe.yaml"] if w.read_type.lower() in PE_READ_TYPES else [])
        )
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

        multiqc -f -n multiqc_report {params.mqc_config} \
            -o $(dirname {output.html}) \
            --file-list $TEMP_LIST

        rm $TEMP_LIST
        """


# -------------------------------------------------------------------------------
# RULES - Read Base Counts
# -------------------------------------------------------------------------------

rule B09_count_read_bases:
    """Count reads and bases of one read file with seqtk size.
    {stage} is 'raw' (centralized file) or 'processed' (trimmed/filtered/corrected)."""
    input:
        reads = _get_count_bases_input
    output:
        counts = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "counts", "{filename}.{stage}.bases.tsv"
        )
    wildcard_constraints:
        read_type = r"[^/]+",
        stage = r"(raw|processed)",
        filename = r"[^/]+"
    threads: cpu_func("count_bases")
    resources:
        mem_mb = mem_func("count_bases"),
        runtime = time_func("count_bases")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs",
            "B09_count_read_bases.{stage}.{filename}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1

        mkdir -p "$(dirname {output.counts})"

        # seqtk size: "<reads>\t<bases>"; a truncated gzip is counted up to the cut
        COUNTS=$(seqtk size {input.reads})
        echo "[GEP2] seqtk size {input.reads}: $COUNTS"

        if ! printf '%s\\n' "$COUNTS" | grep -Eq '^[0-9]+[[:space:]]+[0-9]+$'; then
            echo "[GEP2] ERROR: unexpected seqtk size output: '$COUNTS'"
            exit 1
        fi
        read -r N_READS N_BASES <<< "$COUNTS"
        if [ "$N_BASES" -eq 0 ]; then
            echo "[GEP2] WARNING: no bases in {input.reads}"
        fi

        printf 'file\\tstage\\treads\\tbases\\n%s\\t%s\\t%s\\t%s\\n' \
            "$(basename {input.reads})" "{wildcards.stage}" "$N_READS" "$N_BASES" \
            > {output.counts}
        """


rule B10_summarise_read_bases:
    """Per species and technology: one row per counted read file, plus the total, as a
    TSV and as MultiQC custom content (B08). The EAR yaml adds up only the files of each
    assembly (see Z_reporting.smk)."""
    input:
        counts = _get_bases_summary_inputs
    output:
        summary = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "counts", "{read_type}_bases_summary.tsv"
        ),
        mqc = os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "counts", "{read_type}_bases_mqc.yaml"
        )
    wildcard_constraints:
        read_type = r"[^/]+"
    threads: 1
    resources:
        mem_mb = 1000,
        runtime = 15
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs", "B10_summarise_read_bases.log"
        )
    run:
        _mk_dirs(log[0])

        rows = []
        for path in input.counts:
            with open(path) as fh:
                header = fh.readline().rstrip("\n").split("\t")
                values = fh.readline().rstrip("\n").split("\t")
            rec = dict(zip(header, values))
            rows.append((rec["file"], rec["stage"], int(rec["reads"]), int(rec["bases"])))

        total_reads = sum(r[2] for r in rows)
        total_bases = sum(r[3] for r in rows)
        stages = "/".join(sorted({r[1] for r in rows}))

        with open(output.summary, "w") as out:
            out.write("file\tstage\treads\tbases\tGbp\n")
            for name, stage, n_reads, n_bases in rows:
                out.write(f"{name}\t{stage}\t{n_reads}\t{n_bases}\t{n_bases / 1e9:.3f}\n")
            out.write(f"TOTAL\t{stages}\t{total_reads}\t{total_bases}\t{total_bases / 1e9:.3f}\n")

        # MultiQC table of the per-file rows; the total goes in the description (rows are sorted by name)
        mqc = {
            "id": f"gep2_read_bases_{wildcards.read_type}",
            "section_name": "Read bases",
            "description": (
                f"Reads and bases per file (seqtk size), counted on the {stages} reads; "
                f"paired-end mates are counted separately. <b>Total: {total_reads:,} reads, "
                f"{total_bases / 1e9:,.2f} Gbp ({total_bases:,} bases)</b>"
            ),
            "plot_type": "table",
            "pconfig": {
                "id": f"gep2_read_bases_{wildcards.read_type}_table",
                "title": "GEP2: read bases",
            },
            "headers": {
                "stage": {"title": "Stage", "description": "raw or processed reads"},
                "reads": {"title": "Reads", "description": "Number of reads",
                          "format": "{:,.0f}", "scale": False},
                "bases": {"title": "Bases", "description": "Number of bases",
                          "format": "{:,.0f}", "scale": False},
                "gbp": {"title": "Gbp", "description": "Gigabases", "format": "{:,.3f}",
                        "scale": False},
            },
            "data": {
                name: {"stage": stage, "reads": n_reads, "bases": n_bases,
                       "gbp": round(n_bases / 1e9, 3)}
                for name, stage, n_reads, n_bases in rows
            },
        }
        with open(output.mqc, "w") as out:
            yaml.safe_dump(mqc, out, sort_keys=False)

        with open(log[0], "w") as lf:
            lf.write(f"[GEP2] {len(rows)} file(s), {total_bases} bases ({total_bases / 1e9:.3f} Gbp)\n")
