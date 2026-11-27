# -------------------------------------------------------------------------------
# GEP2 - Read Screening (Kraken2 core_nt)
# -------------------------------------------------------------------------------
#
# Kraken2 core_nt screen of a random subset of every read set (READS_PROC + RUN_KRAKEN),
# shown in the read-QC MultiQC report (B08). Read sets: the files the analyses use;
# a file listed twice is screened once.
#
#   G00  core_nt DB    -> DB_FOLDER/kraken2/core_nt (pinned build, md5-checked, once)
#   G01  seqtk sample  -> 50k reads / 50k pairs per read set
#   G02  Kraken2       -> all subsets in one job (DB loaded once)
#   G03  target order  -> from the DB's own taxonomy
#   G04  summary       -> MultiQC General Stats + table, warnings
# -------------------------------------------------------------------------------

# Pinned core_nt build (https://benlangmead.github.io/aws-indexes/k2) and the md5 sums
# published with it (core_nt_20251015/core_nt.md5)
KRAKEN_DB_BUILD = "20251015"
KRAKEN_DB_URL = f"https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_{KRAKEN_DB_BUILD}.tar.gz"
KRAKEN_DB_ARCHIVE = f"k2_core_nt_{KRAKEN_DB_BUILD}.tar.gz"
KRAKEN_DB_ARCHIVE_MD5 = "0ae08edc9b7b351959a5af6f3ad817f4"
KRAKEN_DB_ARCHIVE_BYTES = 263_961_481_154
KRAKEN_DB_FILE_MD5 = {                  # checked after unpacking / before reuse
    "hash.k2d": "beda141de277d59907a26c65e1411137",
    "opts.k2d": "e94ea0e4b91a8dac829a7dcd1a0e4ba2",
    "taxo.k2d": "4d7ce370780e447e417ee13b5dfd6f95",
    "names.dmp": "684a5b76ddf23d9b8ff7041e44cb45fe",
    "nodes.dmp": "01ec68e720c35030af6777fba6e48ec8",
}
# Unpacked: the files above + two small extras (ktaxonomy.tsv, inspect.txt). Not kept:
# library_report.tsv (22 GB), seqid2taxid.map (5 GB), Bracken *mers.kmer_distrib.
KRAKEN_DB_KEEP = list(KRAKEN_DB_FILE_MD5) + ["ktaxonomy.tsv", "inspect.txt"]
KRAKEN_DB_KEPT_BYTES = 341_000_000_000
KRAKEN_DB_FLAG = os.path.join(KRAKEN_DB_DIR, "core_nt.done")

KRAKEN_SUBSET = 50_000          # reads, or read pairs, per read set
KRAKEN_SEED = 42
KRAKEN_CONFIDENCE = 0.1
KRAKEN_READ_TYPES = ("hifi", "ont", "illumina", "10x", "hic")


def _kraken_enabled():
    return (_as_bool(config.get("READS_PROC", True))
            and _as_bool(config.get("RUN_KRAKEN", False)))


def _kraken_dir(species, read_type):
    return os.path.join(config["OUT_FOLDER"], "GEP2_results", "data", species,
                        "reads", read_type, "kraken2")


def _kraken_read_sets(species):
    """Read sets of a species for the screen: [{'group', 'read_type', 'files', 'paired'}],
    one per long-read file or PE library, in table order. 'group' is
    <read_type>_Path<N>_<base> (also the MultiQC sample name). A set whose original
    files were already listed (same files under another PathN) is skipped."""
    sets, seen = [], set()
    try:
        asm_ids = list(samples_config["sp_name"][species]["asm_id"])
    except (KeyError, TypeError):
        return sets

    for asm_id in asm_ids:
        asm_data = samples_config["sp_name"][species]["asm_id"][asm_id] or {}
        for rt_key, rt_data in (asm_data.get("read_type") or {}).items():
            read_type = normalize_read_type(rt_key)
            if read_type not in KRAKEN_READ_TYPES or not rt_data:
                continue
            paired = read_type in PE_READ_TYPES
            read_files = rt_data.get("read_files") or {}

            for path_key in sorted(read_files, key=_path_num):
                raw_paths = _split_entry(read_files[path_key])
                if not raw_paths or (paired and len(raw_paths) < 2):
                    continue
                m = re.search(r'_Path(\d+)_', os.path.basename(raw_paths[0]))
                if not m:
                    continue
                base = read_base_from_path(raw_paths[0], paired)
                files = _get_processed_read_path(species, read_type, m.group(1), base)
                files = files if isinstance(files, list) else [files]

                sources = tuple(os.path.normpath(str(CENTRALIZE_MAP.get(p, p)))
                                for p in raw_paths[:len(files)])
                if (read_type, sources) in seen:
                    continue
                seen.add((read_type, sources))
                sets.append({"group": f"{read_type}_Path{m.group(1)}_{base}",
                             "read_type": read_type, "files": files, "paired": paired})
    return sets


def _kraken_subset_files(species, read_set):
    d = _kraken_dir(species, read_set["read_type"])
    if read_set["paired"]:
        return [os.path.join(d, f"{read_set['group']}.sub_1.fq.gz"),
                os.path.join(d, f"{read_set['group']}.sub_2.fq.gz")]
    return [os.path.join(d, f"{read_set['group']}.sub.fq.gz")]


def _kraken_report_file(species, read_set):
    return os.path.join(_kraken_dir(species, read_set["read_type"]),
                        f"{read_set['group']}.kraken2.report")


def _kraken_target_file(species):
    return os.path.join(config["OUT_FOLDER"], "GEP2_results", "data", species,
                        "reads", "kraken2_target.json")


def get_kraken_multiqc_files(species, read_type):
    """B08 inputs: Kraken2 report (MultiQC's Kraken module) + GEP2 summaries (custom
    content) of each read set of this species and read type."""
    if not _kraken_enabled():
        return []
    files = []
    for read_set in _kraken_read_sets(species):
        if read_set["read_type"] != read_type.lower():
            continue
        d = _kraken_dir(species, read_set["read_type"])
        files += [_kraken_report_file(species, read_set),
                  os.path.join(d, f"{read_set['group']}_kraken2_mqc.yaml"),
                  os.path.join(d, f"{read_set['group']}_kraken2_table_mqc.yaml")]
    return files


def _kraken_find_set(w):
    for read_set in _kraken_read_sets(w.species):
        if read_set["read_type"] == w.read_type and read_set["group"] == w.group:
            return read_set
    raise ValueError(f"[GEP2] Kraken2: no read set {w.group} for {w.species}")


# All subsets and reports of the run, known at parse time: one G02 job classifies them
# all, so the 340 GB database is loaded once.
KRAKEN_JOBS = []
if _kraken_enabled() and samples_config and "sp_name" in samples_config:
    for _sp in samples_config["sp_name"]:
        for _rs in _kraken_read_sets(_sp):
            KRAKEN_JOBS.append({
                "species": _sp, "group": _rs["group"], "paired": _rs["paired"],
                "subsets": _kraken_subset_files(_sp, _rs),
                "report": _kraken_report_file(_sp, _rs),
            })


# -------------------------------------------------------------------------------
# RULES
# -------------------------------------------------------------------------------

rule G00_kraken_db:
    """core_nt into DB_FOLDER once (skipped if already there with matching md5s): space
    check, resumable aria2c download (12 h per attempt), md5, unpack, md5, delete archive."""
    output:
        flag = KRAKEN_DB_FLAG
    params:
        db_dir = KRAKEN_DB_DIR,
        url = KRAKEN_DB_URL,
        archive = KRAKEN_DB_ARCHIVE,
        archive_md5 = KRAKEN_DB_ARCHIVE_MD5,
        archive_bytes = KRAKEN_DB_ARCHIVE_BYTES,
        kept_bytes = KRAKEN_DB_KEPT_BYTES,
        md5_list = " ".join(f"{f}:{m}" for f, m in KRAKEN_DB_FILE_MD5.items()),
        keep = " ".join(KRAKEN_DB_KEEP),
        build = KRAKEN_DB_BUILD,
        watchdog = WATCHDOG_SCRIPT
    threads: cpu_func("kraken_db_download")
    resources:
        mem_mb = mem_func("kraken_db_download"),
        runtime = time_func("kraken_db_download")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(KRAKEN_DB_DIR, "kraken_db.log")
    shell:
        """
        set -euo pipefail
        mkdir -p {params.db_dir}
        exec > {log} 2>&1
        source {params.watchdog}
        FLAG="$(realpath -m {output.flag})"
        cd {params.db_dir}

        files_ok() {{
            local item f m
            for item in {params.md5_list}; do
                f="${{item%%:*}}"; m="${{item##*:}}"
                [ -s "$f" ] || {{ echo "[GEP2]   $f missing"; return 1; }}
                echo "[GEP2]   md5 $f ..."
                [ "$(md5sum "$f" | cut -d' ' -f1)" = "$m" ] || {{ echo "[GEP2]   $f: md5 mismatch"; return 1; }}
            done
        }}

        echo "[GEP2] Kraken2 core_nt {params.build} in {params.db_dir}"
        if files_ok; then
            echo "[GEP2] Database already present and intact: reusing it"
        else
            # ---- download (resumable) ----
            if [ ! -s "{params.archive}" ] || [ -e "{params.archive}.aria2" ]; then
                HAVE=$( [ -e "{params.archive}" ] && stat -c%s "{params.archive}" || echo 0 )
                NEED=$(( {params.archive_bytes} - HAVE + {params.kept_bytes} + 5000000000 ))
                AVAIL=$(( $(df -k . | awk 'NR==2 {{print $4}}') * 1024 ))
                if [ "$AVAIL" -lt "$NEED" ]; then
                    echo "[GEP2] ERROR: $(( AVAIL / 1000000000 )) GB free in {params.db_dir}," \
                         "need $(( NEED / 1000000000 )) GB (archive + unpacked database)"
                    exit 1
                fi
                ARIA_EXIT=1
                for attempt in 1 2 3; do
                    echo "[GEP2] Download attempt $attempt: {params.url}"
                    ARIA_EXIT=0
                    gep2_download_with_timeout 43200 aria2c \
                        -x 16 -s 16 -k 100M -c \
                        --max-tries=3 --retry-wait=10 \
                        --allow-overwrite=true \
                        --auto-file-renaming=false \
                        --file-allocation=none \
                        --console-log-level=warn --summary-interval=300 \
                        -d . -o "{params.archive}" "{params.url}" || ARIA_EXIT=$?
                    [ "$ARIA_EXIT" -eq 0 ] && break
                    echo "[GEP2] aria2c exit $ARIA_EXIT (a rerun resumes the download)"
                done
                [ "$ARIA_EXIT" -eq 0 ] || exit "$ARIA_EXIT"
            fi

            echo "[GEP2] Checking the archive md5 ..."
            if [ "$(md5sum "{params.archive}" | cut -d' ' -f1)" != "{params.archive_md5}" ]; then
                echo "[GEP2] ERROR: archive md5 mismatch: removing it, rerun to download again"
                rm -f "{params.archive}" "{params.archive}.aria2"
                exit 1
            fi

            # ---- unpack only the files kept ----
            HAVE_KEPT=$( (du -kc {params.keep} 2>/dev/null || true) | awk 'END {{print $1 * 1024}}')   # du fails for files not there yet
            AVAIL=$(( $(df -k . | awk 'NR==2 {{print $4}}') * 1024 ))
            if [ $(( AVAIL + HAVE_KEPT )) -lt {params.kept_bytes} ]; then
                echo "[GEP2] ERROR: $(( AVAIL / 1000000000 )) GB free in {params.db_dir}," \
                     "need $(( {params.kept_bytes} / 1000000000 )) GB to unpack"
                exit 1
            fi
            echo "[GEP2] Unpacking: {params.keep}"
            tar -I pigz -xf "{params.archive}" -C . --wildcards --no-anchored {params.keep}
            for f in {params.keep}; do          # members may sit in a sub-folder
                if [ ! -e "$f" ]; then
                    found=$(find . -mindepth 2 -type f -name "$f" | head -1)
                    [ -n "$found" ] && mv "$found" "./$f"
                fi
            done
            find . -mindepth 1 -type d -empty -delete

            echo "[GEP2] Checking the unpacked files ..."
            if ! files_ok; then
                echo "[GEP2] ERROR: unpacked files do not match the published md5 sums"
                exit 1
            fi
            rm -f "{params.archive}"
            echo "[GEP2] Archive removed"
        fi

        printf 'core_nt %s\\nverified %s\\n' "{params.build}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$FLAG"
        echo "[GEP2] Kraken2 core_nt ready"
        """


rule G01_kraken_subsample:
    """Random subset of one read set (seqtk sample, 2-pass: only the chosen record
    indices are kept in memory). Paired files get the same seed, so the same pairs are
    kept in both; checked by comparing all read names. A smaller set is taken whole."""
    input:
        reads = lambda w: _kraken_find_set(w)["files"][1 if w.sub == "sub_2" else 0]
    output:
        subsets = temp(os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "kraken2", "{group}.{sub}.fq.gz"
        ))
    wildcard_constraints:
        read_type = "|".join(KRAKEN_READ_TYPES),
        group = r"[^/]+_Path\d+_[^/]+",
        sub = r"sub(_[12])?"
    params:
        n = KRAKEN_SUBSET,
        seed = KRAKEN_SEED
    threads: cpu_func("subsample_reads")
    resources:
        mem_mb = mem_func("subsample_reads"),
        runtime = time_func("subsample_reads")
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
            "reads", "{read_type}", "logs", "G01_kraken_subsample.{group}.{sub}.log"
        )
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        echo "[GEP2] {wildcards.group}: random subset of {params.n} reads (seed {params.seed}) from {input.reads}"
        seqtk sample -2 -s{params.seed} {input.reads} {params.n} | pigz -p {threads} > {output.subsets}
        echo "[GEP2] Reads in subset: $(seqtk size {output.subsets} | cut -f1)"
        """


rule G01b_kraken_pair_check:
    """R1 and R2 subsets of a PE library must hold the same pairs in the same order."""
    input:
        r1 = os.path.join(config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
                          "reads", "{read_type}", "kraken2", "{group}.sub_1.fq.gz"),
        r2 = os.path.join(config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
                          "reads", "{read_type}", "kraken2", "{group}.sub_2.fq.gz")
    output:
        check = temp(os.path.join(config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
                                  "reads", "{read_type}", "kraken2", "{group}.pairs_ok"))
    wildcard_constraints:
        read_type = "illumina|10x|hic",
        group = r"[^/]+_Path\d+_[^/]+"
    threads: 1
    resources:
        mem_mb = 1000,
        runtime = 30
    container: CONTAINERS["gep2_base"]
    log:
        os.path.join(config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
                     "reads", "{read_type}", "logs", "G01b_kraken_pair_check.{group}.log")
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1
        names() {{ seqtk seq "$1" | awk 'NR%4==1 {{n = substr($1, 2); sub(/\\/[12]$/, "", n); print n}}'; }}
        read -r N_PAIRS N_BAD < <(paste <(names {input.r1}) <(names {input.r2}) \\
            | awk '{{n++; if ($1 != $2) bad++}} END {{print n+0, bad+0}}')
        echo "[GEP2] {wildcards.group}: $N_PAIRS pairs, $N_BAD R1/R2 name mismatches"
        if [ "$N_PAIRS" -eq 0 ]; then
            echo "[GEP2] ERROR: empty subset"; exit 1
        fi
        if [ "$N_BAD" -gt 0 ]; then
            echo "[GEP2] WARNING: R1 and R2 of {wildcards.group} are not paired record by record" \\
                 "($N_BAD of $N_PAIRS names differ): check the input files"
        fi
        touch {output.check}
        """


if KRAKEN_JOBS:
    rule G02_kraken_classify:
        """Classify every subset with Kraken2 against core_nt in one job, so the database
        is loaded once. If FAST_TEMP_DIR (e.g. /dev/shm) has room, the three .k2d files
        are copied there; otherwise the hash is read once into the page cache. Kraken2
        then memory-maps it for each subset."""
        input:
            db_flag = KRAKEN_DB_FLAG,
            subsets = [s for job in KRAKEN_JOBS for s in job["subsets"]],
            pair_checks = [
                os.path.join(os.path.dirname(job["report"]), f"{job['group']}.pairs_ok")
                for job in KRAKEN_JOBS if job["paired"]
            ]
        output:
            reports = [job["report"] for job in KRAKEN_JOBS]
        params:
            db_dir = KRAKEN_DB_DIR,
            confidence = KRAKEN_CONFIDENCE,
            min_mem_mb = 345_000
        threads: cpu_func("kraken_classify")
        resources:
            mem_mb = mem_func("kraken_classify"),
            runtime = time_func("kraken_classify")
        container: CONTAINERS["kraken"]
        log:
            os.path.join(config["OUT_FOLDER"], "GEP2_results", "logs", "G02_kraken_classify.log")
        run:
            _mk_dirs(log[0])
            staged_note = log[0] + ".staged"
            n_sets = len(KRAKEN_JOBS)
            try:
                shell(
                    """
                    exec >> {log} 2>&1
                    echo "[GEP2] Kraken2 core_nt screen of {n_sets} read set(s)"
                    MEM_KB=$(awk '/^MemTotal:/ {{print $2}}' /proc/meminfo)
                    if [ {resources.mem_mb} -lt {params.min_mem_mb} ] || [ $(( MEM_KB / 1024 )) -lt {params.min_mem_mb} ]; then
                        echo "[GEP2] ERROR: Kraken2 core_nt needs ~360 GB RAM (job: {resources.mem_mb} MB," \
                             "machine: $(( MEM_KB / 1024 )) MB); raise myRAM or set RUN_KRAKEN: off"
                        exit 1
                    fi
                    DB="{params.db_dir}"
                    FAST_TMP="${{GEP2_FAST_TMP:-}}"
                    NEED_KB=$(du -kL "$DB/hash.k2d" "$DB/opts.k2d" "$DB/taxo.k2d" | awk '{{s += $1}} END {{print s + 5242880}}')
                    if [ -n "$FAST_TMP" ] && [ "$FAST_TMP" != "${{GEP2_TMP:-}}" ] && mkdir -p "$FAST_TMP" 2>/dev/null \
                       && [ "$(df -k "$FAST_TMP" | awk 'NR==2 {{print $4}}')" -ge "$NEED_KB" ]; then
                        STAGED="$(mktemp -d "$FAST_TMP/GEP2_kraken2_core_nt_XXXXXX")"
                        echo "$STAGED" > "{staged_note}"      # cleaned up even if the copy fails
                        echo "[GEP2] Copying the database to $STAGED"
                        cp "$DB/hash.k2d" "$DB/opts.k2d" "$DB/taxo.k2d" "$STAGED/"
                    else
                        echo "[GEP2] FAST_TEMP_DIR not set or too small: reading the hash into the page cache"
                        cat "$DB/hash.k2d" > /dev/null
                    fi
                    """
                )
                db = params.db_dir
                if os.path.exists(staged_note):
                    with open(staged_note) as fh:
                        db = fh.read().strip()

                for job in KRAKEN_JOBS:
                    job_name = f"{job['species']}/{job['group']}"
                    paired = "--paired" if job["paired"] else ""
                    files = " ".join(f'"{s}"' for s in job["subsets"])
                    report = job["report"]
                    shell(
                        'echo "[GEP2] {job_name}" >> {log}; '
                        'kraken2 --db "{db}" --memory-mapping --threads {threads} '
                        '--confidence {params.confidence} {paired} '
                        '--report "{report}.tmp" --output /dev/null {files} >> {log} 2>&1; '
                        'mv "{report}.tmp" "{report}"'
                    )
            finally:
                if os.path.exists(staged_note):
                    with open(staged_note) as fh:
                        staged = fh.read().strip()
                    if staged and os.path.basename(staged).startswith("GEP2_kraken2_core_nt_"):
                        shell('rm -rf "{staged}"; echo "[GEP2] Removed the RAM copy {staged}" >> {log}')
                    os.remove(staged_note)


rule G03_kraken_target:
    """Target lineage of a species from the database's own taxonomy (same tax IDs as the
    reports): scientific name from the table key (Genus_species), then the first two
    words, then the genus; its order (class or phylum when a lineage has no order)."""
    input:
        db_flag = KRAKEN_DB_FLAG
    output:
        target = os.path.join(config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
                              "reads", "kraken2_target.json")
    params:
        names = os.path.join(KRAKEN_DB_DIR, "names.dmp"),
        nodes = os.path.join(KRAKEN_DB_DIR, "nodes.dmp")
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 60
    log:
        os.path.join(config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
                     "reads", "logs", "G03_kraken_target.log")
    run:
        import json
        _mk_dirs(log[0])
        words = wildcards.species.replace("_", " ").split()
        candidates = [" ".join(words), " ".join(words[:2]), words[0]] if words else []
        candidates = list(dict.fromkeys(c for c in candidates if c))

        found = {}
        with open(params.names) as fh:
            for line in fh:
                parts = line.split("\t|\t")
                if len(parts) >= 4 and parts[1] in candidates and parts[3].startswith("scientific name"):
                    found.setdefault(parts[1], int(parts[0]))
        query = next((c for c in candidates if c in found), None)

        target = {"species_key": wildcards.species, "query": query, "taxid": None,
                  "level_taxid": None, "level_name": None, "level_rank": None, "lineage": []}
        if query is not None:
            parent, rank, name = {}, {}, {}
            with open(params.nodes) as fh:
                for line in fh:
                    p = line.split("\t|\t")
                    parent[int(p[0])], rank[int(p[0])] = int(p[1]), p[2]
            lineage, t = [], found[query]
            while t not in lineage:
                lineage.append(t)
                if t == 1 or t not in parent:
                    break
                t = parent[t]
            level = next((x for r in ("order", "class", "phylum")
                          for x in lineage if rank.get(x) == r), None)
            target.update({"taxid": found[query], "lineage": lineage,
                           "level_taxid": level, "level_rank": rank.get(level)})
            if level is not None:
                with open(params.names) as fh:
                    for line in fh:
                        parts = line.split("\t|\t")
                        if int(parts[0]) == level and parts[3].startswith("scientific name"):
                            target["level_name"] = parts[1]
                            break

        with open(output.target, "w") as fh:
            json.dump(target, fh, indent=1)
        with open(log[0], "w") as lf:
            if query is None:
                lf.write(f"[GEP2] WARNING: none of {candidates} is in the core_nt taxonomy: "
                         f"no target lineage, reads are only grouped by domain\n")
            else:
                lf.write(f"[GEP2] {query} (taxid {found[query]}): target {target['level_rank']} "
                         f"{target['level_name']} (taxid {target['level_taxid']})\n")


def _kraken_summarise(report_path, target):
    """Read groups of a Kraken2 report, adding up to 100%: unclassified, target (clade of
    the target order), unspecific (above it on its lineage, or root), other eukaryotes,
    human (if outside the target), bacteria, archaea, viruses, other; plus the top
    off-target species."""
    nodes, path = {}, []
    for line in open(report_path):
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 6:
            continue
        clade, direct, code, taxid = int(cols[1]), int(cols[2]), cols[3], int(cols[4])
        name = cols[5]
        depth = (len(name) - len(name.lstrip(" "))) // 2
        if taxid == 0:
            nodes[0] = {"clade": clade, "direct": direct, "code": code, "name": "unclassified", "path": []}
            continue
        path = path[:depth] + [taxid]
        nodes[taxid] = {"clade": clade, "direct": direct, "code": code, "name": name.strip(), "path": list(path)}

    clade = lambda t: nodes[t]["clade"] if t in nodes else 0
    direct = lambda t: nodes[t]["direct"] if t in nodes else 0
    unclassified = clade(0)
    total = unclassified + clade(1)

    level = target.get("level_taxid") if target else None
    lineage = target.get("lineage", []) if target else []
    ancestors = lineage[lineage.index(level) + 1:] if level in lineage else []

    euk, bac, arc, vir = clade(2759), clade(2), clade(2157), clade(10239)
    human_in_target = level is not None and 9606 in nodes and level in nodes[9606]["path"]
    human = 0 if human_in_target else clade(9606)
    target_n = clade(level) if level is not None else None
    # reads stopped on the target's own lineage inside Eukaryota: compatible, not informative
    lineage_direct = sum(direct(t) for t in ancestors if t not in (1, 131567) and t in nodes
                         and 2759 in nodes[t]["path"])
    unspecific = direct(1) + direct(131567) + lineage_direct
    other_euk = euk - (target_n or 0) - human - lineage_direct
    other = clade(1) - euk - bac - arc - vir - direct(1) - direct(131567)

    def pct(n):
        return round(100 * n / total, 2) if total and n is not None else None

    off_target = [(n["clade"], n["name"]) for t, n in nodes.items()
                  if n["code"] == "S" and not (level is not None and level in n["path"])
                  and not (human_in_target and t == 9606)]
    top = max(off_target) if off_target else None

    return {
        "total": total,
        "unclassified_pct": pct(unclassified),
        "target_pct": pct(target_n),
        "unspecific_pct": pct(unspecific),
        "other_eukaryota_pct": pct(other_euk),
        "human_pct": None if human_in_target else pct(human),
        "bacteria_pct": pct(bac),
        "archaea_pct": pct(arc),
        "viruses_pct": pct(vir),
        "other_pct": pct(other),
        "top_off_target": top[1] if top and top[0] > 0 else "none",
        "top_off_target_pct": pct(top[0]) if top and top[0] > 0 else 0.0,
        "_informative_euk": (target_n or 0) + other_euk + human,
        "_target_n": target_n,
    }


rule G04_kraken_summary:
    """Per read set: reads per group as MultiQC custom content (General Statistics columns
    + a table) and warnings in the log. Thresholds are provisional."""
    input:
        report = os.path.join(config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
                              "reads", "{read_type}", "kraken2", "{group}.kraken2.report"),
        target = os.path.join(config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
                              "reads", "kraken2_target.json")
    output:
        gstats = os.path.join(config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
                              "reads", "{read_type}", "kraken2", "{group}_kraken2_mqc.yaml"),
        table = os.path.join(config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
                             "reads", "{read_type}", "kraken2", "{group}_kraken2_table_mqc.yaml")
    wildcard_constraints:
        read_type = "|".join(KRAKEN_READ_TYPES),
        group = r"[^/]+_Path\d+_[^/]+"
    threads: 1
    resources:
        mem_mb = 1000,
        runtime = 15
    log:
        os.path.join(config["OUT_FOLDER"], "GEP2_results", "data", "{species}",
                     "reads", "{read_type}", "logs", "G04_kraken_summary.{group}.log")
    run:
        import json
        _mk_dirs(log[0])
        target = json.load(open(input.target))
        s = _kraken_summarise(input.report, target)
        level = target.get("level_name") or "target"
        rank = target.get("level_rank") or "order"
        group = wildcards.group

        warnings = []
        informative = s["_informative_euk"]
        if s["_target_n"] is not None and s["total"] and informative >= 0.01 * s["total"] \
                and s["_target_n"] < 0.5 * informative:
            warnings.append(f"most eukaryotic reads placed below {rank} level fall outside {level} "
                            f"(top off-target: {s['top_off_target']}): wrong sample or species?")
        for key, label, limit in (("bacteria_pct", "bacteria", 5), ("other_eukaryota_pct", "other eukaryotes", 5),
                                  ("human_pct", "human", 1), ("viruses_pct", "viruses", 5)):
            if s[key] is not None and s[key] >= limit:
                warnings.append(f"{s[key]:.2f}% of reads classified as {label}")
        status = ("WARNING: " + "; ".join(warnings) if warnings else "ok").replace('"', "'")
        s["top_off_target"] = s["top_off_target"].replace('"', "'")

        def num(v):
            return "null" if v is None else v

        with open(output.gstats, "w") as fh:
            fh.write(
                'id: "gep2_kraken2_screen"\n'
                'section_name: "Kraken2 read screen"\n'
                f'description: "Kraken2 core_nt on {KRAKEN_SUBSET:,} reads (pairs) per read set"\n'
                'plot_type: "generalstats"\n'
                'headers:\n'
                '  - k2_target_pct:\n'
                f'      title: "{level} %"\n'
                f'      description: "Reads classified in the target {rank} ({level})"\n'
                '      suffix: "%"\n      min: 0\n      max: 100\n      scale: "RdYlGn"\n'
                '  - k2_other_euk_pct:\n'
                '      title: "Other euk. %"\n'
                f'      description: "Reads classified in eukaryotes outside {level} (human apart)"\n'
                '      suffix: "%"\n      min: 0\n      max: 100\n      scale: "RdYlGn-rev"\n'
                '  - k2_bacteria_pct:\n'
                '      title: "Bacteria %"\n      suffix: "%"\n      min: 0\n      max: 100\n      scale: "RdYlGn-rev"\n'
                '  - k2_human_pct:\n'
                '      title: "Human %"\n'
                '      description: "Reads classified as human (not shown when the target order holds humans)"\n'
                '      suffix: "%"\n      min: 0\n      max: 100\n      scale: "RdYlGn-rev"\n'
                'data:\n'
                f'  {group}:\n'
                f'    k2_target_pct: {num(s["target_pct"])}\n'
                f'    k2_other_euk_pct: {num(s["other_eukaryota_pct"])}\n'
                f'    k2_bacteria_pct: {num(s["bacteria_pct"])}\n'
                f'    k2_human_pct: {num(s["human_pct"])}\n'
            )
        cols = [("unclassified_pct", "Unclassified %"), ("target_pct", f"{level} %"),
                ("unspecific_pct", "Unspecific %"), ("other_eukaryota_pct", "Other euk. %"),
                ("human_pct", "Human %"), ("bacteria_pct", "Bacteria %"), ("archaea_pct", "Archaea %"),
                ("viruses_pct", "Viruses %"), ("other_pct", "Other %")]
        with open(output.table, "w") as fh:
            fh.write(
                'id: "gep2_kraken2_screen_table"\n'
                'section_name: "Kraken2 read screen: read groups"\n'
                f'description: "Share of the subset in each group (adds up to 100%). Unspecific: assigned only '
                f'above {rank} level on the target lineage, or to root. Thresholds for warnings are provisional."\n'
                'plot_type: "table"\n'
                'pconfig:\n  id: "gep2_kraken2_screen_table"\n  title: "Kraken2 read groups"\n'
                'headers:\n'
                + "".join(f'  {k}:\n    title: "{t}"\n    suffix: "%"\n    min: 0\n    max: 100\n' for k, t in cols)
                + '  top_off_target:\n    title: "Top off-target species"\n'
                + '  top_off_target_pct:\n    title: "Top off-target %"\n    suffix: "%"\n'
                + '  status:\n    title: "Status"\n'
                + 'data:\n'
                + f'  {group}:\n'
                + "".join(f'    {k}: {num(s[k])}\n' for k, _ in cols)
                + f'    top_off_target: "{s["top_off_target"]}"\n'
                + f'    top_off_target_pct: {s["top_off_target_pct"]}\n'
                + f'    status: "{status}"\n'
            )
        with open(log[0], "w") as lf:
            lf.write(f"[GEP2] {wildcards.species}/{group}: {s['total']} reads; "
                     + ", ".join(f"{t} {s[k]}" for k, t in cols) + "\n")
            lf.write(f"[GEP2] {status}\n")
        if warnings:
            print(f"[GEP2] Kraken2 {wildcards.species}/{group}: {status}")
