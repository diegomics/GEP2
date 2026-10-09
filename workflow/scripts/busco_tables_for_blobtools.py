#!/usr/bin/env python3
"""Write the gene tables of a GEP2 compleasm/BUSCO *_results.tar.gz in the form
`blobtools add --busco` (BlobToolKit >= 4.5.4) reads, one file per lineage.

  BUSCO:     ./run_<lineage>/full_table.tsv               -> written unchanged
  compleasm: ./<lineage>/full_table_busco_format.tsv      -> two changes:
     - adds the two header lines blobtools needs (compleasm does not write them)
     - keeps one row per Retrocopy gene (blobtools would count it as duplicated)
  Tables nested deeper (e.g. BUSCO's auto_lineage/ trials) are ignored.

The paths written are printed (stdout) with the most specific lineage first and
eukaryota last: the first table added to a BlobDir is the one in the snail plot.
A one-line summary of what was done to each table goes to stderr (the job log).

Usage: busco_tables_for_blobtools.py ARCHIVE OUTDIR COMPLEASM_VERSION
"""
import os
import re
import sys
import tarfile

LINEAGE = r"[A-Za-z0-9_]+_odb\d+(?:\.\d+)?"
BUSCO_TABLE = re.compile(rf"^(?:\./)?run_({LINEAGE})/full_table\.tsv$")
COMPLEASM_TABLE = re.compile(rf"^(?:\./)?({LINEAGE})/full_table_busco_format\.tsv$")


def fix_compleasm(text, lineage, version):
    """Add the BUSCO-style header lines and keep one row per Retrocopy gene.
    Returns the new text and a short note of what was changed."""
    if "lineage dataset is:" in text:          # header already present: leave it
        header = []
        note = "header already present"
    else:
        genes = {line.split("\t")[0] for line in text.splitlines()
                 if line.strip() and not line.startswith("#")}
        header = [f"# compleasm version is: {version}",
                  f"# The lineage dataset is: {lineage} (Creation date: NA, "
                  f"number of genomes: NA, number of BUSCOs: {len(genes)})"]
        note = f"added BUSCO-style header (compleasm {version}, {len(genes)} genes)"
    rows, retrocopies, retro_rows = [], set(), 0
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) > 1 and parts[1] == "Retrocopy":
            retro_rows += 1
            if parts[0] in retrocopies:
                continue
            retrocopies.add(parts[0])
        rows.append(line)
    if retro_rows:
        note += f"; Retrocopy rows {retro_rows} -> {len(retrocopies)} (one per gene)"
    else:
        note += "; no Retrocopy rows"
    return "\n".join(header + rows) + "\n", note


def main(archive, outdir, version):
    os.makedirs(outdir, exist_ok=True)
    written = []
    with tarfile.open(archive, "r|gz") as tar:      # reads the archive once
        for member in tar:
            busco = BUSCO_TABLE.match(member.name)
            compleasm = COMPLEASM_TABLE.match(member.name)
            if not member.isfile() or not (busco or compleasm):
                continue
            lineage = (busco or compleasm).group(1)
            text = tar.extractfile(member).read().decode("utf-8", "replace")
            if compleasm:
                text, note = fix_compleasm(text, lineage, version)
                print(f"[GEP2] compleasm {lineage}: {note}", file=sys.stderr)
            else:
                print(f"[GEP2] BUSCO {lineage}: table used as is", file=sys.stderr)
            path = os.path.join(outdir, f"{lineage}_full_table.tsv")
            with open(path, "w") as fh:
                fh.write(text)
            written.append((lineage.startswith("eukaryota"), lineage, path))
    for _, _, path in sorted(written):
        print(path)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
