#!/usr/bin/env python3

# GEP2 ERGA Assembly Report (EAR) YAML Template Generator
# by Diego De Panis, 2026
# This script is part of the GEP2 pipeline
# note: AI tools may have been used to improve, clean and/or comment this version of the code

# Fills an EAR YAML template with the analysis outputs that are available:
#   - species        (from --species, underscores -> spaces)
#   - genomescope    (GenomeScope2 linear plot; shared, under PROFILING)
#   - gfastats_nstar (gfastats stats file; per haplotype)
#   - busco_short    (BUSCO or compleasm summary.txt; per haplotype)
#   - merqury_folder (Merqury/MerquryFK output folder; one run covers all haps, so hap1 only)
#   - fcs_gx         (FCS-GX action report, discovered inside the fcs dir; per haplotype)
#   - blobplot       (blobtools *.blob.circle.png, discovered inside the blob dir; per haplotype)
#   - DATA           (read amount in Gbp per technology, summed from GEP2's per-file
#                     seqtk size counts of this assembly's reads: processed reads when
#                     processing was on for that technology, raw reads otherwise)
# Not auto-filled (kept as placeholders): smudgeplot, pretext_temp_link, PIPELINE and NOTES.

import argparse
import glob
import os
import re
import sys

__version__ = '0.0.3'

# make_EAR.py template version mirrored below (ERGA-consortium/EARs scripts/EAR_basic_template.yaml)
EAR_TEMPLATE_VERSION = '26.11.16'


# -------------------------------------------------------------------------------
# EAR TEMPLATE. Anything not auto-filled stays as-is for manual editing
# (make_EAR.py shows any leftover <Insert ...> placeholder as NA).
# -------------------------------------------------------------------------------
EAR_TEMPLATE = '''# ERGA Assembly Report (EAR) - YAML template for make_EAR.py (v26.11.16)
# Paths can be absolute or relative. Relative paths are read from the folder you run
# make_EAR.py in. When in doubt, use full paths.
# Anything left empty (""), written as NA (also n/a, none, "-" in quotes), or left with its
# <Insert ...> placeholder shows as NA in the report.
# Plots can be PNG or SVG.

# SAMPLE INFORMATION
tolid: <Insert ToLID> # (e.g., "mEleMax1" or empty "")
species: <Insert species name> # MANDATORY! (e.g., "Homo sapiens")
sex: <Insert sample sex> # (e.g., "XX", "X0", "ZW", "unknown", "NA"... or empty "")
submitter: <Insert submitter full name> # (e.g., "John Doe" or empty "")
affiliation: <Insert submitter affiliation> # (e.g., "Sanger" or empty "")
tags: <Insert tag> # (e.g., "ERGA-BGE", "ERGA-Community"... or empty "")

# SEQUENCING DATA
DATA: # one line per data type with its amount in Gbp; the type is required, the amount can be "NA" or empty
  - <Insert data type>: <Insert read amount in Gbp> # (e.g., HiFi: 30.23 or Hi-C: "")

# GENOME PROFILING DATA
PROFILING:
  genomescope: <Insert /path/to/genomescope_linear_plot.png> # ...or empty ""
  smudgeplot: <Insert /path/to/smudgeplot.png> # ...or empty ""

# ASSEMBLY DATA
ASSEMBLIES:
  hap1: # one block per assembly, named hap1, hap2, hap3, hap4 (e.g. a primary/alternate pair is hap1 and hap2)
    gfastats_nstar: <Insert /path/to/gfastats--nstar-report.txt> # MANDATORY!
    busco_short: <Insert /path/to/busco_short_summary.txt> # MANDATORY! BUSCO short summary or compleasm summary.txt
    merqury_folder: <Insert /path/to/Merqury_results_folder> # Merqury or MerquryFK output folder of one run covering all haplotypes; giving it in hap1 is enough ...or empty ""
    fcs_gx: <Insert /path/to/fcs-gx_action_report.txt> # FCS-GX action report (not the taxonomy report) ...or empty ""
    blobplot: <Insert /path/to/blobplot.png> # ...or empty ""
    pretext_temp_link: <Insert .pretext file web link> # can also be a web folder or empty ""

# METHODS DATA
PIPELINE: <Insert methods notes> # free text in quotes "" (e.g., "Hifiasm v0.25 Hi-C mode, YaHS default. All assembled in Galaxy"), or the path to a .txt file (keeps line breaks), or empty ""

# CURATION NOTES
NOTES:
  obs_ploidy: <Insert observed ploidy number> # integer; also gives the coverage per haplotype
  obs_haploid_num: <Insert observed haploid number> # MANDATORY! the haploid number n as a whole number (e.g., 23)
  obs_sex: <Insert observed sex> # (e.g., "XX", "X0", "ZW", "unknown", "NA"... or empty "")
  curation_notes: <Insert curation notes> # free text in quotes "" about the curation, plastids, sample (e.g., low input, pooled)..., or the path to a .txt file (keeps line breaks), or empty "". Shown as "Comments" in the report
  alignment_png: <Insert /path/to/alignment.png> # like d-genies, syri, etc... or empty ""
'''

MAX_HAPS = 4

# Keys that describe one run for all haplotypes: written in hap1 only, dropped from
# hap2+ (make_EAR.py reads them from whichever haplotype block has them).
RUN_LEVEL_KEYS = {'merqury_folder'}

# GEP2 read types -> DATA type names. make_EAR.py matches HiFi, ONT and Hi-C to the
# EBP minimum coverage per haplotype; Illumina and 10x are shown in Gbp only.
DATA_TYPE_NAMES = {'hifi': 'HiFi', 'ont': 'ONT', 'hic': 'Hi-C', 'illumina': 'Illumina', '10x': '10x'}
DATA_TYPE_ORDER = ['hifi', 'ont', 'hic', 'illumina', '10x']


# -------------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------------
def _abspath_or_none(path):
    """Absolute path if the target exists, else None (-> placeholder kept)."""
    if not path:
        return None
    if os.path.exists(path):
        return os.path.abspath(path)
    print(f"Warning: expected file not found, leaving placeholder: {path}")
    return None


def _find_one(directory, pattern):
    """First file matching 'pattern' inside 'directory' (abspath), else None.
    Used for FCS-GX and blobtools, whose filenames carry unpredictable bits
    (taxid / BlobDir prefix) and so are discovered at runtime rather than
    hard-wired as rule inputs."""
    if not directory or not os.path.isdir(directory):
        return None
    hits = sorted(glob.glob(os.path.join(directory, pattern)))
    if hits:
        return os.path.abspath(hits[0])
    print(f"Warning: no '{pattern}' found in {directory}, leaving placeholder")
    return None


def _fill_line(line, values):
    """Replace 'key: <Insert ...> # comment' with 'key: "value"' when we have a
    non-empty value for that key. Any line we don't manage - or manage but have no
    value for - is returned untouched, so its placeholder and inline comment are
    preserved for manual completion."""
    m = re.match(r'^(?P<indent>\s*)(?P<key>[A-Za-z_]\w*):\s*<[^>]*>.*$', line)
    if not m:
        return line
    key = m.group('key')
    if key not in values:
        return line
    val = values[key]
    if not val:
        return line
    return f'{m.group("indent")}{key}: "{val}"'


def _split_template(template):
    """Split the template into (preamble, hap_header, hap_fields, postamble).

    preamble   : everything up to (not including) the 'hap1:' line
    hap_header : the 'hap1:' line (kept as-is for hap1)
    hap_fields : the field lines of the hap block (up to the METHODS DATA comment,
                 without the blank separator lines, so they are not repeated per hap)
    postamble  : the blank separator line(s) + from the METHODS DATA comment to the end
    """
    lines = template.splitlines()

    try:
        hap_idx = next(i for i, ln in enumerate(lines) if re.match(r'^\s*hap1:', ln))
        post_idx = next(i for i, ln in enumerate(lines) if ln.strip() == '# METHODS DATA')
    except StopIteration:
        raise ValueError("EAR_TEMPLATE is missing its 'hap1:' or '# METHODS DATA' anchor")

    # Blank lines just before '# METHODS DATA' separate sections; keep them once
    while post_idx > hap_idx + 1 and not lines[post_idx - 1].strip():
        post_idx -= 1

    preamble = lines[:hap_idx]
    hap_header = lines[hap_idx]
    hap_fields = lines[hap_idx + 1:post_idx]
    postamble = lines[post_idx:]
    return preamble, hap_header, hap_fields, postamble


def _line_key(line):
    """YAML key of a 'key: ...' line, else None."""
    m = re.match(r'^\s*([A-Za-z_]\w*):', line)
    return m.group(1) if m else None


def read_data_totals(specs):
    """Add up GEP2 per-file base counts per read type.

    specs: ['<read_type>=<counts.tsv>', ...]; each counts file is the B09 table
    (header 'file stage reads bases' + one row). The same counts file given twice
    is counted once. Returns {read_type: {'bases', 'files', 'stages'}} in DATA order.
    """
    totals, seen = {}, set()
    for spec in specs:
        if '=' not in spec:
            raise ValueError(f"--read-counts expects <read_type>=<file>, got: {spec}")
        read_type, path = spec.split('=', 1)
        read_type = read_type.strip().lower()
        path = os.path.abspath(path)
        if path in seen:
            continue
        seen.add(path)

        with open(path) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            values = fh.readline().rstrip('\n').split('\t')
        rec = dict(zip(header, values))
        if 'bases' not in rec:
            raise ValueError(f"No 'bases' column in read counts file: {path}")

        entry = totals.setdefault(read_type, {'bases': 0, 'files': 0, 'stages': set()})
        entry['bases'] += int(rec['bases'])
        entry['files'] += 1
        entry['stages'].add(rec.get('stage', 'raw'))

    order = {rt: i for i, rt in enumerate(DATA_TYPE_ORDER)}
    return dict(sorted(totals.items(), key=lambda kv: (order.get(kv[0], len(order)), kv[0])))


def _format_gbp(bases):
    """Gbp for the DATA value: 3 decimals (1 Mbp); below 1 Mbp, plain decimals down to
    single bases (never 1e-05, which YAML would read as text)."""
    gbp = bases / 1e9
    return f"{gbp:.3f}" if gbp >= 0.001 else f"{gbp:.9f}".rstrip("0")


def _data_lines(template_line, totals):
    """Replace the DATA placeholder line by one line per read type."""
    indent = template_line[:len(template_line) - len(template_line.lstrip())]
    lines = []
    for read_type, entry in totals.items():
        name = DATA_TYPE_NAMES.get(read_type, read_type)
        stage = '/'.join(sorted(entry['stages']))
        files = f"{entry['files']} file{'s' if entry['files'] != 1 else ''}"
        if entry['bases'] > 0:
            lines.append(f'{indent}- {name}: {_format_gbp(entry["bases"])}'
                         f' # GEP2: {stage} reads, {files}')
        else:
            print(f"Warning: 0 bases counted for {name}; leaving its amount as a placeholder")
            lines.append(f'{indent}- {name}: <Insert read amount in Gbp>'
                         f' # GEP2 counted 0 bases ({stage} reads, {files})')
    return lines


def build_template(species, genomescope, gfastats, busco_short,
                   merqury_dir, fcs_dirs, blob_dirs, data_totals=None):
    """Assemble the filled EAR text. 'species' is already cleaned (spaces).
    data_totals: output of read_data_totals(); empty/None keeps the DATA placeholder."""
    num_haps = len(gfastats)
    if num_haps < 1:
        raise ValueError("At least one gfastats file is required (one per haplotype)")
    if num_haps > MAX_HAPS:
        raise ValueError(f"EAR supports up to {MAX_HAPS} haplotypes, got {num_haps}")

    preamble, hap_header, hap_fields, postamble = _split_template(EAR_TEMPLATE)

    # Shared (top-level / PROFILING) fields
    shared_values = {
        'species': species,                            # always present
        'genomescope': _abspath_or_none(genomescope),  # shared across haps
    }

    # Merqury folder is resolved once and reused identically for every haplotype
    merqury_folder = None
    if merqury_dir and os.path.isdir(merqury_dir):
        merqury_folder = os.path.abspath(merqury_dir)
    elif merqury_dir:
        print(f"Warning: Merqury folder not found, leaving placeholder: {merqury_dir}")

    out = []

    # Preamble: fills species + genomescope (+ DATA lines); every other line stays
    data_filled = False
    for ln in preamble:
        if data_totals and re.match(r'^\s*-\s*<Insert data type>', ln):
            out.extend(_data_lines(ln, data_totals))
            data_filled = True
        else:
            out.append(_fill_line(ln, shared_values))
    if data_totals and not data_filled:
        raise ValueError("EAR_TEMPLATE is missing its '- <Insert data type>: ...' DATA line")

    # One block per haplotype
    for i in range(num_haps):
        # hap1 keeps the template's header (comment included); extras are bare keys
        out.append(hap_header if i == 0 else f'  hap{i + 1}:')

        hap_values = {
            'gfastats_nstar': _abspath_or_none(gfastats[i]),
            'busco_short': _abspath_or_none(busco_short[i]) if i < len(busco_short) else None,
            'merqury_folder': merqury_folder,
            # FCS-GX *action* report ({prefix}.{taxid}.fcs_gx_report.txt), not *.taxonomy.rpt
            'fcs_gx': _find_one(fcs_dirs[i], '*fcs_gx_report.txt') if i < len(fcs_dirs) else None,
            'blobplot': _find_one(blob_dirs[i], '*.blob.circle.png') if i < len(blob_dirs) else None,
        }
        fields = hap_fields if i == 0 else [ln for ln in hap_fields
                                            if _line_key(ln) not in RUN_LEVEL_KEYS]
        out.extend(_fill_line(ln, hap_values) for ln in fields)

    # Postamble: PIPELINE + NOTES, kept
    out.extend(postamble)

    return '\n'.join(out) + '\n'


# -------------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description='Fill an ERGA EAR YAML template from available GEP2 outputs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single assembly (hap1 only)
  python make_EAR_template.py -s Rosalia_alpina -a hifiasm_l2_h_yahs \\
    -g stats.hap1.txt \\
    --busco-short hap1_summary.txt \\
    --genomescope hap1_linear_plot.png \\
    --merqury-dir merqury/ \\
    --fcs-dir decontamination/fcs-gx/hap1/ \\
    --blob-dir decontamination/blobtools/hap1/ \\
    -o hifiasm_l2_h_yahs_EAR.yaml

  # Two assemblies (hap1 + hap2; --merqury-dir is shared)
  python make_EAR_template.py -s Rosalia_alpina -a hifiasm_l2_h_yahs \\
    -g stats.hap1.txt stats.hap2.txt \\
    --busco-short hap1_summary.txt hap2_summary.txt \\
    --genomescope hap1_linear_plot.png \\
    --merqury-dir merqury/ \\
    --fcs-dir fcs/hap1/ fcs/hap2/ \\
    --blob-dir blob/hap1/ blob/hap2/ \\
    --read-counts hifi=counts/hifi_Path1_x_filtered.processed.bases.tsv \\
                  hic=counts/hic_Path1_y_1_trimmed.processed.bases.tsv \\
                  hic=counts/hic_Path1_y_2_trimmed.processed.bases.tsv \\
    -o hifiasm_l2_h_yahs_EAR.yaml
        """
    )

    parser.add_argument('-v', '--version', action='version',
                        version=f'%(prog)s {__version__}')
    parser.add_argument('-s', '--species', required=True,
                        help='Species name (with underscore, e.g., Rosalia_alpina)')
    parser.add_argument('-a', '--assembly', required=True,
                        help='Assembly ID (used for logging)')
    parser.add_argument('-g', '--gfastats', required=True, nargs='+',
                        help='gfastats stats file(s) - one per haplotype (defines hap count)')
    parser.add_argument('--busco-short', required=False, nargs='+', default=[],
                        help='BUSCO or compleasm summary.txt file(s) - one per haplotype')
    parser.add_argument('--genomescope', required=False,
                        help='GenomeScope2 linear plot PNG (shared)')
    parser.add_argument('--merqury-dir', required=False,
                        help='Merqury/MerquryFK output folder (shared across haplotypes)')
    parser.add_argument('--fcs-dir', required=False, nargs='+', default=[],
                        help='FCS-GX output dir(s) - one per haplotype (report globbed at runtime)')
    parser.add_argument('--blob-dir', required=False, nargs='+', default=[],
                        help='blobtools output dir(s) - one per haplotype (png globbed at runtime)')
    parser.add_argument('--read-counts', required=False, nargs='+', default=[],
                        help='<read_type>=<counts.tsv> per read file of this assembly (GEP2 '
                             'seqtk size tables); summed per read type into DATA (Gbp)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output EAR YAML file')

    args = parser.parse_args()

    try:
        species_clean = args.species.replace('_', ' ')
        num_haps = len(args.gfastats)

        print(f"Building EAR template (make_EAR.py template v{EAR_TEMPLATE_VERSION}) "
              f"for: {species_clean} / {args.assembly}")
        print(f"Haplotypes (gfastats files): {num_haps}")

        # Per-hap lists should be empty (analysis off) or length == num_haps.
        # Anything else is matched by position, mirroring the report's behaviour.
        for name, files in (('busco-short', args.busco_short),
                            ('fcs-dir', args.fcs_dir),
                            ('blob-dir', args.blob_dir)):
            if files and len(files) != num_haps:
                print(f"Warning: {len(files)} --{name} given but {num_haps} haplotype(s); "
                      f"matched by position (extras ignored, missing -> placeholder)")

        data_totals = read_data_totals(args.read_counts)
        for read_type, entry in data_totals.items():
            print(f"DATA {DATA_TYPE_NAMES.get(read_type, read_type)}: {entry['bases']} bases "
                  f"({_format_gbp(entry['bases'])} Gbp) from {entry['files']} file(s), "
                  f"{'/'.join(sorted(entry['stages']))} reads")

        content = build_template(
            species_clean,
            args.genomescope,
            args.gfastats,
            args.busco_short,
            args.merqury_dir,
            args.fcs_dir,
            args.blob_dir,
            data_totals,
        )

        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        with open(args.output, 'w') as fh:
            fh.write(content)

        print(f"EAR template written: {args.output}")
        print("Remember to complete the remaining <Insert ...> fields by hand.")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
