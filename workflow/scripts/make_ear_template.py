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
#   - merqury_folder (Merqury/MerquryFK output folder; per haplotype, same value)
#   - fcs_gx         (FCS-GX report, discovered inside the fcs dir; per haplotype)
#   - blobplot       (blobtools *.blob.circle.png, discovered inside the blob dir; per haplotype)

import argparse
import glob
import os
import re
import sys

__version__ = '0.0.1'


# -------------------------------------------------------------------------------
# EAR TEMPLATE. Anything not auto-filled stays as-is for manual editing.
# -------------------------------------------------------------------------------
EAR_TEMPLATE = '''# SAMPLE INFORMATION
tolid: <Insert ToLID> # (e.g., "mEleMax1" or empty "")
species: <Insert species name> # MANDATORY! (e.g., "Homo sapiens")
sex: <Insert sample sex> # (e.g., "XX", "X0", "ZW", "unknown", "NA"... or empty "")
submitter: <Insert submitter full name> # (e.g., "John Doe" or empty "")
affiliation: <Insert submitter affiliation> # (e.g., "Sanger" or empty "")
tags: <Insert tag> # (e.g., "ERGA-BGE", "ERGA-Community"... or empty "")
# SEQUENCING DATA
DATA: # MANDATORY! Add below name of available data, one below the other
  - <Insert data type>: <insert data coverage> # (e.g., HiFi: 30). Coverage can be "NA" or empty (e.g., Hi-C: "")
# GENOME PROFILING DATA
PROFILING:
  genomescope: <Insert /path/to/genomescope_linear_plot.png> # ...or empty ""
  smudgeplot: <Insert /path/to/smudgeplot.png> # ...or empty ""
# ASSEMBLY DATA
ASSEMBLIES:
  hap1: # repeat this block if more than one haplotype (up to 4), call it hap2...
    gfastats_nstar: <Insert /path/to/gfastats--nstar-report.txt> # MANDATORY!
    busco_short: <Insert /path/to/busco_short_summary.txt> # MANDATORY!
    merqury_folder: <Insert /path/to/Merqury_results_foder> # ...or empty ""
    fcs_gx: <Insert /path/to/fcs-gx_action_report.txt> # ...or empty ""
    blobplot: <Insert /path/to/blobplot.png> # ...or empty ""
# METHODS DATA
PIPELINE: <Insert methods notes> # text block in quotes "", can also be path to txt file or empty ""
# CURATION NOTES
NOTES:
  obs_ploidy: <Insert observed ploidy number> # integer
  obs_haploid_num: <Insert observed haploid number> # MANDATORY! integer
  obs_sex: <Insert observed sex> # (e.g., "XX", "X0", "ZW", "unknown", "NA"... or empty "")
  curation_notes: <Insert curation notes> # text block in quotes "" related to the curation, plastids, sample... can also be path to txt file or empty ""
  alignment_png: <Insert /path/to/alignment.png> # like d-genies, syri, etc... or empty ""
'''

MAX_HAPS = 4


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
    hap_fields : the field lines of the hap block (up to the METHODS DATA comment)
    postamble  : from the METHODS DATA comment to the end
    """
    lines = template.splitlines()

    try:
        hap_idx = next(i for i, ln in enumerate(lines) if re.match(r'^\s*hap1:', ln))
        post_idx = next(i for i, ln in enumerate(lines) if ln.strip() == '# METHODS DATA')
    except StopIteration:
        raise ValueError("EAR_TEMPLATE is missing its 'hap1:' or '# METHODS DATA' anchor")

    preamble = lines[:hap_idx]
    hap_header = lines[hap_idx]
    hap_fields = lines[hap_idx + 1:post_idx]
    postamble = lines[post_idx:]
    return preamble, hap_header, hap_fields, postamble


def build_template(species, genomescope, gfastats, busco_short,
                   merqury_dir, fcs_dirs, blob_dirs):
    """Assemble the filled EAR text. 'species' is already cleaned (spaces)."""
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

    # Preamble: fills species + genomescope; every other line stays
    out.extend(_fill_line(ln, shared_values) for ln in preamble)

    # One block per haplotype
    for i in range(num_haps):
        # hap1 keeps the template's header (comment included); extras are bare keys
        out.append(hap_header if i == 0 else f'  hap{i + 1}:')

        hap_values = {
            'gfastats_nstar': _abspath_or_none(gfastats[i]),
            'busco_short': _abspath_or_none(busco_short[i]) if i < len(busco_short) else None,
            'merqury_folder': merqury_folder,
            'fcs_gx': _find_one(fcs_dirs[i], '*fcs_gx_report.txt') if i < len(fcs_dirs) else None,
            'blobplot': _find_one(blob_dirs[i], '*.blob.circle.png') if i < len(blob_dirs) else None,
        }
        out.extend(_fill_line(ln, hap_values) for ln in hap_fields)

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
    parser.add_argument('-o', '--output', required=True,
                        help='Output EAR YAML file')

    args = parser.parse_args()

    try:
        species_clean = args.species.replace('_', ' ')
        num_haps = len(args.gfastats)

        print(f"Building EAR template for: {species_clean} / {args.assembly}")
        print(f"Haplotypes (gfastats files): {num_haps}")

        # Per-hap lists should be empty (analysis off) or length == num_haps.
        # Anything else is matched by position, mirroring the report's behaviour.
        for name, files in (('busco-short', args.busco_short),
                            ('fcs-dir', args.fcs_dir),
                            ('blob-dir', args.blob_dir)):
            if files and len(files) != num_haps:
                print(f"Warning: {len(files)} --{name} given but {num_haps} haplotype(s); "
                      f"matched by position (extras ignored, missing -> placeholder)")

        content = build_template(
            species_clean,
            args.genomescope,
            args.gfastats,
            args.busco_short,
            args.merqury_dir,
            args.fcs_dir,
            args.blob_dir,
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