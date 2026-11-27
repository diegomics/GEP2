#!/usr/bin/env python3

# GEP2 Genome Stats Report Generator
# by Diego De Panis, 2025
# This script is part of the GEP2 pipeline
# note: AI tools may have been used to improve, clean and/or comment this version of the code

# Generates a markdown report from: gfastats, compleasm/BUSCO, Merqury/MerquryFK,
# GenomeScope2, Inspector, BlobToolKit, FCS-GX and the Hi-C stats.

import argparse
import math
import re
import struct
import sys
import os
import glob
import subprocess
import shutil
import requests
from urllib.parse import quote
from pathlib import Path

__version__ = '0.5'

# ---- Star-rating calibration -------------------------------------------------
# Adapts EBP's N50 targets to species with small chromosomes...
# The idea here is to cap absolute tiers at what chromosome can physically hold.
# expected_chr_size (E) = genome_size / haploid_number (mean chromosome).
# CHR_SCALE_FRAC (0.90) estimates chromosome N50 from the mean and keeps
# karyotype variance in the assembly's favor.

CHR_SCALE_FRAC = 0.90   # '****' bar: Chromosome-scale (EBP "C")
CHR_HALF_FRAC  = 0.50   # '***-' bar when the absolute tier is out of reach
CHR_TENTH_FRAC = 0.10   # '**--' bar likewise

# Absolute (uncapped) tiers, as (****, ***-, **--) in bp.
N50_BASE_TIERS = {
    'scaffold_n50': (100_000_000, 10_000_000, 100_000),
    'contig_n50':   ( 10_000_000,  1_000_000, 100_000),
}


def n50_tiers(metric_type, expected_chr_size=None):
    """(****, ***-, **--) bp thresholds for N50: the default tiers, each capped at a
    fraction of the expected chromosome size (a contig/scaffold can't be longer)."""
    t4, t3, t2 = N50_BASE_TIERS[metric_type]
    if not expected_chr_size:
        return t4, t3, t2
    return (min(t4, CHR_SCALE_FRAC * expected_chr_size),
            min(t3, CHR_HALF_FRAC  * expected_chr_size),
            min(t2, CHR_TENTH_FRAC * expected_chr_size))


def l90_tiers(haploid_number):
    """(****, ***-, **--) L90 thresholds for haploid number n: n, max(n+1, 1.5n), 5n."""
    if not haploid_number:
        return None
    return (haploid_number,
            max(haploid_number + 1, (haploid_number * 3) // 2),
            haploid_number * 5)


def _bp_short(value):
    """bp as 100 Mb / 1.72 Mb / 857.7 kb / 900 bp (legend and notes)."""
    for unit, size, digits in (("Gb", 1e9, 2), ("Mb", 1e6, 2), ("kb", 1e3, 1)):
        if value >= size:
            text = f"{value / size:.{digits}f}".rstrip('0').rstrip('.')
            return f"{text} {unit}"
    return f"{value:.0f} bp"


# This is a crap and isn't working yet, will work on it soon...
def convert_md_to_pdf(md_file, pdf_file=None):
    """Markdown -> PDF with pandoc + weasyprint."""
    if pdf_file is None:
        pdf_file = os.path.splitext(md_file)[0] + '.pdf'

    # Check for pandoc
    if not shutil.which('pandoc'):
        print("Error: pandoc not found. Install with: conda install pandoc", file=sys.stderr)
        return False

    # Check for weasyprint
    try:
        subprocess.run(['weasyprint', '--version'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: weasyprint not found. Install with: pip install weasyprint", file=sys.stderr)
        return False

    # CSS for better PDF styling
    css_content = """
    body {
        font-family: "DejaVu Sans", sans-serif;
        font-size: 11pt;
        line-height: 1.4;
        max-width: 100%;
        margin: 0 auto;
        padding: 20px;
    }
    h1, h2, h3, h4 { color: #2c3e50; margin-top: 1.5em; }
    h3 { border-bottom: 1px solid #ccc; padding-bottom: 0.3em; }
    table {
        border-collapse: collapse;
        width: 100%;
        margin: 1em 0;
        font-size: 10pt;
    }
    th, td {
        border: 1px solid #ddd;
        padding: 8px;
        text-align: left;
    }
    th { background-color: #f5f5f5; font-weight: bold; }
    tr:nth-child(even) { background-color: #fafafa; }
    img { max-width: 100%; height: auto; margin: 1em 0; }
    code {
        background-color: #f4f4f4;
        padding: 2px 6px;
        border-radius: 3px;
        font-family: "DejaVu Sans Mono", monospace;
    }
    a { color: #3498db; }
    hr { border: none; border-top: 1px solid #ccc; margin: 2em 0; }
    @page { margin: 1in; size: A4; }
    """

    # Write temporary CSS file
    css_file = md_file + '.tmp.css'
    try:
        with open(css_file, 'w') as f:
            f.write(css_content)

        # Run pandoc with weasyprint
        cmd = [
            'pandoc', md_file,
            '-o', pdf_file,
            '--pdf-engine=weasyprint',
            '--css=' + css_file,
            '--standalone',
            '--from', 'markdown-yaml_metadata_block',
        ]

        print(f"Converting to PDF: {pdf_file}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"Error converting to PDF: {result.stderr}", file=sys.stderr)
            return False

        print(f"PDF generated successfully: {pdf_file}")
        return True

    finally:
        # Clean up temp CSS
        if os.path.exists(css_file):
            os.remove(css_file)


def get_species_genomic_data_from_goat(species):
    """Get taxon ID, family/order lineage, haploid number and source from GoaT API based on species name.

    Looks up the species name on GoaT. If the species is not present, falls back to the genus.

    Returns a dict with keys: taxon_id, family, haploid_number, haploid_source, resolution_level, error.
    'family' holds the family-rank name when available, otherwise the order-rank name
    (or None if neither is present).
    'resolution_level' is 'species' or 'genus', indicating which name matched.
    """
    result = {
        'taxon_id': None,
        'family': None,
        'haploid_number': None,
        'haploid_source': None,
        'genome_size': None,
        'genome_size_source': None,
        'resolution_level': None,
        'error': None
    }

    try:
        # Try the species name; if GoaT has no match, fall back to the genus.
        parts = species.split()
        candidates = [(species, 'species')]
        if len(parts) >= 2:
            candidates.append((parts[0], 'genus'))

        first_result = None
        for name, level in candidates:
            name_encoded = quote(name)
            search_url = f'https://goat.genomehubs.org/api/v2/search?query=tax_name%28{name_encoded}%29&result=taxon'
            response = requests.get(search_url, timeout=30)
            response.raise_for_status()
            search_data = response.json()
            if search_data.get('results'):
                first_result = search_data['results'][0]['result']
                result['resolution_level'] = level
                break

        if first_result is None:
            result['error'] = "Species not found"
            return result

        taxon_id = first_result['taxon_id']
        result['taxon_id'] = taxon_id

        # Traverse the lineage array: prefer family, fall back to order
        if 'lineage' in first_result:
            family_name = None
            order_name = None
            for node in first_result['lineage']:
                rank = node.get('taxon_rank') or node.get('rank')
                if rank == 'family':
                    family_name = node.get('scientific_name')
                elif rank == 'order':
                    order_name = node.get('scientific_name')
            result['family'] = family_name or order_name

        record_url = f'https://goat.genomehubs.org/api/v2/record?recordId={taxon_id}&result=taxon&taxonomy=ncbi'
        record_response = requests.get(record_url, timeout=30)
        record_response.raise_for_status()
        record_data = record_response.json()

        if not record_data.get('records'):
            result['error'] = "No genomic records found"
            return result

        attributes = record_data['records'][0]['record'].get('attributes', {})
        haploid_info = attributes.get('haploid_number', {})
        haploid_number = haploid_info.get('value')
        haploid_source = haploid_info.get('aggregation_source')

        if haploid_number is not None:
            result['haploid_number'] = int(haploid_number)
            result['haploid_source'] = haploid_source
        else:
            result['error'] = "Haploid number not available"

        # Genome size, used only to estimate mean chromosome size (genome_size /
        # haploid_number) for the small-genome star-rating carve-out. Optional:
        # if absent, the carve-out simply never fires and scoring is unchanged.
        size_info = attributes.get('genome_size', {})
        if size_info.get('value') is not None:
            result['genome_size'] = size_info['value']
            result['genome_size_source'] = size_info.get('aggregation_source')

        return result

    except requests.exceptions.RequestException as e:
        result['error'] = f"API request failed: {e}"
        return result
    except (KeyError, ValueError, TypeError) as e:
        result['error'] = f"Data parsing error: {e}"
        return result
    except Exception as e:
        result['error'] = f"Unexpected error: {e}"
        return result


def parse_gfastats(filepath):
    """Parse gfastats output file and extract relevant metrics."""
    metrics = {}

    with open(filepath, 'r') as f:
        content = f.read()

    patterns = {
        'total_bp': r'Total scaffold length:\s*(\d+)',
        'gc_percent': r'GC content %:\s*([\d.]+)',
        'scaffolds': r'# scaffolds:\s*(\d+)',
        'scaffold_n50': r'Scaffold N50:\s*(\d+)',
        'scaffold_l50': r'Scaffold L50:\s*(\d+)',
        'scaffold_l90': r'Scaffold L90:\s*(\d+)',
        'contigs': r'# contigs:\s*(\d+)',
        'contig_n50': r'Contig N50:\s*(\d+)',
        'contig_l50': r'Contig L50:\s*(\d+)',
        'contig_l90': r'Contig L90:\s*(\d+)',
        'gaps': r'# gaps in scaffolds:\s*(\d+)',
        'total_gap_bp': r'Total gap length in scaffolds:\s*(\d+)'
    }

    # Track fields whose regex found nothing: they silently become "N/A" with a
    # '····' rating in the report, so a gfastats wording change would otherwise
    # drop metrics without anyone noticing. Reported in the run log instead.
    missing = []
    for key, pattern in patterns.items():
        match = re.search(pattern, content)
        if match:
            metrics[key] = match.group(1)
        else:
            missing.append(key)
    metrics['_missing'] = missing
    metrics['_source'] = filepath

    # Calculate gaps per Gbp
    if 'total_bp' in metrics and 'gaps' in metrics:
        total_bp = int(metrics['total_bp'])
        gaps = int(metrics['gaps'])
        gaps_per_gbp = (gaps / total_bp) * 1_000_000_000 if total_bp > 0 else 0
        metrics['gaps_per_gbp'] = gaps_per_gbp

    return metrics


def parse_compleasm(filepath):
    """Parse compleasm output file and extract completeness metrics.
    Duplicated = D + R: with --retrocopy compleasm splits genes with several complete
    copies (one with introns, the others without) out of D into R. Adding them back
    keeps the row equal to BUSCO's Duplicated (as make_EAR.py does)."""
    metrics = {
        'eukaryota_single': None,
        'eukaryota_dupl': None,
        'eukaryota_lineage': None,
        'other_lineage': None,
        'other_single': None,
        'other_dupl': None
    }

    with open(filepath, 'r') as f:
        content = f.read()

    sections = content.split('## lineage:')

    for section in sections[1:]:
        lines = section.strip().split('\n')
        lineage_name = lines[0].strip()

        s_match = re.search(r'S:([\d.]+)%', section)
        d_match = re.search(r'D:([\d.]+)%', section)
        r_match = re.search(r'R:([\d.]+)%', section)   # only with --retrocopy

        if s_match and d_match:
            s_percent = float(s_match.group(1))
            d_percent = float(d_match.group(1)) + (float(r_match.group(1)) if r_match else 0.0)

            if 'eukaryota' in lineage_name.lower():
                metrics['eukaryota_single'] = s_percent
                metrics['eukaryota_dupl'] = d_percent
                # Format nicely by targeting the _odb tag only (see parse_busco_summary),
                # so multi-underscore names like stramenopiles_alveolata_odb10 stay intact.
                metrics['eukaryota_lineage'] = lineage_name.replace('_odb', ' (odb') + ')' if '_odb' in lineage_name else lineage_name
            else:
                # Target the _odb tag only, matching parse_busco_summary. The previous
                # replace('_', ' (') broke multi-underscore lineages (double parens).
                other_name = lineage_name.replace('_odb', ' (odb') + ')' if '_odb' in lineage_name else lineage_name
                metrics['other_lineage'] = other_name
                metrics['other_single'] = s_percent
                metrics['other_dupl'] = d_percent

    return metrics


def parse_compleasm_full(filepath):
    """
    Parse compleasm full_table.tsv file and extract completeness metrics + frameshift rate.
    Lineage is derived from the parent directory name (e.g., rosales_odb12/).
    Only one lineage per file.
    """
    metrics = {
        'eukaryota_single': None,
        'eukaryota_dupl': None,
        'eukaryota_lineage': None,
        'other_lineage': None,
        'other_single': None,
        'other_dupl': None,
        'frameshift_rate': None,
        'frameshift_is_eukaryota': None
    }

    # Get lineage from directory name (e.g., rosales_odb12)
    lineage = os.path.basename(os.path.dirname(os.path.abspath(filepath)))

    seen = {}
    counts = {'Single': 0, 'Duplicated': 0, 'Retrocopy': 0,
              'Fragmented': 0, 'Interspaced': 0, 'Missing': 0}
    total_complete_copies = 0  # All Single/Duplicated rows (all physical copies)
    frameshifted_copies = 0    # Single/Duplicated rows with column 11 > 0

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue

            gene_id = parts[0]
            status = parts[1]

            # Count unique genes per category (for S/D percentages)
            if gene_id not in seen:
                seen[gene_id] = True
                if status in counts:
                    counts[status] += 1

            # Count all physical copies for frameshift rate
            if status in ('Single', 'Duplicated'):
                total_complete_copies += 1
                if len(parts) >= 11:
                    try:
                        if int(parts[10]) > 0:
                            frameshifted_copies += 1
                    except (ValueError, IndexError):
                        pass

    total = sum(counts.values())
    if total == 0:
        return metrics

    s_percent = (counts['Single'] / total) * 100
    # Duplicated = D + R (see parse_compleasm): a 'Retrocopy' gene has several complete
    # copies too, compleasm --retrocopy only labels it apart
    d_percent = ((counts['Duplicated'] + counts['Retrocopy']) / total) * 100

    # Frameshift rate
    if total_complete_copies > 0:
        frameshift_rate = (frameshifted_copies / total_complete_copies) * 100
    else:
        frameshift_rate = 0.0

    is_eukaryota = 'eukaryota' in lineage.lower()

    if is_eukaryota:
        metrics['eukaryota_single'] = s_percent
        metrics['eukaryota_dupl'] = d_percent
        metrics['eukaryota_lineage'] = lineage.replace('_odb', ' (odb') + ')' if '_odb' in lineage else lineage
    else:
        lineage_name = lineage.replace('_odb', ' (odb') + ')'
        metrics['other_lineage'] = lineage_name
        metrics['other_single'] = s_percent
        metrics['other_dupl'] = d_percent

    metrics['frameshift_rate'] = frameshift_rate
    metrics['frameshift_is_eukaryota'] = is_eukaryota

    # The report only ever shows the %, these are the numbers are for sanity-check
    metrics['_lineage_dir'] = lineage
    metrics['_counts'] = dict(counts)
    metrics['_n_markers'] = total
    metrics['_complete_copies'] = total_complete_copies
    metrics['_frameshifted_copies'] = frameshifted_copies

    return metrics


def parse_busco_summary(filepath):
    """
    Parse BUSCO short_summary.txt file.
    BUSCO only outputs one lineage per file, and we omit frameshift rates.
    """

    # Initialize dictionary exactly like parse_compleasm (no frameshift keys)
    metrics = {
        'eukaryota_single': None,
        'eukaryota_dupl': None,
        'eukaryota_lineage': None,
        'other_lineage': None,
        'other_single': None,
        'other_dupl': None
    }

    with open(filepath, 'r') as f:
        content = f.read()

    # Extract the single lineage name
    lineage_match = re.search(r'# The lineage dataset is:\s*([^\s(]+)', content)
    if not lineage_match:
        return metrics # Return empty metrics if file is malformed

    lineage_name = lineage_match.group(1)

    # Extract Single and Duplicated percentages
    s_match = re.search(r'S:([\d.]+)%', content)
    d_match = re.search(r'D:([\d.]+)%', content)

    if s_match and d_match:
        s_percent = float(s_match.group(1))
        d_percent = float(d_match.group(1))

        # Assign to eukaryota or 'other' based on the single lineage
        if 'eukaryota' in lineage_name.lower():
            metrics['eukaryota_single'] = s_percent
            metrics['eukaryota_dupl'] = d_percent
            metrics['eukaryota_lineage'] = lineage_name.replace('_odb', ' (odb') + ')' if '_odb' in lineage_name else lineage_name
        else:
            # Format nicely: e.g., "purpureocillium_takamizusanense (odb12)"
            other_name = lineage_name.replace('_odb', ' (odb') + ')' if '_odb' in lineage_name else lineage_name
            metrics['other_lineage'] = other_name
            metrics['other_single'] = s_percent
            metrics['other_dupl'] = d_percent

    return metrics


def detect_kmer_tool(qv_file, completeness_file):
    """Return the QV/completeness table label by peeking at the file header.
    MerquryFK writes a header row starting with 'Assembly'; Merqury has none.
    Defaults to 'MERQ' (Merqury) when neither file is readable."""
    for path in (qv_file, completeness_file):
        if not path:
            continue
        try:
            with open(path) as f:
                for line in f:
                    if line.strip():
                        first = line.split("\t", 1)[0].strip().lower()
                        return "MERQ.FK" if first == "assembly" else "MERQ"
        except OSError:
            continue
    return "MERQ"


def _order_merqury_rows(rows, num_assemblies):
    """Merqury rows [(label, value)] -> values in assembly order (asm1, asm2, ...);
    other labels keep their file order."""
    by_label = {label.strip().lower(): value for label, value in rows}
    wanted = [f"asm{i + 1}" for i in range(num_assemblies)]
    if all(w in by_label for w in wanted):
        return [by_label[w] for w in wanted]
    values = [value for _, value in rows][:num_assemblies]
    return values + [None] * (num_assemblies - len(values))


def parse_merqury_qv(filepath, num_assemblies):
    """Parse a Merqury or MerquryFK .qv summary file.

    Merqury (no header):   <asm>  <asm-only>  <total>  <QV>  <error rate>   -> QV at col 3
    MerquryFK (header):    Assembly  No Support  Total  Error %  QV         -> QV by name
    Returns QV per assembly, in order; ignores the trailing 'both' row.
    """
    qv_values = []

    try:
        with open(filepath) as f:
            rows = [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]
        if not rows:
            return [None] * num_assemblies

        qv_idx, start = 3, 0
        if rows[0] and rows[0][0].strip().lower() == "assembly":   # MerquryFK header
            start = 1
            qv_idx = next((i for i, c in enumerate(rows[0])
                           if c.strip().upper() == "QV"), len(rows[0]) - 1)

        labelled = []
        for parts in rows[start:]:
            label = parts[0].strip().lower() if parts else ""
            if label == "both" or "+" in label:      # skip combined row
                continue
            if len(parts) > qv_idx:
                try:
                    labelled.append((label, float(parts[qv_idx])))
                except ValueError:
                    pass

        qv_values = _order_merqury_rows(labelled, num_assemblies)

    except Exception as e:
        print(f"Warning: Could not parse QV file: {e}")
        qv_values = [None] * num_assemblies

    return qv_values


def parse_merqury_completeness(filepath, num_assemblies):
    """Parse a Merqury or MerquryFK completeness file.

    Merqury (no header):   <asm>  all  <found>  <total>  <completeness %>   -> last col
    MerquryFK (header):    Assembly  Region  Found  Total  % Covered        -> by name
    Both put the percentage last; ignores the trailing combined row.
    """
    completeness_values = []

    try:
        with open(filepath) as f:
            rows = [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]
        if not rows:
            return [None] * num_assemblies

        pct_idx, start = 4, 0
        if rows[0] and rows[0][0].strip().lower() == "assembly":   # MerquryFK header
            start = 1
            pct_idx = next((i for i, c in enumerate(rows[0])
                            if c.strip().lower() in ("% covered", "completeness")),
                           len(rows[0]) - 1)

        labelled = []
        for parts in rows[start:]:
            label = parts[0].strip().lower() if parts else ""
            if label == "both" or "+" in label:      # skip asm1+asm2 row
                continue
            if len(parts) > pct_idx:
                try:
                    labelled.append((label, float(parts[pct_idx])))
                except ValueError:
                    pass

        completeness_values = _order_merqury_rows(labelled, num_assemblies)

    except Exception as e:
        print(f"Warning: Could not parse completeness file: {e}")
        completeness_values = [None] * num_assemblies

    return completeness_values


def parse_merqury_both(filepath, column):
    """Return the combined asm1+asm2 value from a Merqury/MerquryFK .qv
    (column="qv") or .completeness.stats (column="completeness") file.
    The combined row is labelled 'both'/'Both' (or 'asm1+asm2' in MerquryFK
    completeness); returns None if the file has no such row."""
    try:
        with open(filepath) as f:
            rows = [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]
        idx = 3 if column == "qv" else 4                   # Merqury (no header)
        if rows and rows[0][0].strip().lower() == "assembly":   # MerquryFK header
            names = ("qv",) if column == "qv" else ("% covered", "completeness")
            idx = next((i for i, c in enumerate(rows[0])
                        if c.strip().lower() in names), len(rows[0]) - 1)
        for parts in rows:
            label = parts[0].strip().lower()
            if (label == "both" or "+" in label) and len(parts) > idx:
                return float(parts[idx])
    except (OSError, ValueError):
        pass
    return None


# ---- K-mer false duplications (as ERGA's make_EAR.py / Merqury's false_duplications.sh)
# From each assembly's copy-number spectrum: Merqury .spectra-cn.hist, or MerquryFK
# .spectra-cn.cni (only kept with -k)
CN_CLASSES = ["read-only", "1", "2", "3", "4", ">4"]


def read_cn_hist(path):
    """Copy-number spectrum of a Merqury .spectra-cn.hist as
    ({copy class: {read multiplicity: distinct k-mers}}, None); None = no upper
    limit (see read_cn_cni). Raises OSError / ValueError."""
    spectrum = {label: {} for label in CN_CLASSES}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) != 3 or parts[0] not in spectrum:
                continue                      # header ('Copies ...') or foreign line
            try:
                spectrum[parts[0]][int(parts[1])] = int(parts[2])
            except ValueError:
                continue
    if not spectrum["1"]:
        raise ValueError("no single-copy rows; is it a <prefix>.<asm>.spectra-cn.hist?")
    return spectrum, None


def read_cn_cni(path):
    """(spectrum, limit) from a MerquryFK .spectra-cn.cni: one FastK histogram per copy
    class (little-endian int32 k, low, high, 2x int64, then high-low+1 int64 counts).
    The last bin may hold everything above it, so it is dropped and 'high' returned as
    the limit. Raises OSError / ValueError."""
    with open(path, "rb") as fh:
        data = fh.read()
    histograms, offset = [], 0
    while offset < len(data):
        if len(data) - offset < 28:
            raise ValueError("truncated histogram header; is it a MerquryFK .cni?")
        _k, low, high = struct.unpack_from("<3i", data, offset)
        n = high - low + 1
        if low < 1 or n < 1 or offset + 28 + 8 * n > len(data):
            raise ValueError("unexpected histogram layout; is it a MerquryFK .cni?")
        histograms.append((low, high, struct.unpack_from(f"<{n}q", data, offset + 28)))
        offset += 28 + 8 * n
    if len(histograms) < len(CN_CLASSES):
        raise ValueError(f"{len(histograms)} histograms where at least {len(CN_CLASSES)} were expected")
    limits = {high for _, high, _ in histograms[:len(CN_CLASSES)]}
    if len(limits) != 1:
        raise ValueError("the copy classes cover different ranges")
    spectrum = {label: {low + i: count for i, count in enumerate(counts) if low + i < high}
                for label, (low, high, counts) in zip(CN_CLASSES, histograms)}
    return spectrum, limits.pop()


def read_cn_spectrum(path):
    """(spectrum, limit, kind) from either file type, told apart by content (a .cni
    is binary and starts with NUL bytes; a .hist is text), so any file name works."""
    with open(path, "rb") as fh:
        head = fh.read(64)
    if b"\x00" in head:
        return read_cn_cni(path) + ("MerquryFK .cni",)
    return read_cn_hist(path) + ("Merqury .hist",)


def false_duplications(spectrum, limit=None):
    """False duplications (%) as Merqury's false_duplications.sh: below 1.5x the
    single-copy peak, the share of assembly k-mers found 2+ times among those found
    1+ times. Returns (percent, cutoff); ValueError if it can't be computed."""
    single = spectrum.get("1", {})
    peak_count, peak_mult = 0, None
    for mult in sorted(single):
        if single[mult] > peak_count:
            peak_count, peak_mult = single[mult], mult
    if peak_mult is None:
        raise ValueError("no single-copy k-mers")
    cutoff = round(peak_mult * 1.5)       # half to even, like awk's printf "%.0f"
    if limit is not None and cutoff > limit:
        raise ValueError(f"the cutoff ({cutoff}) lies beyond the histogram, which stops at {limit}")

    def below(label):
        return sum(count for mult, count in spectrum.get(label, {}).items() if mult < cutoff)
    one_copy = below("1")
    more_copies = sum(below(label) for label in ("2", "3", "4", ">4"))
    if one_copy + more_copies == 0:
        raise ValueError("no assembly k-mers below the cutoff")
    return 100 * more_copies / (one_copy + more_copies), cutoff


def find_cn_spectra(merqury_dir, prefix, num_assemblies):
    """Per-assembly copy-number spectra written by GEP2's Merqury/MerquryFK runs:
    <prefix>.asm<N>.spectra-cn.hist (Merqury) or .cni (MerquryFK -k), N = 1..n, in
    the column order (GEP2 links the assemblies as asm1, asm2). A haploid run without
    the assembly name (<prefix>.spectra-cn.*) is accepted too. Missing -> None."""
    paths = [None] * num_assemblies
    if not merqury_dir or not os.path.isdir(merqury_dir) or not prefix:
        return paths
    for i in range(num_assemblies):
        stems = [f"{prefix}.asm{i + 1}.spectra-cn"]
        if num_assemblies == 1:
            stems.append(f"{prefix}.spectra-cn")
        for stem in stems:
            found = [os.path.join(merqury_dir, stem + ext) for ext in (".hist", ".cni")
                     if os.path.isfile(os.path.join(merqury_dir, stem + ext))]
            if found:
                paths[i] = found[0]
                break
    return paths


def compute_false_duplications(cn_paths):
    """False duplications (%) per assembly from find_cn_spectra(); None where the
    spectrum is missing or unusable (the reason is printed)."""
    values = []
    for i, path in enumerate(cn_paths):
        if not path:
            values.append(None)
            continue
        try:
            spectrum, limit, kind = read_cn_spectrum(path)
            percent, cutoff = false_duplications(spectrum, limit)
        except (OSError, ValueError, struct.error) as e:
            print(f"Warning: could not compute k-mer false duplications for asm{i + 1} "
                  f"from {path}: {e}")
            values.append(None)
            continue
        print(f"K-mer false duplications asm{i + 1}: {percent:.2f}% "
              f"(read multiplicity cutoff {cutoff}; {kind} {path})")
        values.append(percent)
    return values


def find_merqury_plots(merqury_dir, prefix=None):
    """Merqury/MerquryFK spectra plots (.fl preferred, then .ln/.st). Matched by exact
    name from the run prefix: <prefix>.spectra-cn.* (combined), <prefix>.asm<N>.spectra-cn.*
    (assembly N), so a prefix like ilXxx1_hap1 isn't misread; else by name pattern."""
    plots = {
        'spectra_cn': [],
        'spectra_asm': None,
        'spectra_cn_combined': None
    }

    if not merqury_dir or not os.path.isdir(merqury_dir):
        return plots

    def _first_glob(*patterns):
        for pat in patterns:
            hits = sorted(glob.glob(os.path.join(merqury_dir, pat)))
            if hits:
                return hits
        return []

    styles = ("fl", "ln", "st")

    # Exact names from the run prefix
    if prefix:
        for style in styles:
            per_asm = {}
            combined = None
            for path in glob.glob(os.path.join(merqury_dir, f"*.spectra-cn.{style}.png")):
                stem = os.path.basename(path)[:-len(f".spectra-cn.{style}.png")]
                if stem == prefix:
                    combined = path
                else:
                    m = re.fullmatch(re.escape(prefix) + r"\.asm(\d+)", stem)
                    if m:
                        per_asm[int(m.group(1))] = path
            if per_asm or combined:
                plots['spectra_cn'] = [per_asm[n] for n in sorted(per_asm)]
                plots['spectra_cn_combined'] = combined
                break
        for style in styles:
            exact = os.path.join(merqury_dir, f"{prefix}.spectra-asm.{style}.png")
            if os.path.isfile(exact):
                plots['spectra_asm'] = exact
                break

    # Fallback: name-pattern rule (a tag like .asm1. / _hap1_ marks a per-assembly plot)
    if not plots['spectra_cn'] and not plots['spectra_cn_combined']:
        cn_files = _first_glob("*.spectra-cn.fl.png", "*.spectra-cn.ln.png", "*.spectra-cn.st.png")
        for cn_file in cn_files:
            if re.search(r'[\._](asm|hap)\d+[\._]', os.path.basename(cn_file)):
                plots['spectra_cn'].append(cn_file)
            else:
                plots['spectra_cn_combined'] = cn_file
    if not plots['spectra_asm']:
        asm_files = _first_glob("*.spectra-asm.fl.png", "*.spectra-asm.ln.png", "*.spectra-asm.st.png")
        if asm_files:
            plots['spectra_asm'] = asm_files[0]

    # A single copy-number plot with no assembly tag is the assembly's own plot
    if not plots['spectra_cn'] and plots['spectra_cn_combined']:
        plots['spectra_cn'].append(plots['spectra_cn_combined'])
        plots['spectra_cn_combined'] = None

    return plots


# Star tiers per metric: (direction, (****, ***-, **--)); 'ge' higher is better,
# 'le' lower is better, past '**--' is '*---'. N50/L90 tiers: rating_tiers().
# The '#Scores legend' is printed from these tiers.
RATING_TIERS = {
    'gaps_per_gbp':          ('le', (200, 1000, 10000)),
    'compl_single':          ('ge', (95, 90, 80)),
    'compl_dupl':            ('le', (2, 5, 10)),
    'kmer_dupl':             ('le', (2, 5, 10)),      # k-mer false duplications; the EAR flags > 5%
    'compl_frameshift':      ('le', (2, 5, 15)),
    'merqury_qv':            ('ge', (50, 40, 30)),
    'merqury_completeness':  ('ge', (95, 90, 80)),
    # ---- Hi-C ----
    'hic_mapped':            ('ge', (80, 60, 40)),    # % of read pairs mapped (both reads, chromap)
    'hic_dupl':              ('le', (10, 20, 35)),    # % duplicate pairs among mapped pairs; provisional
    'hic_cis_short':         ('le', (10, 20, 30)),    # % of cis pairs under 1 kb (self-ligation, noise)
    # Cis distance shares, of cis pairs (provisional): former all-pairs tiers / 0.8
    # (assumed cis share), rounded up
    'hic_cis_10kb':          ('ge', (60, 50, 40)),    # % of cis pairs at least 10 kb apart
    'hic_cis_40kb':          ('ge', (50, 40, 30)),    # % of cis pairs at least 40 kb apart
}

# Small chromosomes (this should be double-checked!!): the cis-distance tiers assume
# ~HIC_REF_CHR chromosomes and are scaled down by hic_cis_scale(). Trans is never
# rated (small chromosomes, fragmentation and random ligation all raise it).
HIC_REF_CHR = 50_000_000
HIC_CIS_DIST = {'hic_cis_10kb': 10_000, 'hic_cis_40kb': 40_000}


def _cis_share_beyond(chr_len, dist, min_dist=1_000):
    """Share of cis contacts at least `dist` apart, among those at least `min_dist`
    apart, on a chromosome of chr_len: contact probability ~ 1/s, and (chr_len - s)
    position pairs at separation s."""
    def mass(a, b):
        return chr_len * math.log(b / a) - (b - a)
    if chr_len <= dist:
        return 0.0
    return mass(dist, chr_len) / mass(min_dist, chr_len)


def hic_cis_scale(metric_type, expected_chr_size):
    """Factor (<= 1) for the cis-distance tiers of a species with mean chromosome size
    expected_chr_size; 1 when unknown or not below HIC_REF_CHR, None when the
    chromosomes are too small (under 10x the distance) to rate it at all."""
    dist = HIC_CIS_DIST[metric_type]
    if not expected_chr_size or expected_chr_size >= HIC_REF_CHR:
        return 1.0
    if expected_chr_size < 10 * dist:
        return None
    return _cis_share_beyond(expected_chr_size, dist) / _cis_share_beyond(HIC_REF_CHR, dist)


# Rated metrics seen while building the report, in order of first use:
# metric_type -> (spec or None, note). Printed as the '#Scores legend' block.
RATINGS_USED = {}


def rating_tiers(metric_type, haploid_number=None, expected_chr_size=None):
    """((direction, (****, ***-, **--)) or None, note) for a metric. The spec is None
    when the metric can't be rated here (unknown metric, L90 without a haploid number,
    Hi-C contact metrics on very small chromosomes); the note says why, or how the
    tiers were adapted to the species ('' when they weren't)."""
    if metric_type in ('scaffold_n50', 'contig_n50'):
        # Tiers are capped at chromosome scale (EBP "C"); see N50_BASE_TIERS
        tiers = n50_tiers(metric_type, expected_chr_size)
        note = ""
        if tuple(tiers) != tuple(N50_BASE_TIERS[metric_type]):
            note = f"capped by mean chromosome {_bp_short(expected_chr_size)}"
        return ('ge', tiers), note
    if metric_type == 'l90_haploid':
        tiers = l90_tiers(haploid_number)
        if tiers is None:
            return None, "no haploid number"
        return ('le', tiers), f"haploid n = {haploid_number}"
    if metric_type in HIC_CIS_DIST:
        scale = hic_cis_scale(metric_type, expected_chr_size)
        if scale is None:
            return None, f"mean chromosome {_bp_short(expected_chr_size)} too small"
        direction, tiers = RATING_TIERS[metric_type]
        if scale < 1:
            return ((direction, tuple(round(t * scale, 1) for t in tiers)),
                    f"scaled to mean chromosome {_bp_short(expected_chr_size)}")
        return (direction, tiers), ""
    return RATING_TIERS.get(metric_type), ""


def get_rating(value, metric_type, haploid_number=None, expected_chr_size=None):
    """Star rating for a metric value. expected_chr_size (bp or None, mean chromosome
    size) caps the N50 tiers on sub-megabase chromosomes and scales the Hi-C
    cis-distance tiers (see HIC_REF_CHR)."""
    spec, note = rating_tiers(metric_type, haploid_number, expected_chr_size)
    if metric_type not in RATINGS_USED and (spec is not None or note):
        RATINGS_USED[metric_type] = (spec, note)
    if value is None or spec is None:
        return '····'

    direction, (t4, t3, t2) = spec
    better = (lambda v, t: v >= t) if direction == 'ge' else (lambda v, t: v <= t)
    if better(value, t4):
        return '****'
    elif better(value, t3):
        return '***-'
    elif better(value, t2):
        return '**--'
    else:
        return '*---'


def format_number(value):
    """Format large numbers with commas."""
    if value is None:
        return "N/A"
    if isinstance(value, (int, float)):
        if isinstance(value, float):
            if value != value or value in (float("inf"), float("-inf")):
                return str(value)          # nan/inf (e.g. Merqury QV with no errors)
            if value < 10:
                return f"{value:.4f}"
            elif value < 100:
                return f"{value:.2f}"
            else:
                return f"{int(value):,}"
        return f"{value:,}"
    return str(value)


def parse_fcs_gx(filepath):
    """
    Parse FCS-GX report file and count flagged sequences.
    Skips the first 2 header lines and counts the remaining lines.
    Originally was "awk 'NR > 2' report.txt | wc -l", but let's python
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        return max(0, len(lines) - 2)
    except Exception as e:
        print(f"Warning: Could not parse FCS-GX file {filepath}: {e}")
        return None


def parse_inspector(filepath):
    """
    Parse Inspector summary_statistics file and extract the structural error count.
    Looks for a line like: Structural error\t0
    """
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith("Structural error"):
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        return int(parts[-1])
        print(f"Warning: 'Structural error' line not found in {filepath}")
        return None
    except Exception as e:
        print(f"Warning: Could not parse Inspector file {filepath}: {e}")
        return None


class _RandomTransCounter:
    """Sums chrom_freq lines as they are read (millions on fragmented assemblies)."""

    def __init__(self):
        self.ends, self.pairs, self.usable = {}, 0, True

    def add(self, key, value):
        """Count one stats line; returns False when it is not a chrom_freq line."""
        if not key.startswith("chrom_freq/"):
            return False
        names = key[len("chrom_freq/"):].split("/")
        try:
            count = int(value)
        except ValueError:
            count = None
        if len(names) != 2 or count is None:
            self.usable = False
            return True
        self.pairs += count
        for name in names:
            self.ends[name] = self.ends.get(name, 0) + count
        return True

    def percent(self):
        if not self.usable or not self.pairs:
            return None
        return (1 - sum((n / (2 * self.pairs)) ** 2 for n in self.ends.values())) * 100


def expected_random_trans(stats):
    """% trans expected for randomly joined pairs: 1 - sum(p_i^2), p_i = each sequence's
    share of read ends (pairtools chrom_freq). None without chrom_freq or if a name has '/'."""
    counter = _RandomTransCounter()
    for key, value in stats.items():
        counter.add(key, value)
    return counter.percent()


def parse_pairtools_stats(filepath):
    """Hi-C metrics from a pairtools stats file (dedup --mark-dups --output-stats, E02/E15).
    Duplicates are over all mapped pairs; cis/trans over non-duplicate pairs. *_pct of
    cis pairs: cis_10kb, cis_40kb, junk_cis (<1 kb); trans_pct of all non-duplicate
    pairs, trans_wo_short_pct without cis <1 kb. Without duplicate marking,
    dup_rate/complexity are None."""
    metrics = {
        'total_pairs': None,
        'uu_pairs': None,
        'cis': None,
        'trans': None,
        'valid_pair_rate': None,
        'cis_trans_ratio': None,
        'trans_pct': None,
        'trans_wo_short_pct': None,
        'random_trans_pct': None,
        'cis_10kb_pct': None,
        'cis_40kb_pct': None,
        'frac_cis': None,
        'longrange_cis': None,
        'longrange_cis_40k': None,
        'junk_cis': None,
        'total_mapped': None,
        'total_dups': None,
        'total_nodups': None,
        'dup_rate': None,
        'complexity': None,
    }

    try:
        kv = {}
        random_trans = _RandomTransCounter()     # chrom_freq lines are counted, not kept
        with open(filepath, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                if not line or '\t' not in line:
                    continue
                key, _, value = line.partition('\t')
                key, value = key.strip(), value.strip()
                if not random_trans.add(key, value):
                    kv[key] = value
    except Exception as e:
        print(f"Warning: Could not parse pairtools stats file {filepath}: {e}")
        return metrics

    def _int(key):
        try:
            return int(kv[key])
        except (KeyError, ValueError):
            return None

    def _float(key):
        try:
            return float(kv[key])
        except (KeyError, ValueError):
            return None

    metrics['total_pairs'] = _int('total')
    metrics['uu_pairs'] = _int('pair_types/UU')
    metrics['cis'] = _int('cis')
    metrics['trans'] = _int('trans')
    metrics['frac_cis'] = _float('summary/frac_cis')

    # valid (unique-unique) pair rate, as % of all pairs
    if metrics['uu_pairs'] is not None and metrics['total_pairs']:
        metrics['valid_pair_rate'] = (metrics['uu_pairs'] / metrics['total_pairs']) * 100

    # cis-to-trans ratio, and trans as % of all (non-duplicate) cis + trans pairs
    if metrics['cis'] is not None and metrics['trans']:
        metrics['cis_trans_ratio'] = metrics['cis'] / metrics['trans']
    if metrics['cis'] is not None and metrics['trans'] is not None and \
            metrics['cis'] + metrics['trans']:
        metrics['trans_pct'] = metrics['trans'] / (metrics['cis'] + metrics['trans']) * 100

    metrics['random_trans_pct'] = random_trans.percent()
    cis_1kb_plus = _int('cis_1kb+')
    if metrics['trans'] is not None and cis_1kb_plus is not None and \
            metrics['trans'] + cis_1kb_plus:
        metrics['trans_wo_short_pct'] = metrics['trans'] / (metrics['trans'] + cis_1kb_plus) * 100

    # cis distance classes as % of cis pairs
    for key, column in (('cis_10kb_pct', 'cis_10kb+'), ('cis_40kb_pct', 'cis_40kb+')):
        count = _int(column)
        if metrics['cis'] and count is not None:
            metrics[key] = count / metrics['cis'] * 100

    # long-range cis (>=10kb) as % of all contacts
    frac_10kb = _float('summary/frac_cis_10kb+')
    if frac_10kb is not None:
        metrics['longrange_cis'] = frac_10kb * 100

    # very-long-range cis (>=40kb) as % of all contacts
    frac_40kb = _float('summary/frac_cis_40kb+')
    if frac_40kb is not None:
        metrics['longrange_cis_40k'] = frac_40kb * 100

    # "Junk" cis fraction = cis contacts below 1 kb (self-ligation / dangling ends),
    # as a % of all cis contacts: (cis - cis_1kb+) / cis. High -> noisier library.
    cis_1kb = _int('cis_1kb+')
    if metrics['cis'] and cis_1kb is not None:
        metrics['junk_cis'] = ((metrics['cis'] - cis_1kb) / metrics['cis']) * 100

    # Duplicates (only meaningful when the file comes from pairtools dedup)
    metrics['total_mapped'] = _int('total_mapped')
    metrics['total_dups'] = _int('total_dups')
    metrics['total_nodups'] = _int('total_nodups')
    dedup_marked = 'pair_types/DD' in kv or (metrics['total_dups'] or 0) > 0
    if dedup_marked and metrics['total_mapped']:
        metrics['dup_rate'] = metrics['total_dups'] / metrics['total_mapped'] * 100
        complexity = _float('summary/complexity_naive')
        if complexity is not None and complexity == complexity and complexity not in (float('inf'), 0):
            metrics['complexity'] = complexity

    return metrics


def parse_chromap_log(filepath):
    """Hi-C mapping from a chromap log: MAPQ threshold, read pairs ('Number of reads' / 2),
    mapped_pct (output pairs / read pairs), unique_yield (pair-level uni-mappings / read pairs)."""
    metrics = {
        'mapq_threshold': None,
        'read_pairs': None,
        'uni_pairs': None,
        'output_pairs': None,
        'unique_yield': None,
        'mapped_pct': None,
    }

    try:
        with open(filepath, 'r') as f:
            content = f.read()
    except Exception as e:
        print(f"Warning: Could not read chromap log {filepath}: {e}")
        return metrics

    mapq_match = re.search(r'MAPQ-threshold:\s*(\d+)', content)
    if mapq_match:
        metrics['mapq_threshold'] = int(mapq_match.group(1))

    # "Number of reads: N" counts individual reads; a Hi-C pair is two reads.
    reads_match = re.search(r'Number of reads:\s*(\d+)', content)
    if reads_match:
        metrics['read_pairs'] = int(reads_match.group(1)) // 2

    # Pair-level counts appear as: "# uni-mappings: X, # multi-mappings: Y, total: Z"
    pair_match = re.search(
        r'#\s*uni-mappings:\s*(\d+),\s*#\s*multi-mappings:\s*(\d+),\s*total:\s*(\d+)',
        content)
    if pair_match:
        metrics['uni_pairs'] = int(pair_match.group(1))
        metrics['output_pairs'] = int(pair_match.group(3))
    else:
        # Fallback to the explicit "Number of output mappings (passed filters): N"
        out_match = re.search(r'Number of output mappings \(passed filters\):\s*(\d+)', content)
        if out_match:
            metrics['output_pairs'] = int(out_match.group(1))

    # pair-level unique yield = unique output pairs / input read pairs
    if metrics['uni_pairs'] is not None and metrics['read_pairs']:
        metrics['unique_yield'] = (metrics['uni_pairs'] / metrics['read_pairs']) * 100

    # mapped = every output pair (unique or multi-mapped) / input read pairs
    if metrics['output_pairs'] is not None and metrics['read_pairs']:
        metrics['mapped_pct'] = (metrics['output_pairs'] / metrics['read_pairs']) * 100

    return metrics


def parse_hic_library_qc(stats_files):
    """Per-library rows from E15's <lib>.<subset>.dedup_stats.txt (+ .chromap.log), by
    library then subset size; a subset holding the whole library is shown once, as 'all'."""
    by_lib = {}
    for path in stats_files:
        name = os.path.basename(path)
        m = re.fullmatch(r"(.+)\.(\d+(?:\.\d+)?)([kKmM]?)\.dedup_stats\.txt", name)
        if not m:
            print(f"Warning: Hi-C library stats file not named <library>.<subset>.dedup_stats.txt: {path}")
            continue
        lib, number, unit = m.group(1), float(m.group(2)), m.group(3).lower()
        nominal = number * {"": 1, "k": 1e3, "m": 1e6}[unit]
        pt = parse_pairtools_stats(path)
        log = path[:-len(".dedup_stats.txt")] + ".chromap.log"
        cm = parse_chromap_log(log) if os.path.isfile(log) else {}
        if not cm:
            print(f"Warning: no chromap log beside {path}")
        pairs = cm.get('read_pairs')
        row = dict(pt)
        row.update({
            'lib': lib,
            'nominal': nominal,
            'subset': f"{m.group(2)}{m.group(3)}",
            'read_pairs': pairs,
            'unique_yield': cm.get('unique_yield'),
            'mapped_pct': cm.get('mapped_pct'),
        })
        if pairs is not None and pairs < nominal * 0.999:
            row['subset'] = "all"
        by_lib.setdefault(lib, []).append(row)

    rows = []
    for lib, lib_rows in by_lib.items():
        seen_all = False
        for row in sorted(lib_rows, key=lambda r: r['nominal']):
            if row['subset'] == "all":
                if seen_all:
                    continue
                seen_all = True
            rows.append(row)
    return rows


# Hi-C library summary: status flags (provisional)
HIC_LIB_SATURATED_DUP = 20     # % duplicates at the library's full depth
HIC_LIB_LOW_MAPPED = 60        # % of read pairs mapped
HIC_LIB_HIGH_JUNK = 20         # % of cis pairs under 1 kb (the '***-' bar of Hi-C Cis <1 kb)
HIC_LIB_SUBSET_DIFF = 2        # percentage points between subsets of one library
HIC_LIB_TARGET_DUP = 20        # the library line gives the depth at this % duplicates


def _unique_molecules(sequenced, complexity):
    """Distinct molecules expected among `sequenced` mapped pairs of a library with
    `complexity` distinct molecules (Lander-Waterman, the model pairtools uses)."""
    return complexity * (1 - math.exp(-sequenced / complexity))


def _dup_pct(sequenced, complexity):
    """% duplicates expected among `sequenced` mapped pairs (same model)."""
    if not sequenced:
        return None
    return (1 - _unique_molecules(sequenced, complexity) / sequenced) * 100


def _depth_ratio_at_dup(dup_pct):
    """Mapped pairs / complexity at which `dup_pct` % duplicates are expected (same model;
    about 0.47 for 20%). Solved by bisection: the duplicate rate rises with depth."""
    low, high = 0.0, 1.0
    while _dup_pct(high, 1.0) < dup_pct:
        high *= 2
    for _ in range(100):
        mid = (low + high) / 2
        if _dup_pct(mid, 1.0) < dup_pct:
            low = mid
        else:
            high = mid
    return (low + high) / 2


def read_hic_library_pairs(count_args):
    """{library: read pairs} from 'library=<R1 counts file>' arguments; the counts file
    is GEP2's per-file seqtk size table (file, stage, reads, bases), so R1's reads are
    the library's pairs."""
    pairs = {}
    for arg in count_args or []:
        lib, sep, path = arg.partition("=")
        if not sep:
            print(f"Warning: --hic-library-counts expects library=file, got: {arg}")
            continue
        try:
            with open(path) as fh:
                rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
            header = rows[0]
            pairs[lib] = int(rows[1][header.index("reads")])
        except (OSError, ValueError, IndexError) as e:
            print(f"Warning: could not read the read count of Hi-C library {lib} from {path}: {e}")
    return pairs


def summarise_hic_libraries(lib_rows, library_pairs=None):
    """One summary per library, from its largest subset or the whole library: lib,
    read_pairs, measured, dup_pct, pairs_at_target (for HIC_LIB_TARGET_DUP %), status."""
    library_pairs = library_pairs or {}
    by_lib = {}
    for row in lib_rows:
        by_lib.setdefault(row['lib'], []).append(row)

    summaries = []
    for lib, rows in by_lib.items():
        top = rows[-1]                      # rows come sorted by subset size, 'all' last
        measured = top['subset'] == "all"
        pairs = top.get('read_pairs') if measured else library_pairs.get(lib)
        mapped = top.get('total_mapped')
        complexity = top.get('complexity')
        # share of read pairs that map, to turn mapped pairs into read pairs and back
        mapped_share = mapped / top['read_pairs'] if mapped and top.get('read_pairs') else None
        depth = pairs * mapped_share if pairs and mapped_share else None   # mapped, whole library

        dup_pct = top.get('dup_rate') if measured else None
        lower_bound = False
        if dup_pct is None and complexity and depth:
            dup_pct = _dup_pct(depth, complexity)
        if dup_pct is None and top.get('dup_rate') is not None:
            # whole-library depth unknown (no read count): the library holds at least
            # its largest subset, so that subset's measured rate is a floor
            dup_pct, lower_bound = top['dup_rate'], True

        pairs_at_target = None
        if complexity and mapped_share:
            pairs_at_target = _depth_ratio_at_dup(HIC_LIB_TARGET_DUP) * complexity / mapped_share

        status = []
        checked = any(v is not None for v in (dup_pct, top.get('mapped_pct'), top.get('junk_cis')))
        if dup_pct is not None and dup_pct > HIC_LIB_SATURATED_DUP:
            status.append("saturated")
        if top.get('mapped_pct') is not None and top['mapped_pct'] < HIC_LIB_LOW_MAPPED:
            status.append("low mapped")
        if top.get('junk_cis') is not None and top['junk_cis'] > HIC_LIB_HIGH_JUNK:
            status.append("high junk")
        # rates counted after removing duplicates (Mapped includes them, so it should not move)
        for key in ('trans_pct', 'junk_cis', 'cis_10kb_pct'):
            values = [r[key] for r in rows if r.get(key) is not None]
            if len(values) > 1 and max(values) - min(values) > HIC_LIB_SUBSET_DIFF:
                status.append("subsets disagree")
                break

        summaries.append({
            'lib': lib, 'read_pairs': pairs, 'measured': measured, 'dup_pct': dup_pct,
            'lower_bound': lower_bound, 'subset': top['subset'],
            'pairs_at_target': pairs_at_target,
            'status': status or (["ok"] if checked else ["N/A"]),
        })
    return summaries


def _pairs_short(value, digits=3):
    """Read pairs as 2.97M / 66M / 850k (`digits` significant digits)."""
    value = float(f"{value:.{digits}g}")     # round first, so 999,999 becomes 1M, not 1000k
    for unit, size in (("G", 1e9), ("M", 1e6), ("k", 1e3)):
        if value >= size:
            scaled = value / size
            decimals = max(0, digits - len(str(int(scaled))))
            text = f"{scaled:.{decimals}f}"
            if "." in text:
                text = text.rstrip('0').rstrip('.')
            return text + unit
    return f"{value:.0f}"


def hic_library_line(sm):
    """One plain line per library: a **PROJECTION** of how deep it can be sequenced, then
    any **WARNINGS** as `tags`. A saturated library is already past that depth, so its
    line is the warnings, with its duplicate rate."""
    if sm['status'] == ["N/A"]:
        return f"{sm['lib']}: no data"
    warnings = [label for label in sm['status'] if label != "ok"]
    now = _pairs_short(sm['read_pairs']) if sm['read_pairs'] else None

    if "saturated" in warnings and sm.get('lower_bound'):
        detail = (f" (at least {sm['dup_pct']:.0f}% duplicates, already in its "
                  f"{sm['subset']} subset)")
        labels = [f"`{label}`" + (detail if label == "saturated" else "") for label in warnings]
        return f"{sm['lib']}: **WARNINGS** " + ", ".join(labels)

    if "saturated" in warnings and sm['dup_pct'] is not None:
        approx, projected = ("", "") if sm['measured'] else ("~", ", projected")
        where = f" at its {now} read pairs" if now else ""
        detail = f" ({approx}{sm['dup_pct']:.0f}% duplicates{where}{projected})"
        labels = [f"`{label}`" + (detail if label == "saturated" else "") for label in warnings]
        return f"{sm['lib']}: **WARNINGS** " + ", ".join(labels)

    now_text = f" ({now} now)" if now else ""
    if sm['pairs_at_target'] is None:
        line = f"{sm['lib']}: too few duplicates to project how deep it can be sequenced{now_text}"
    else:
        line = (f"{sm['lib']}: **PROJECTION** ~{_pairs_short(sm['pairs_at_target'], 2)} read pairs "
                f"would give {HIC_LIB_TARGET_DUP}% duplicates{now_text}")
    if warnings:
        line += "; **WARNINGS** " + ", ".join(f"`{label}`" for label in warnings)
    return line


PERCENT_METRICS = {
    'compl_single', 'compl_dupl', 'kmer_dupl', 'compl_frameshift', 'merqury_completeness',
    'hic_mapped', 'hic_dupl', 'hic_cis_short', 'hic_cis_10kb', 'hic_cis_40kb',
}


def _legend_value(metric_type, value):
    """A tier value as shown in the legend: 100 Mb / 1.53 Mb / 912 kb, 1,000, 95%, 0.5."""
    if metric_type in ('scaffold_n50', 'contig_n50'):
        return _bp_short(value)
    text = f"{value:,.2f}".rstrip('0').rstrip('.')
    return text + ("%" if metric_type in PERCENT_METRICS else "")


def scores_legend_lines(labels, haploid_number=None, expected_chr_size=None):
    """'#Scores legend' block: the star tiers of every metric rated in this report.
    labels: metric_type -> row label, in the order the legend lists them (metrics
    missing from labels come last)."""
    if not RATINGS_USED:
        return []
    order = [m for m in labels if m in RATINGS_USED] + \
            [m for m in RATINGS_USED if m not in labels]
    rows = []
    for metric_type in order:
        spec, note = RATINGS_USED[metric_type]
        label = labels.get(metric_type, metric_type)
        if spec is None:
            rows.append((label, [f"not rated ({note})"]))
            continue
        direction, tiers = spec
        better, worse = ('>=', '<') if direction == 'ge' else ('<=', '>')
        cells = [f"{better} {_legend_value(metric_type, t)}" for t in tiers]
        cells.append(f"{worse} {_legend_value(metric_type, tiers[-1])}")
        if note:
            cells.append(f"({note})")
        rows.append((label, cells))

    header = ["metric", "****", "***-", "**--", "*---"]
    tiered = [cells for _, cells in rows if len(cells) >= 4]
    widths = [max([len(header[0])] + [len(label) for label, _ in rows])]
    widths += [max([len(header[c + 1])] + [len(cells[c]) for cells in tiered]) for c in range(4)]
    lines = ["#Scores legend",
             "  ".join(f"{h:<{w}}" for h, w in zip(header, widths)).rstrip()]
    for label, cells in rows:
        if len(cells) == 1:
            lines.append(f"{label:<{widths[0]}}  {cells[0]}")
            continue
        parts = [f"{label:<{widths[0]}}"] + [f"{c:<{w}}" for c, w in zip(cells[:4], widths[1:])]
        parts += cells[4:]
        lines.append("  ".join(parts).rstrip())
    lines.append("···· = not rated")
    return lines


def generate_report(species_name, assembly_id, gfastats_list, compleasm_list, busco_list,
                   merqury_qv_values, merqury_completeness_values,
                   haploid_number, haploid_source, taxon_id, family, resolution_level,
                   genomescope_plot,
                   merqury_plots, hic_plots, blob_plots, fcs_gx_files,
                   inspector_values, output_file, kmer_tool="MERQ",
                   pairtools_list=None, chromap_list=None, expected_chr_size=None,
                   merqury_both=(None, None), kmer_false_dup=None,
                   hic_library_qc=None, version_lines=None, hic_library_pairs=None):
    """Generate the markdown report supporting 1 or 2 assemblies.

    version_lines: pipeline/container version lines; printed with the '#Scores legend'
    in one code block at the end of the report."""

    RATINGS_USED.clear()
    num_assemblies = len(gfastats_list)
    is_diploid = num_assemblies == 2
    tool_display = "MerquryFK" if kmer_tool == "MERQ.FK" else "Merqury"

    # Build table data: list of tuples (metric, [(value1, rating1), (value2, rating2), ...])
    table_data = []

    # Helper to add metric for all assemblies
    def add_metric(metric_name, values, ratings):
        row = [(format_number(v), r) for v, r in zip(values, ratings)]
        table_data.append((metric_name, row))

    # ---- gfastats metrics ----

    # Total bp
    values = [int(g.get('total_bp', 0)) if g.get('total_bp') else None for g in gfastats_list]
    ratings = ['····'] * num_assemblies
    add_metric("Total bp", values, ratings)

    # GC %
    values = [float(g.get('gc_percent', 0)) if g.get('gc_percent') else None for g in gfastats_list]
    ratings = ['····'] * num_assemblies
    add_metric("GC %", values, ratings)

    # Gaps/Gbp
    values = [g.get('gaps_per_gbp') for g in gfastats_list]
    ratings = [get_rating(v, 'gaps_per_gbp') for v in values]
    add_metric("Gaps/Gbp", values, ratings)

    # Total gap bp
    values = [int(g.get('total_gap_bp', 0)) if g.get('total_gap_bp') else None for g in gfastats_list]
    ratings = ['····'] * num_assemblies
    add_metric("Total gap bp", values, ratings)

    # Scaffolds
    values = [int(g.get('scaffolds', 0)) if g.get('scaffolds') else None for g in gfastats_list]
    ratings = ['····'] * num_assemblies
    add_metric("Scaffolds", values, ratings)

    # Scaffold N50
    values = [int(g.get('scaffold_n50', 0)) if g.get('scaffold_n50') else None for g in gfastats_list]
    ratings = [get_rating(v, 'scaffold_n50', expected_chr_size=expected_chr_size) for v in values]
    add_metric("Scaffold N50", values, ratings)

    # Scaffold L50
    values = [int(g.get('scaffold_l50', 0)) if g.get('scaffold_l50') else None for g in gfastats_list]
    ratings = ['····'] * num_assemblies
    add_metric("Scaffold L50", values, ratings)

    # Scaffold L90
    values = [int(g.get('scaffold_l90', 0)) if g.get('scaffold_l90') else None for g in gfastats_list]
    ratings = [get_rating(v, 'l90_haploid', haploid_number) for v in values]
    add_metric("Scaffold L90 (‡)", values, ratings)

    # Contigs
    values = [int(g.get('contigs', 0)) if g.get('contigs') else None for g in gfastats_list]
    ratings = ['····'] * num_assemblies
    add_metric("Contigs", values, ratings)

    # Contig N50
    values = [int(g.get('contig_n50', 0)) if g.get('contig_n50') else None for g in gfastats_list]
    ratings = [get_rating(v, 'contig_n50', expected_chr_size=expected_chr_size) for v in values]
    add_metric("Contig N50", values, ratings)

    # Contig L50
    values = [int(g.get('contig_l50', 0)) if g.get('contig_l50') else None for g in gfastats_list]
    ratings = ['····'] * num_assemblies
    add_metric("Contig L50", values, ratings)

    # Contig L90
    values = [int(g.get('contig_l90', 0)) if g.get('contig_l90') else None for g in gfastats_list]
    ratings = [get_rating(v, 'l90_haploid', haploid_number) for v in values]
    add_metric("Contig L90   (‡)", values, ratings)

    # ---- compleasm metrics (skip if not available) ----

    # Check if any compleasm data is available
    has_compleasm = any(
        c.get('eukaryota_single') is not None or c.get('other_single') is not None
        for c in compleasm_list
    )

    if has_compleasm:
        # COMPL Sing. (e) - eukaryota
        values = [c.get('eukaryota_single') for c in compleasm_list]
        if any(v is not None for v in values):
            ratings = [get_rating(v, 'compl_single') for v in values]
            formatted_values = [f"{v:.2f}%" if v is not None else "N/A" for v in values]
            add_metric("COMPL Sing.  (e)", formatted_values, ratings)

        # COMPL Dupl. (e)
        values = [c.get('eukaryota_dupl') for c in compleasm_list]
        if any(v is not None for v in values):
            ratings = [get_rating(v, 'compl_dupl') for v in values]
            formatted_values = [f"{v:.2f}%" if v is not None else "N/A" for v in values]
            add_metric("COMPL Dupl.  (e)", formatted_values, ratings)

        # COMPL Sing. (x) - other lineage
        values = [c.get('other_single') for c in compleasm_list]
        if any(v is not None for v in values):
            ratings = [get_rating(v, 'compl_single') for v in values]
            formatted_values = [f"{v:.2f}%" if v is not None else "N/A" for v in values]
            add_metric("COMPL Sing.  (x)", formatted_values, ratings)

        # COMPL Dupl. (x)
        values = [c.get('other_dupl') for c in compleasm_list]
        if any(v is not None for v in values):
            ratings = [get_rating(v, 'compl_dupl') for v in values]
            formatted_values = [f"{v:.2f}%" if v is not None else "N/A" for v in values]
            add_metric("COMPL Dupl.  (x)", formatted_values, ratings)

        # COMPL Frame. - frameshift rate (only from full_table.tsv)
        values = [c.get('frameshift_rate') for c in compleasm_list]
        if any(v is not None for v in values):
            # Determine label based on lineage
            is_euk = any(c.get('frameshift_is_eukaryota') for c in compleasm_list)
            frame_label = "COMPL Frame. (e)" if is_euk else "COMPL Frame. (x)"
            ratings = [get_rating(v, 'compl_frameshift') for v in values]
            formatted_values = [f"{v:.2f}%" if v is not None else "N/A" for v in values]
            add_metric(frame_label, formatted_values, ratings)

    # ---- BUSCO metrics (skip if not available) ----

    # Check if any BUSCO data is available
    has_busco = any(
        b.get('eukaryota_single') is not None or b.get('other_single') is not None
        for b in busco_list
    )

    if has_busco:
        # BUSCO Sing. (e) - eukaryota
        values = [b.get('eukaryota_single') for b in busco_list]
        if any(v is not None for v in values):
            ratings = [get_rating(v, 'compl_single') for v in values]
            formatted_values = [f"{v:.2f}%" if v is not None else "N/A" for v in values]
            add_metric("BUSCO Sing.  (e)", formatted_values, ratings)

        # BUSCO Dupl. (e)
        values = [b.get('eukaryota_dupl') for b in busco_list]
        if any(v is not None for v in values):
            ratings = [get_rating(v, 'compl_dupl') for v in values]
            formatted_values = [f"{v:.2f}%" if v is not None else "N/A" for v in values]
            add_metric("BUSCO Dupl.  (e)", formatted_values, ratings)

        # BUSCO Sing. (x) - other lineage
        values = [b.get('other_single') for b in busco_list]
        if any(v is not None for v in values):
            ratings = [get_rating(v, 'compl_single') for v in values]
            formatted_values = [f"{v:.2f}%" if v is not None else "N/A" for v in values]
            add_metric("BUSCO Sing.  (x)", formatted_values, ratings)

        # BUSCO Dupl. (x)
        values = [b.get('other_dupl') for b in busco_list]
        if any(v is not None for v in values):
            ratings = [get_rating(v, 'compl_dupl') for v in values]
            formatted_values = [f"{v:.2f}%" if v is not None else "N/A" for v in values]
            add_metric("BUSCO Dupl.  (x)", formatted_values, ratings)

    # ---- Merqury metrics ----

    if merqury_qv_values and any(v is not None for v in merqury_qv_values):
        ratings = [get_rating(v, 'merqury_qv') for v in merqury_qv_values]
        add_metric(f"{kmer_tool} QV", merqury_qv_values, ratings)

    if merqury_completeness_values and any(v is not None for v in merqury_completeness_values):
        if is_diploid:
            # One haplotype lacks the other's heterozygous k-mers (~85% is normal), so only the
            # combined value is rated (as make_EAR.py)
            ratings = ['····'] * num_assemblies
        else:
            ratings = [get_rating(v, 'merqury_completeness') for v in merqury_completeness_values]
        add_metric(f"{kmer_tool} Compl.", merqury_completeness_values, ratings)

    # K-mer false duplications (Merqury's false_duplications.sh method)
    if kmer_false_dup and any(v is not None for v in kmer_false_dup):
        ratings = [get_rating(v, 'kmer_dupl') for v in kmer_false_dup]
        formatted_values = [f"{v:.2f}%" if v is not None else "N/A" for v in kmer_false_dup]
        add_metric(f"{kmer_tool} Dupl.", formatted_values, ratings)

    # ---- Inspector metrics ----

    if inspector_values and any(v is not None for v in inspector_values):
        ratings = ['····'] * num_assemblies
        add_metric("INSP Str. Error", inspector_values, ratings)

    # ---- Build table ----

    # Calculate column widths
    metric_width = max(len("metric"), max(len(row[0]) for row in table_data))

    # Value and rating widths for each assembly
    value_widths = []
    for i in range(num_assemblies):
        max_val_width = max(len(f"asm{i+1} value"),
                          max(len(str(row[1][i][0])) for row in table_data))
        value_widths.append(max_val_width)

    rating_width = len("asm1 rating")

    # Build header
    header_parts = [f"{'metric':<{metric_width}}"]
    separator_parts = ['-' * metric_width]

    for i in range(num_assemblies):
        header_parts.append(f"{'asm' + str(i+1) + ' value':>{value_widths[i]}}")
        header_parts.append(f"{'asm' + str(i+1) + ' rating':<{rating_width}}")
        separator_parts.append('-' * (value_widths[i] - 1) + ':')  # Right align values
        separator_parts.append('-' * rating_width)

    header = "| " + " | ".join(header_parts) + " |"
    separator = "| " + " | ".join(separator_parts) + " |"

    # Build rows
    table_rows = []
    for metric, asm_data in table_data:
        row_parts = [f"{metric:<{metric_width}}"]
        for i, (value, rating) in enumerate(asm_data):
            row_parts.append(f"{str(value):>{value_widths[i]}}")
            rating_str = f"`{rating}`"
            row_parts.append(f"{rating_str:<{rating_width}}")  # Left-align and pad rating
        table_rows.append("| " + " | ".join(row_parts) + " |")



    # Get other lineage name
    other_lineage = None
    for c in compleasm_list:
        if c.get('other_lineage'):
            other_lineage = c['other_lineage']
            break

    # Prepare haploid number information
    if haploid_number is not None:
        if resolution_level == 'genus':
            haploid_info = f"‡ = Haploid number is {haploid_number} (genus-level estimate, {haploid_source}, [GoaT](https://goat.genomehubs.org))<br>"
        else:
            haploid_info = f"‡ = Haploid number is {haploid_number} ({haploid_source}, [GoaT](https://goat.genomehubs.org))<br>"
    else:
        haploid_info = "‡ = Haploid number not found on [GoaT](https://goat.genomehubs.org)<br>"

    # Prepare header info
    if resolution_level == 'genus':
        taxon_id_str = f"{taxon_id}, genus-level" if taxon_id is not None else "not found"
    else:
        taxon_id_str = str(taxon_id) if taxon_id is not None else "not found"
    family_str = family if family is not None else "not found"

    # Build the complete report
    report_lines = [
        "GEP2 genome stats report",
        "---",
        f"### {species_name} (ID: {taxon_id_str})",
        f"##### {family_str}",
        f"#### {assembly_id}",
        "",
        header,
        separator,
    ]
    report_lines.extend(table_rows)

    # Add legend only for metrics that were included
    report_lines.append("")

    # Check if Eukaryota data was actually generated by either tool
    has_eukaryota = any(
        m.get('eukaryota_single') is not None
        for m in compleasm_list + busco_list
    )

    # Grab dynamic names (check compleasm first, then busco, mirroring other_lineage below)
    euka_lineage = next((m.get('eukaryota_lineage') for m in compleasm_list + busco_list if m.get('eukaryota_lineage')), "eukaryota (odb12.2)")

    other_lineage = next((m.get('other_lineage') for m in compleasm_list if m.get('other_lineage')), None)
    if not other_lineage:
        other_lineage = next((m.get('other_lineage') for m in busco_list if m.get('other_lineage')), None)

    # Print legend ONLY for the lineages that were actually used
    if has_compleasm or has_busco:
        if has_eukaryota:
            report_lines.append(f"e = {euka_lineage}<br>")
        if other_lineage:
            report_lines.append(f"x = {other_lineage}<br>")

        # Different lineages per assembly (e.g. LINEAGE: auto): say so, columns aren't comparable
        gene_list = compleasm_list if has_compleasm else busco_list
        per_asm_lineages = []
        for m in gene_list[:num_assemblies]:
            names = [n for n in (m.get('eukaryota_lineage'), m.get('other_lineage')) if n]
            per_asm_lineages.append(" + ".join(names) if names else "N/A")
        if len(set(per_asm_lineages)) > 1:
            details = ", ".join(f"asm{i + 1}: {n}" for i, n in enumerate(per_asm_lineages))
            report_lines.append(f"**Warning!** Gene-completeness lineages differ across "
                                f"assemblies ({details}); their values are not directly comparable<br>")
            print(f"WARNING: gene-completeness lineages differ across assemblies ({details})")

    report_lines.append(haploid_info)

    # K-mer tool, and the combined asm1+asm2 values (only present with 2 assemblies).
    # With two assemblies the combined completeness is the one rated (see the table).
    has_kmer_rows = any(v is not None for v in (merqury_qv_values or []) + (merqury_completeness_values or [])
                        + (kmer_false_dup or []))
    if has_kmer_rows:
        report_lines.append(f"{kmer_tool} = {tool_display} (QV, k-mer completeness"
                            + (", k-mer false duplications" if kmer_false_dup and
                               any(v is not None for v in kmer_false_dup) else "")
                            + ")<br>")
    both_qv, both_compl = merqury_both
    if both_qv is not None:
        report_lines.append(f"{kmer_tool} QV for both: {format_number(both_qv)} "
                            f"`{get_rating(both_qv, 'merqury_qv')}`<br>")
    if both_compl is not None:
        report_lines.append(f"{kmer_tool} Compl. for both: {format_number(both_compl)} "
                            f"`{get_rating(both_compl, 'merqury_completeness')}`<br>")

    # Add GenomeScope2 section if plot exists
    if genomescope_plot and os.path.exists(genomescope_plot):
        # Calculate relative path from report location
        report_dir = os.path.dirname(os.path.abspath(output_file))
        gs_rel_path = os.path.relpath(genomescope_plot, report_dir)

        report_lines.extend([
            "",
            "---",
            "### GenomeScope2 Profiling",
            "",
            f"![GenomeScope2 Linear Plot]({gs_rel_path})",
        ])

    # Add Merqury plots section
    if merqury_plots:
        has_plots = (merqury_plots.get('spectra_asm') or
                    merqury_plots.get('spectra_cn') or
                    merqury_plots.get('spectra_cn_combined'))

        if has_plots:
            report_dir = os.path.dirname(os.path.abspath(output_file))

            report_lines.extend([
                "",
                "---",
                f"### {tool_display} Plots",
            ])

            # Spectra-asm plot (combined assembly spectra)
            if merqury_plots.get('spectra_asm'):
                rel_path = os.path.relpath(merqury_plots['spectra_asm'], report_dir)
                report_lines.extend([
                    "",
                    "#### Assembly Spectra",
                    f"![Spectra ASM]({rel_path})",
                ])

            # Individual spectra-cn plots
            if merqury_plots.get('spectra_cn'):
                report_lines.extend([
                    "",
                    "#### Copy Number Spectra",
                ])
                for i, cn_plot in enumerate(merqury_plots['spectra_cn']):
                    rel_path = os.path.relpath(cn_plot, report_dir)
                    label = f"asm{i+1}" if len(merqury_plots['spectra_cn']) > 1 else "Assembly"
                    report_lines.append(f"![Spectra CN - {label}]({rel_path})")
                    report_lines.append("")

            # Combined spectra-cn plot (diploid only)
            if merqury_plots.get('spectra_cn_combined'):
                rel_path = os.path.relpath(merqury_plots['spectra_cn_combined'], report_dir)
                report_lines.extend([
                    "#### Combined Copy Number Spectra",
                    f"![Spectra CN Combined]({rel_path})",
                    "",
                ])

    # Add Hi-C Metrics section
    pt_list = pairtools_list or []
    cm_list = chromap_list or []

    def _hic_get_entry(lst, i):
        return lst[i] if i < len(lst) and lst[i] else {}

    def _hic_get(lst, i, key):
        return _hic_get_entry(lst, i).get(key)

    def _millions(v):
        return v / 1e6 if v is not None else None

    # Funnel, each row over the one above: read pairs -> mapped -> non-duplicate ->
    # trans/cis -> cis distance. (label, metric_type, [(shown, rated)] per assembly)
    def _count_pct(count, pct, pct_fmt='{:.1f}%'):
        if pct is None:
            return (None, None)
        text = pct_fmt.format(pct)
        return (f"{count:,} ({text})" if count is not None else text, pct)

    def _plain(v, fmt):
        return (fmt.format(v), v) if v is not None else (None, None)

    def _trans_cell(trans, wo_short, random):
        if trans is None:
            return (None, None)
        extra = []
        if wo_short is not None:
            extra.append(f"w/o <1kb: {wo_short:.1f}%")
        if random is not None:
            extra.append(f"random: ~{random:.0f}%")
        text = f"{trans:.2f}%" + (f" ({'; '.join(extra)})" if extra else "")
        return (text, trans)

    hic_specs = [
        ("Total Pairs", None,
         [_plain(_hic_get(cm_list, i, 'read_pairs'), '{:,}') for i in range(num_assemblies)]),
        ("Mapped", 'hic_mapped',
         [_count_pct(_hic_get(cm_list, i, 'output_pairs'), _hic_get(cm_list, i, 'mapped_pct'))
          for i in range(num_assemblies)]),
        ("Dupl.", 'hic_dupl',
         [_count_pct(_hic_get(pt_list, i, 'total_dups'), _hic_get(pt_list, i, 'dup_rate'), '{:.2f}%')
          for i in range(num_assemblies)]),
        ("Trans", None,
         [_trans_cell(_hic_get(pt_list, i, 'trans_pct'), _hic_get(pt_list, i, 'trans_wo_short_pct'),
                      _hic_get(pt_list, i, 'random_trans_pct'))
          for i in range(num_assemblies)]),
        ("Cis <1 kb", 'hic_cis_short',
         [_plain(_hic_get(pt_list, i, 'junk_cis'), '{:.2f}%') for i in range(num_assemblies)]),
        ("Cis ≥10 kb", 'hic_cis_10kb',
         [_plain(_hic_get(pt_list, i, 'cis_10kb_pct'), '{:.2f}%') for i in range(num_assemblies)]),
        ("Cis ≥40 kb", 'hic_cis_40kb',
         [_plain(_hic_get(pt_list, i, 'cis_40kb_pct'), '{:.2f}%') for i in range(num_assemblies)]),
        ("Complexity", None,
         [_plain(_millions(_hic_get(pt_list, i, 'complexity')), '{:,.1f}M')
          for i in range(num_assemblies)]),
    ]

    hic_rows = []
    for label, mtype, cells in hic_specs:
        if any(shown is not None for shown, _ in cells):
            row = [(shown if shown is not None else "N/A",
                    get_rating(value, mtype, expected_chr_size=expected_chr_size))
                   for shown, value in cells]
            hic_rows.append((label, row))

    if hic_rows:
        hic_metric_width = max(len("metric"), max(len(r[0]) for r in hic_rows))
        hic_value_widths = [
            max(len(f"asm{i+1} value"), max(len(str(r[1][i][0])) for r in hic_rows))
            for i in range(num_assemblies)
        ]
        hic_rating_width = len("asm1 rating")

        hic_hdr = [f"{'metric':<{hic_metric_width}}"]
        hic_sep = ['-' * hic_metric_width]
        for i in range(num_assemblies):
            hic_hdr.append(f"{'asm' + str(i+1) + ' value':>{hic_value_widths[i]}}")
            hic_hdr.append(f"{'asm' + str(i+1) + ' rating':<{hic_rating_width}}")
            hic_sep.append('-' * (hic_value_widths[i] - 1) + ':')
            hic_sep.append('-' * hic_rating_width)

        report_lines.extend([
            "",
            "---",
            "### Hi-C Metrics",
            "",
            "| " + " | ".join(hic_hdr) + " |",
            "| " + " | ".join(hic_sep) + " |",
        ])
        for label, asm_data in hic_rows:
            hic_parts = [f"{label:<{hic_metric_width}}"]
            for i, (value, rating) in enumerate(asm_data):
                hic_parts.append(f"{str(value):>{hic_value_widths[i]}}")
                hic_parts.append(f"{('`' + rating + '`'):<{hic_rating_width}}")
            report_lines.append("| " + " | ".join(hic_parts) + " |")

        # Short help, one line per row, as a blockquote; star thresholds are in the
        # '#Scores legend' only
        mapqs = sorted({m for m in (_hic_get(cm_list, i, 'mapq_threshold')
                                    for i in range(num_assemblies)) if m})   # MAPQ 0 = no filter
        mapq_text = f" (MAPQ≥{', '.join(str(m) for m in mapqs)})" if mapqs else ""
        help_lines = [
            f"Mapped = both reads{mapq_text}, of total pairs; low values point to contamination, an "
            "incomplete reference or poor reads",
            "Dupl. = PCR/optical duplicates, of mapped pairs; rises with depth, a high value "
            "points to a library close to exhausted",
            "Trans = non-duplicate read pairs on different sequences (cis=the rest); "
            "many small chromosomes, a fragmented assembly or random ligation can raise it. "
            "(w/o <1kb = same without cis under 1kb, which pulls trans down; "
            "random = expected if pairs were joined at random, given how the reads spread "
            "over the sequences. The further w/o <1kb is below random, the more real contact "
            "signal)",
            "Cis <1kb = of cis pairs; noise, mostly self-ligation and dangling ends",
            "Cis ≥10kb = of cis pairs; fewer on small chromosomes (same for ≥40kb)",
            "Complexity = total distinct molecules estimated to be in the library, based on "
            "the % of duplicates among mapped pairs; complexity/2 mapped pairs gives about "
            "20% duplicates (a read is a duplicate when it picks a molecule already read, so "
            "this chance grows from 0% at the start to ~40% by complexity/2, averaging ~20% "
            "duplicates over all reads)",
        ]

        # Small chromosomes: say why the cis-distance stars use lower tiers (or none)
        not_rated, scaled = [], []
        for label, mtype in (("Cis ≥10 kb", 'hic_cis_10kb'), ("Cis ≥40 kb", 'hic_cis_40kb')):
            spec, note = rating_tiers(mtype, expected_chr_size=expected_chr_size)
            if spec is None:
                not_rated.append(label)
            elif note:
                scaled.append(label)
        if not_rated or scaled:
            parts = []
            if scaled:
                parts.append(" and ".join(scaled) + " stars on lower tiers")
            if not_rated:
                parts.append(" and ".join(not_rated) + " not rated")
            help_lines.append(f"{'; '.join(parts)}: mean chromosome "
                              f"{_bp_short(expected_chr_size)} (see #Scores legend)")

        report_lines.append("")
        report_lines.extend(f"> {line}<br>" for line in help_lines)

    # Per-library check: random read-pair subsets of each library, mapped to asm1
    lib_rows = list(hic_library_qc or [])
    # Single library: asm1's main mapping is that whole library (same reads and
    # settings), so it becomes the 'all' row: measured, not projected
    if len({r['lib'] for r in lib_rows}) == 1 and \
            not any(r['subset'] == "all" for r in lib_rows):
        pt_main, cm_main = _hic_get_entry(pt_list, 0), _hic_get_entry(cm_list, 0)
        if pt_main.get('total_mapped') is not None and cm_main.get('read_pairs'):
            whole = dict(pt_main)
            whole.update({
                'lib': lib_rows[0]['lib'],
                'nominal': float('inf'),
                'subset': "all",
                'read_pairs': cm_main.get('read_pairs'),
                'unique_yield': cm_main.get('unique_yield'),
                'mapped_pct': cm_main.get('mapped_pct'),
            })
            lib_rows.append(whole)
    if lib_rows:
        def _fmt(v, f):
            return f.format(v) if v is not None else "N/A"

        # same order as the rows of the Hi-C Metrics table
        header = ["library", "subset", "total pairs", "mapped", "dupl.",
                  "trans", "cis <1 kb", "cis ≥10 kb", "cis ≥40 kb", "complexity"]
        body = []
        for r in lib_rows:
            body.append([
                r['lib'], r['subset'], _fmt(r.get('read_pairs'), '{:,}'),
                _fmt(r.get('mapped_pct'), '{:.2f}%'), _fmt(r.get('dup_rate'), '{:.2f}%'),
                _fmt(r.get('trans_pct'), '{:.2f}%'), _fmt(r.get('junk_cis'), '{:.2f}%'),
                _fmt(r.get('cis_10kb_pct'), '{:.2f}%'), _fmt(r.get('cis_40kb_pct'), '{:.2f}%'),
                _fmt(r['complexity'] / 1e6 if r.get('complexity') is not None else None, '{:,.1f}M'),
            ])
        widths = [max(len(h), max(len(row[c]) for row in body)) for c, h in enumerate(header)]
        report_lines.extend([
            "",
            "#### Hi-C Library Check",
            "",
            "| " + " | ".join(f"{h:<{widths[c]}}" if c < 2 else f"{h:>{widths[c]}}"
                             for c, h in enumerate(header)) + " |",
            "| " + " | ".join('-' * widths[c] if c < 2 else '-' * (widths[c] - 1) + ':'
                             for c in range(len(header))) + " |",
        ])
        for row in body:
            report_lines.append("| " + " | ".join(f"{v:<{widths[c]}}" if c < 2 else f"{v:>{widths[c]}}"
                                                  for c, v in enumerate(row)) + " |")
        # One plain line per library (status, and how deep it can be sequenced), then how
        # to read the table, as a blockquote like the Hi-C Metrics help
        summaries = summarise_hic_libraries(lib_rows, hic_library_pairs)
        report_lines.append("")
        report_lines.extend(f"{hic_library_line(sm)}<br>" for sm in summaries)
        lines = [
            "Random read-pair subsets of each library, mapped to asm1; 'all' = the whole "
            "library (a library smaller than a subset, or the only library of this assembly, "
            "taken from the mapping above)",
            "Mapped includes duplicates, so it should be the same in all subsets (random "
            "noise at 1M pairs should be below 0.1 points)",
            "Dupl. rises with depth, while complexity should stay about the same (it often "
            "gradually rises with depth, so projections from subsets tend to inflate duplicates)",
            "Trans and cis are counted after removing duplicates, so a shift with depth "
            "points to uneven duplicates, e.g. an over-amplified library mixed into the same "
            "files (give such libraries as separate Path entries to check them one by one)",
            f"WARNINGS: `saturated` = over {HIC_LIB_SATURATED_DUP}% duplicates, points to a "
            f"library close to exhausted; `low mapped` = under {HIC_LIB_LOW_MAPPED}% of read "
            "pairs map, points to contamination, incomplete assembly or poor reads (see "
            f"the Kraken2 read screen in MultiQC, if run); `high junk` = sub-1kb pairs over "
            f"{HIC_LIB_HIGH_JUNK}% of cis pairs, points to poor digestion or ligation; "
            "`subsets disagree` = Trans, Cis <1kb or Cis ≥10kb differs by more than "
            f"{HIC_LIB_SUBSET_DIFF} points between subsets, points to uneven duplicates "
            "(see above)",
        ]
        report_lines.append("")
        report_lines.extend(f"> {line}<br>" for line in lines)

    # Add Hi-C Contact Maps section
    if hic_plots:
        # Filter to ensure files actually exist
        valid_hic_plots = [f for f in hic_plots if os.path.exists(f)]

        if valid_hic_plots:
            report_dir = os.path.dirname(os.path.abspath(output_file))

            report_lines.extend([
                "",
                "---",
                "### Hi-C Contact Maps",
            ])

            for i, hic_plot in enumerate(valid_hic_plots):
                rel_path = os.path.relpath(hic_plot, report_dir)

                # Determine label: use "asmX" if multiple files, otherwise generic "Assembly"
                if len(valid_hic_plots) > 1:
                    label = f"asm{i+1}"
                else:
                    label = "Assembly"

                report_lines.extend([
                    "",
                    f"#### {label}",
                    f"![Hi-C Map - {label}]({rel_path})",
                ])

    # Add Blob Plots and FCS-GX section
    valid_blob_plots = [f for f in blob_plots if os.path.exists(f)] if blob_plots else []
    valid_fcs_gx = [f for f in fcs_gx_files if os.path.exists(f)] if fcs_gx_files else []

    if valid_blob_plots or valid_fcs_gx:
        report_dir = os.path.dirname(os.path.abspath(output_file))

        report_lines.extend([
            "",
            "---",
            "### Contamination Screening",
        ])

        # Blob plots
        if valid_blob_plots:
            for i, blob_plot in enumerate(valid_blob_plots):
                rel_path = os.path.relpath(blob_plot, report_dir)

                # Determine label: use "asmX" if multiple files, otherwise generic "Assembly"
                if len(valid_blob_plots) > 1:
                    label = f"asm{i+1}"
                else:
                    label = "Assembly"

                report_lines.extend([
                    "",
                    f"#### {label}",
                    f"![Blob Plot - {label}]({rel_path})",
                ])

        # FCS-GX flagged sequences
        if valid_fcs_gx:
            report_lines.append("")
            for i, fcs_file in enumerate(valid_fcs_gx):
                if len(valid_fcs_gx) > 1:
                    label = f"asm{i+1}"
                else:
                    label = "Assembly"

                count = parse_fcs_gx(fcs_file)
                if count is not None:
                    report_lines.append(f"FCS-GX flagged sequences {label}: {count} <br>")
                else:
                    report_lines.append(f"FCS-GX flagged sequences {label}: N/A <br>")


    # Footer: pipeline/container versions and the star tiers used in this report
    compl_tools = [name for name, present in (("COMPL", has_compleasm), ("BUSCO", has_busco))
                   if present]
    compl_name = "/".join(compl_tools) or "COMPL"
    legend_labels = {
        'gaps_per_gbp': "Gaps/Gbp",
        'scaffold_n50': "Scaffold N50",
        'contig_n50': "Contig N50",
        'l90_haploid': "Scaffold/Contig L90",
        'compl_single': f"{compl_name} Sing.",
        'compl_dupl': f"{compl_name} Dupl.",
        'compl_frameshift': "COMPL Frame.",
        'merqury_qv': f"{kmer_tool} QV",
        'merqury_completeness': f"{kmer_tool} Compl." + (" (both)" if is_diploid else ""),
        'kmer_dupl': f"{kmer_tool} Dupl.",
        'hic_mapped': "Hi-C Mapped",
        'hic_dupl': "Hi-C Dupl.",
        'hic_cis_short': "Hi-C Cis <1 kb",
        'hic_cis_10kb': "Hi-C Cis ≥10 kb",
        'hic_cis_40kb': "Hi-C Cis ≥40 kb",
    }
    footer = [line.rstrip() for line in (version_lines or []) if line.strip()]
    legend = scores_legend_lines(legend_labels, haploid_number, expected_chr_size)
    if footer and legend:
        footer.append("")
    footer += legend
    if footer:
        report_lines.extend(["", "<br>", "<br>", "", "```"] + footer + ["```"])
    report_lines.append("")

    # Write to output file
    report = "\n".join(report_lines)

    with open(output_file, 'w') as f:
        f.write(report)

    print(f"Report generated successfully: {output_file}")


def _fmt_bp(value):
    """Compact bp for log lines: 8.21 Mb / 912 kb / 430 bp."""
    if value is None:
        return "n/a"
    value = float(value)
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f} Mb"
    if value >= 1_000:
        return f"{value / 1_000:.0f} kb"
    return f"{value:.0f} bp"


def log_run_diagnostics(goat_data, expected_chr_size, gfastats_list,
                        compleasm_list, busco_list):
    """Print everything that shapes the ratings but never appears in the report.
    """
    num_assemblies = len(gfastats_list)

    # ---- what GoaT gave us ----
    level = goat_data.get('resolution_level')
    if level == 'genus':
        print(f"[goat]    WARNING: species not found, resolved at GENUS level - "
              f"haploid number and genome size are genus estimates")
    print(f"[goat]    taxon_id={goat_data.get('taxon_id')} "
          f"lineage={goat_data.get('family')} (resolved at {level or 'n/a'} level)")
    print(f"[goat]    haploid_number={goat_data.get('haploid_number')} "
          f"(source: {goat_data.get('haploid_source') or 'n/a'})")
    gs = goat_data.get('genome_size')
    print(f"[goat]    genome_size={_fmt_bp(gs)} "
          f"(source: {goat_data.get('genome_size_source') or 'n/a'})")

    # ---- the derived scale and the thresholds it produces ----
    if expected_chr_size:
        print(f"[scaling] expected chromosome = genome_size / haploid_number "
              f"= {_fmt_bp(expected_chr_size)}")
    else:
        print("[scaling] WARNING: no expected chromosome size (GoaT genome_size or "
              "haploid_number missing) - N50 tiers stay at their absolute values "
              "and may be unreachable for a small genome")

    for mtype, label in (('scaffold_n50', 'scaffold N50'), ('contig_n50', 'contig N50  ')):
        t4, t3, t2 = n50_tiers(mtype, expected_chr_size)
        base = N50_BASE_TIERS[mtype]
        capped = " [capped at chromosome scale]" if (t4, t3, t2) != base else ""
        print(f"[scaling] {label} tiers: **** >= {_fmt_bp(t4)} | "
              f"***- >= {_fmt_bp(t3)} | **-- >= {_fmt_bp(t2)}{capped}")

    l90 = l90_tiers(goat_data.get('haploid_number'))
    if l90:
        print(f"[scaling] L90 tiers (n={goat_data.get('haploid_number')}): "
              f"**** <= {l90[0]} | ***- <= {l90[1]} | **-- <= {l90[2]}")

    # ---- per-assembly sanity checks ----
    for i, g in enumerate(gfastats_list):
        tag = f"asm{i+1}"
        total_bp = int(g['total_bp']) if g.get('total_bp') else None

        if g.get('_missing'):
            print(f"[check]   {tag}: WARNING: gfastats fields not matched "
                  f"({', '.join(g['_missing'])}) - these show as N/A with a '····' "
                  f"rating; check the gfastats version/wording in {g.get('_source')}")

        # Assembly size vs GoaT. A large mismatch means the taxon match, the
        # ploidy, or the assembly itself is off - and expected_chr_size, which
        # sets the N50 bars, is derived from that same GoaT size.
        if total_bp and gs:
            delta = (total_bp - float(gs)) / float(gs) * 100
            flag = "  <-- WARNING: >30% off, check taxon match / ploidy" if abs(delta) > 30 else ""
            print(f"[check]   {tag}: assembly {_fmt_bp(total_bp)} vs GoaT genome size "
                  f"{_fmt_bp(gs)} ({delta:+.1f}%){flag}")

        # Reachability. A tier above the assembly's own total length can never be
        # awarded, whatever the assembly quality - this check catches
        # a mis-scaled ladder even when GoaT has no genome size at all.
        if total_bp:
            for mtype, label in (('scaffold_n50', 'scaffold N50'), ('contig_n50', 'contig N50')):
                t4, t3, t2 = n50_tiers(mtype, expected_chr_size)
                for stars, bar in (('****', t4), ('***-', t3)):
                    if bar > total_bp:
                        print(f"[check]   {tag}: WARNING: {label} '{stars}' bar "
                              f"{_fmt_bp(bar)} exceeds the whole assembly "
                              f"({_fmt_bp(total_bp)}) - tier unreachable")
                if expected_chr_size and t4 > expected_chr_size:
                    print(f"[check]   {tag}: WARNING: {label} '****' bar {_fmt_bp(t4)} "
                          f"exceeds one chromosome ({_fmt_bp(expected_chr_size)})")

        gaps = g.get('gaps')
        if gaps is not None:
            print(f"[{tag}]    gfastats: {g.get('scaffolds', '?')} scaffolds / "
                  f"{g.get('contigs', '?')} contigs, {gaps} gaps "
                  f"({_fmt_bp(g.get('total_gap_bp'))} total)")

    # ---- raw counts behind the completeness percentages ----
    for i, c in enumerate(compleasm_list[:num_assemblies]):
        if not c or c.get('_counts') is None:
            continue
        counts = ", ".join(f"{k}={v}" for k, v in c['_counts'].items() if v)
        print(f"[asm{i+1}]    compleasm {c.get('_lineage_dir')}: {counts} "
              f"(n={c.get('_n_markers')})")
        if c.get('_complete_copies'):
            print(f"[asm{i+1}]    frameshifts: {c['_frameshifted_copies']}/"
                  f"{c['_complete_copies']} complete copies "
                  f"({c.get('frameshift_rate', 0):.2f}%)")


def main():
    parser = argparse.ArgumentParser(
        description='Generate GEP2 genome assembly stats report',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single assembly (haploid mode)
  python make_gep2_report.py -s Elephas_maximus -a mEleMax1_new \\
    -g stats1.txt -c compl1.txt \\
    -q merqury.qv -m merqury.completeness.stats \\
    --genomescope-plot linear_plot.png \\
    --merqury-dir kmer_stats/ \\
    -o report.md --also-pdf

  # Two assemblies (diploid mode)
  python make_gep2_report.py -s Elephas_maximus -a mEleMax1_new \\
    -g stats1.txt stats2.txt -c compl1.txt compl2.txt \\
    -q merqury.qv -m merqury.completeness.stats \\
    -o report.md
        """
    )

    parser.add_argument('-v', '--version', action='version',
                       version=f'%(prog)s {__version__}')
    parser.add_argument('-s', '--species', required=True,
                       help='Species name (with underscore, e.g., Homo_sapiens)')
    parser.add_argument('-a', '--assembly', required=True,
                       help='Assembly ID')
    parser.add_argument('-g', '--gfastats', required=True, nargs='+',
                       help='gfastats output file(s) - one per assembly')
    parser.add_argument('-c', '--compleasm-summary', required=False, nargs='+', default=[],
                       help='compleasm summary.txt file(s) - one per assembly')
    parser.add_argument('--compleasm-full', required=False, nargs='+', default=[],
                       help='compleasm {lineage}_odb{version}/full_table.tsv file(s) - one per assembly')
    parser.add_argument('--busco-summary', required=False, nargs='+', default=[],
                       help='BUSCO short_summary.txt file(s) - one per assembly')
    parser.add_argument('-q', '--merqury-qv', required=False,
                       help='Merqury QV file (*.qv)')
    parser.add_argument('-m', '--merqury-completeness', required=False,
                       help='Merqury completeness file (*.completeness.stats)')
    parser.add_argument('--merqury-dir', required=False,
                       help='Merqury output directory (for finding plots)')
    parser.add_argument('--genomescope-plot', required=False,
                       help='GenomeScope2 linear plot PNG')
    parser.add_argument('--pairtools-stats', required=False, nargs='+', default=[],
                       help='pairtools dedup stats file(s) (pairtools_dedup_stats.txt) - one per assembly')
    parser.add_argument('--chromap-log', required=False, nargs='+', default=[],
                       help='chromap log file(s) - one per assembly')
    parser.add_argument('--hic-library-stats', required=False, nargs='+', default=[],
                       help='per-library Hi-C check: <library>.<subset>.dedup_stats.txt file(s), '
                            'each with its <library>.<subset>.chromap.log beside it')
    parser.add_argument('--hic-library-counts', required=False, nargs='+', default=[],
                       help='<library>=<R1 read counts file> (seqtk size table) per Hi-C '
                            'library, for its depth in the library summary')
    parser.add_argument('--hic', required=False, nargs='+', default=[],
                       help='Hi-C contact-map png file(s) - one per assembly')
    parser.add_argument('--blob', required=False, nargs='+', default=[],
                       help='Blob Plot png file(s) - one per assembly')
    parser.add_argument('--fcs-gx', required=False, nargs='+', default=[],
                       help='FCS-GX report.txt file(s) - one per assembly')
    parser.add_argument('--Inspector', required=False, nargs='+', default=[],
                       help='Inspector summary_statistics file(s) - one per assembly')
    parser.add_argument('--also-pdf', action='store_true',
                       help='Also generate PDF (requires pandoc and weasyprint)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output markdown file')
    parser.add_argument('--versions-file', required=False,
                       help='Text file with pipeline/container version lines, printed in the '
                            'closing code block above the scores legend')

    args = parser.parse_args()

    try:
        num_assemblies = len(args.gfastats)

        if num_assemblies > 2:
            print("Error: Maximum 2 assembly files supported", file=sys.stderr)
            sys.exit(1)

        # Validate compleasm count matches gfastats count (if provided)
        compleasm_source = None
        if args.compleasm_full:
            compleasm_source = 'full'
            if len(args.compleasm_full) != num_assemblies:
                print(f"Error: Number of compleasm-full files ({len(args.compleasm_full)}) must match "
                      f"number of gfastats files ({num_assemblies})", file=sys.stderr)
                sys.exit(1)
            if args.compleasm_summary:
                print("Note: --compleasm-full provided, ignoring --compleasm-summary")
        elif args.compleasm_summary:
            compleasm_source = 'summary'
            if len(args.compleasm_summary) != num_assemblies:
                print(f"Error: Number of compleasm-summary files ({len(args.compleasm_summary)}) must match "
                      f"number of gfastats files ({num_assemblies})", file=sys.stderr)
                sys.exit(1)

        # Clean species name
        species_clean = args.species.replace('_', ' ')
        print(f"Processing species: {species_clean}")
        print(f"Number of assemblies: {num_assemblies}")

        # Get haploid number, taxon ID and family/order from GoaT
        print("Fetching data from GoaT...")
        goat_data = get_species_genomic_data_from_goat(species_clean)
        haploid_number = goat_data['haploid_number']
        haploid_source = goat_data['haploid_source']
        taxon_id = goat_data['taxon_id']
        family = goat_data['family']
        resolution_level = goat_data['resolution_level']
        genome_size = goat_data['genome_size']

        # Mean chromosome size (genome_size / haploid_number) drives the small-genome
        # star carve-out. None => carve-out disabled and scoring is unchanged.
        expected_chr_size = None
        if genome_size and haploid_number and haploid_number > 0:
            try:
                expected_chr_size = float(genome_size) / haploid_number
            except (TypeError, ValueError):
                expected_chr_size = None
        if goat_data['error']:
            print(f"Warning: Could not retrieve all GoaT data: {goat_data['error']}")
        else:
            print(f"Retrieved haploid number: {haploid_number} (source: {haploid_source})")

        if taxon_id:
            print(f"Retrieved taxon ID: {taxon_id}")
        if family:
            print(f"Retrieved family/order: {family}")

        # Parse gfastats files
        print("Parsing gfastats output(s)...")
        gfastats_list = []
        for gf_file in args.gfastats:
            gfastats_list.append(parse_gfastats(gf_file))

        # Parse compleasm files (if provided) - full takes priority over summary
        print("Parsing compleasm output(s)...")
        compleasm_list = []
        if compleasm_source == 'full':
            for comp_file in args.compleasm_full:
                compleasm_list.append(parse_compleasm_full(comp_file))
        elif compleasm_source == 'summary':
            for comp_file in args.compleasm_summary:
                compleasm_list.append(parse_compleasm(comp_file))
        else:
            # Empty compleasm results
            compleasm_list = [{'eukaryota_single': None, 'eukaryota_dupl': None,
                              'eukaryota_lineage': None,
                              'other_lineage': None, 'other_single': None,
                              'other_dupl': None}] * num_assemblies

        # Parse BUSCO files (if provided) independently
        print("Parsing BUSCO output(s)...")
        busco_list = []
        if args.busco_summary:
            for busco_file in args.busco_summary:
                busco_list.append(parse_busco_summary(busco_file))
        else:
            # Empty BUSCO results
            busco_list = [{'eukaryota_single': None, 'eukaryota_dupl': None,
                           'eukaryota_lineage': None,
                           'other_lineage': None, 'other_single': None,
                           'other_dupl': None}] * num_assemblies

        # Parse Merqury files (if provided)
        merqury_qv_values = [None] * num_assemblies
        merqury_completeness_values = [None] * num_assemblies

        if args.merqury_qv:
            print("Parsing Merqury QV file...")
            merqury_qv_values = parse_merqury_qv(args.merqury_qv, num_assemblies)

        if args.merqury_completeness:
            print("Parsing Merqury completeness file...")
            merqury_completeness_values = parse_merqury_completeness(
                args.merqury_completeness, num_assemblies)

        # Combined asm1+asm2 values ('both' rows), shown below the table
        merqury_both = (None, None)
        if num_assemblies == 2:
            merqury_both = (
                parse_merqury_both(args.merqury_qv, "qv") if args.merqury_qv else None,
                parse_merqury_both(args.merqury_completeness, "completeness")
                if args.merqury_completeness else None)

        # Merqury run prefix: <prefix>.qv names it (GEP2: the assembly ID)
        merqury_prefix = args.assembly
        for kmer_file, suffix in ((args.merqury_qv, ".qv"),
                                  (args.merqury_completeness, ".completeness.stats")):
            if kmer_file and os.path.basename(kmer_file).endswith(suffix):
                merqury_prefix = os.path.basename(kmer_file)[:-len(suffix)]
                break

        # Find Merqury plots
        merqury_plots = find_merqury_plots(args.merqury_dir, merqury_prefix)

        # K-mer false duplications from the per-assembly copy-number spectra
        cn_paths = find_cn_spectra(args.merqury_dir, merqury_prefix, num_assemblies)
        kmer_false_dup = compute_false_duplications(cn_paths)
        if args.merqury_dir and not any(cn_paths):
            print(f"Note: no {merqury_prefix}.asm<N>.spectra-cn.hist/.cni in {args.merqury_dir} "
                  f"(MerquryFK keeps .cni only with -k); k-mer false duplications not shown")

        # Identify the k-mer QV tool from the file header (Merqury vs MerquryFK)
        kmer_tool = detect_kmer_tool(args.merqury_qv, args.merqury_completeness)

        # Parse Inspector files (if provided)
        inspector_values = [None] * num_assemblies
        if args.Inspector:
            print("Parsing Inspector output(s)...")
            for i, insp_file in enumerate(args.Inspector[:num_assemblies]):
                inspector_values[i] = parse_inspector(insp_file)

        # Hi-C metrics inputs. The two files are complementary (unique yield comes
        # from chromap, the pair/cis metrics from pairtools), so they're normally given
        # together - but either alone is allowed and just yields a partial section.
        if bool(args.pairtools_stats) != bool(args.chromap_log):
            print("Note: only one of --pairtools-stats / --chromap-log provided; "
                  "the Hi-C Metrics section will be partial")
        for _name, _files in (("pairtools-stats", args.pairtools_stats),
                              ("chromap-log", args.chromap_log)):
            if _files and len(_files) != num_assemblies:
                print(f"Warning: {len(_files)} --{_name} file(s) given but {num_assemblies} "
                      f"assembly(ies); files are matched by position "
                      f"(extras ignored, missing shown as N/A)")

        # Parse Hi-C metrics files (if provided) - one per assembly
        pairtools_list = [None] * num_assemblies
        if args.pairtools_stats:
            print("Parsing pairtools stats...")
            for i, pt_file in enumerate(args.pairtools_stats[:num_assemblies]):
                pairtools_list[i] = parse_pairtools_stats(pt_file)

        chromap_list = [None] * num_assemblies
        if args.chromap_log:
            print("Parsing chromap log(s)...")
            for i, cm_file in enumerate(args.chromap_log[:num_assemblies]):
                chromap_list[i] = parse_chromap_log(cm_file)

        # Per-library Hi-C check (subsets mapped to asm1)
        hic_library_qc = []
        hic_library_pairs = read_hic_library_pairs(args.hic_library_counts)
        if args.hic_library_stats:
            print("Parsing per-library Hi-C check...")
            hic_library_qc = parse_hic_library_qc(args.hic_library_stats)

        # Everything that shapes the ratings but never reaches the report
        print("-" * 72)
        log_run_diagnostics(goat_data, expected_chr_size, gfastats_list,
                            compleasm_list, busco_list)
        print("-" * 72)

        version_lines = []
        if args.versions_file:
            with open(args.versions_file) as fh:
                version_lines = fh.read().splitlines()

        # Create output directory if needed
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

        # Generate report
        print("Generating report...")
        generate_report(
            species_clean,
            args.assembly,
            gfastats_list,
            compleasm_list,
            busco_list,
            merqury_qv_values,
            merqury_completeness_values,
            haploid_number,
            haploid_source,
            taxon_id,
            family,
            resolution_level,
            args.genomescope_plot,
            merqury_plots,
            args.hic,
            args.blob,
            args.fcs_gx,
            inspector_values,
            args.output,
            kmer_tool=kmer_tool,
            pairtools_list=pairtools_list,
            chromap_list=chromap_list,
            expected_chr_size=expected_chr_size,
            merqury_both=merqury_both,
            kmer_false_dup=kmer_false_dup,
            hic_library_qc=hic_library_qc,
            version_lines=version_lines,
            hic_library_pairs=hic_library_pairs
        )

        # Convert to PDF if requested
        if args.also_pdf:
            convert_md_to_pdf(args.output)

    except FileNotFoundError as e:
        print(f"Error: File not found - {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
