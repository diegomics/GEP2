#!/usr/bin/env python3

# GEP2 Genome Stats Report Generator
# by Diego De Panis, 2025
# This script is part of the GEP2 pipeline
# note: AI tools may have been used to improve, clean and/or comment this version of the code

# Generates a markdown report aggregating results from:
# - gfastats (assembly metrics)
# - compleasm/busco (gene completeness)
# - busco (gene completeness)
# - Merqury/MerquryFK (k-mer QV and completeness)
# - GenomeScope2 (genome profiling)
# - Inspector (structural errors)
# - Blobtools (contamination screening)
# - FCS-GX (contamination screening)

import argparse
import re
import sys
import os
import glob
import subprocess
import shutil
import requests
from urllib.parse import quote
from pathlib import Path

__version__ = '0.1.9'

# This is a crap and isn't working yet, will work on it soon...
def convert_md_to_pdf(md_file, pdf_file=None):
    """
    Convert markdown file to PDF using pandoc with weasyprint.
    
    Requirements:
        weasyprint --break-system-packages
        pandoc
    """
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
    
    for key, pattern in patterns.items():
        match = re.search(pattern, content)
        if match:
            metrics[key] = match.group(1)
    
    # Calculate gaps per Gbp
    if 'total_bp' in metrics and 'gaps' in metrics:
        total_bp = int(metrics['total_bp'])
        gaps = int(metrics['gaps'])
        gaps_per_gbp = (gaps / total_bp) * 1_000_000_000 if total_bp > 0 else 0
        metrics['gaps_per_gbp'] = gaps_per_gbp
    
    return metrics


def parse_compleasm(filepath):
    """Parse compleasm output file and extract completeness metrics."""
    metrics = {
        'eukaryota_single': None,
        'eukaryota_dupl': None,
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
        
        if s_match and d_match:
            s_percent = float(s_match.group(1))
            d_percent = float(d_match.group(1))
            
            if 'eukaryota' in lineage_name.lower():
                metrics['eukaryota_single'] = s_percent
                metrics['eukaryota_dupl'] = d_percent
            else:
                other_name = lineage_name.replace('_', ' (') + ')'
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
              'Fragmented': 0, 'Interspersed': 0, 'Missing': 0}
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
    d_percent = (counts['Duplicated'] / total) * 100
    
    # Frameshift rate
    if total_complete_copies > 0:
        frameshift_rate = (frameshifted_copies / total_complete_copies) * 100
    else:
        frameshift_rate = 0.0
    
    is_eukaryota = 'eukaryota' in lineage.lower()
    
    if is_eukaryota:
        metrics['eukaryota_single'] = s_percent
        metrics['eukaryota_dupl'] = d_percent
    else:
        lineage_name = lineage.replace('_odb', ' (odb') + ')'
        metrics['other_lineage'] = lineage_name
        metrics['other_single'] = s_percent
        metrics['other_dupl'] = d_percent
    
    metrics['frameshift_rate'] = frameshift_rate
    metrics['frameshift_is_eukaryota'] = is_eukaryota
    
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

        for parts in rows[start:]:
            if len(qv_values) >= num_assemblies:
                break
            label = parts[0].strip().lower() if parts else ""
            if label == "both" or "+" in label:      # skip combined row
                continue

            if len(parts) > qv_idx:
                try:
                    qv_values.append(float(parts[qv_idx]))
                    continue

                except ValueError:
                    pass

        # Pad with None if fewer lines than expected
        while len(qv_values) < num_assemblies:
            qv_values.append(None)

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

        for parts in rows[start:]:
            if len(completeness_values) >= num_assemblies:
                break
            label = parts[0].strip().lower() if parts else ""
            if label == "both" or "+" in label:      # skip asm1+asm2 row
                continue
            if len(parts) > pct_idx:
                try:
                    completeness_values.append(float(parts[pct_idx]))
                    continue
                except ValueError:
                    pass

        while len(completeness_values) < num_assemblies:
            completeness_values.append(None)

    except Exception as e:
        print(f"Warning: Could not parse completeness file: {e}")
        completeness_values = [None] * num_assemblies

    return completeness_values


def find_merqury_plots(merqury_dir, asm_id):
    """Find Merqury/MerquryFK spectra plots. Both use spectra-cn/spectra-asm with
    .fl/.ln/.st suffixes; prefer filled (.fl), fall back to line/stack."""
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
    
    # Find ALL spectra-cn plots
    cn_files = _first_glob("*.spectra-cn.fl.png", "*.spectra-cn.ln.png", "*.spectra-cn.st.png")

    for cn_file in cn_files:
        basename = os.path.basename(cn_file)

        # Use Regex instead of dot counting.
        # This correctly handles prefixes with dots (e.g., gfArmOsto1.1)
        # It looks for .asmX. or .hapX. or _asmX_ tags.   
        if re.search(r'[\._](asm|hap)\d+[\._]', basename):
            # It has a tag -> Per-assembly plot
            plots['spectra_cn'].append(cn_file)
        else:
            # No tag -> Combined plot
            plots['spectra_cn_combined'] = cn_file

    # For haploid, the only spectra-cn is the "combined" one, move it to spectra_cn list
    # (Because in haploid mode, the 'combined' plot IS the assembly plot)
    if not plots['spectra_cn'] and plots['spectra_cn_combined']:
        plots['spectra_cn'].append(plots['spectra_cn_combined'])
        plots['spectra_cn_combined'] = None

    # Find spectra-asm plot
    asm_files = _first_glob("*.spectra-asm.fl.png", "*.spectra-asm.ln.png", "*.spectra-asm.st.png")
    if asm_files:
        plots['spectra_asm'] = asm_files[0]

    return plots


def get_rating(value, metric_type, haploid_number=None):
    """Get star rating based on metric type and value."""
    if value is None:
        return '····'
    
    if metric_type == 'gaps_per_gbp':
        if value < 200:
            return '****'
        elif value < 1000:
            return '***-'
        elif value < 10000:
            return '**--'
        else:
            return '*---'
    
    elif metric_type == 'scaffold_n50':
        if value > 100_000_000:
            return '****'
        elif value > 10_000_000:
            return '***-'
        elif value > 100_000:
            return '**--'
        else:
            return '*---'
    
    elif metric_type == 'contig_n50':
        if value > 10_000_000:
            return '****'
        elif value > 1_000_000:
            return '***-'
        elif value > 100_000:
            return '**--'
        else:
            return '*---'
    
    elif metric_type == 'compl_single':
        if value > 95:
            return '****'
        elif value > 90:
            return '***-'
        elif value > 80:
            return '**--'
        else:
            return '*---'
    
    elif metric_type == 'compl_dupl':
        if value < 2:
            return '****'
        elif value < 5:
            return '***-'
        elif value < 7:
            return '**--'
        else:
            return '*---'
    
    elif metric_type == 'l90_haploid':
        if haploid_number is None:
            return '····'
        if value <= haploid_number:
            return '****'
        elif value <= haploid_number + 100:
            return '***-'
        elif value <= haploid_number + 1000:
            return '**--'
        else:
            return '*---'
    
    elif metric_type == 'merqury_qv':
        if value > 50:
            return '****'
        elif value > 40:
            return '***-'
        elif value > 30:
            return '**--'
        else:
            return '*---'
    
    elif metric_type == 'merqury_completeness':
        if value > 95:
            return '****'
        elif value > 90:
            return '***-'
        elif value > 80:
            return '**--'
        else:
            return '*---'
    
    elif metric_type == 'compl_frameshift':
        if value < 2:
            return '****'
        elif value < 5:
            return '***-'
        elif value < 15:
            return '**--'
        else:
            return '*---'
    
    return '····'


def format_number(value):
    """Format large numbers with commas."""
    if value is None:
        return "N/A"
    if isinstance(value, (int, float)):
        if isinstance(value, float):
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


def generate_report(species_name, assembly_id, gfastats_list, compleasm_list, busco_list, 
                   merqury_qv_values, merqury_completeness_values,
                   haploid_number, haploid_source, taxon_id, family, resolution_level,
                   genomescope_plot, 
                   merqury_plots, hic_plots, blob_plots, fcs_gx_files,
                   inspector_values, output_file, kmer_tool="MERQ"):
    """Generate the markdown report supporting 1 or 2 assemblies."""
    
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
    ratings = [get_rating(v, 'scaffold_n50') for v in values]
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
    ratings = [get_rating(v, 'contig_n50') for v in values]
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
        ratings = [get_rating(v, 'merqury_completeness') for v in merqury_completeness_values]
        add_metric(f"{kmer_tool} Compl.", merqury_completeness_values, ratings)
    
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

    # Grab dynamic names
    euka_lineage = next((m.get('eukaryota_lineage') for m in busco_list if m.get('eukaryota_lineage')), "eukaryota (odb12)")
    
    other_lineage = next((m.get('other_lineage') for m in compleasm_list if m.get('other_lineage')), None)
    if not other_lineage:
        other_lineage = next((m.get('other_lineage') for m in busco_list if m.get('other_lineage')), None)

    # Print legend ONLY for the lineages that were actually used
    if has_compleasm or has_busco:
        if has_eukaryota:
            report_lines.append(f"e = {euka_lineage}<br>")
        if other_lineage:
            report_lines.append(f"x = {other_lineage}<br>")
    
    report_lines.append(haploid_info)
    
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


    # Footer
    report_lines.extend([
        "",
    ])
    
    # Write to output file
    report = "\n".join(report_lines)
    
    with open(output_file, 'w') as f:
        f.write(report)
    
    print(f"Report generated successfully: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate GEP2 genome assembly stats report',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single assembly (haploid mode)
  python make_report.py -s Homo_sapiens -a hg38 \\
    -g stats1.txt -c compl1.txt \\
    -q merqury.qv -m merqury.completeness.stats \\
    --genomescope-plot linear_plot.png \\
    --merqury-dir kmer_stats/ \\
    -o report.md --also-pdf

  # Two assemblies (diploid mode)
  python make_report.py -s Homo_sapiens -a hg38 \\
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
        
        # Find Merqury plots
        merqury_plots = find_merqury_plots(args.merqury_dir, args.assembly)

        # Identify the k-mer QV tool from the file header (Merqury vs MerquryFK)
        kmer_tool = detect_kmer_tool(args.merqury_qv, args.merqury_completeness)
        
        # Parse Inspector files (if provided)
        inspector_values = [None] * num_assemblies
        if args.Inspector:
            print("Parsing Inspector output(s)...")
            for i, insp_file in enumerate(args.Inspector[:num_assemblies]):
                inspector_values[i] = parse_inspector(insp_file)
        
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
            kmer_tool=kmer_tool
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
    