#!/usr/bin/env python3

# GEP2 Genome Stats Report Generator
# by Diego De Panis, 2025
# This script is part of the GEP2 pipeline

# Generates a markdown report aggregating results from:
# - gfastats (assembly metrics)
# - compleasm (gene completeness)
# - Merqury (k-mer QV and completeness)
# - GenomeScope2 (genome profiling)

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


def convert_md_to_pdf(md_file, pdf_file=None):
    """
    Convert markdown file to PDF using pandoc with weasyprint.
    
    Requirements:
        pip install weasyprint --break-system-packages
        apt install pandoc (or conda install pandoc)
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
    """Get haploid number and source from GoaT API based on species name."""
    try:
        species_encoded = quote(species)
        search_url = f'https://goat.genomehubs.org/api/v2/search?query=tax_name%28{species_encoded}%29&result=taxon'
        response = requests.get(search_url, timeout=30)
        response.raise_for_status()
        
        search_data = response.json()
        
        if not search_data.get('results'):
            return None, None, "Species not found"
        
        first_result = search_data['results'][0]['result']
        taxon_id = first_result['taxon_id']
        
        record_url = f'https://goat.genomehubs.org/api/v2/record?recordId={taxon_id}&result=taxon&taxonomy=ncbi'
        record_response = requests.get(record_url, timeout=30)
        record_response.raise_for_status()
        
        record_data = record_response.json()
        
        if not record_data.get('records'):
            return None, None, "No genomic records found"
        
        attributes = record_data['records'][0]['record'].get('attributes', {})
        haploid_info = attributes.get('haploid_number', {})
        haploid_number = haploid_info.get('value')
        haploid_source = haploid_info.get('aggregation_source')
        
        if haploid_number is None:
            return None, None, "Haploid number not available"
        
        return int(haploid_number), haploid_source, None
        
    except requests.exceptions.RequestException as e:
        return None, None, f"API request failed: {e}"
    except (KeyError, ValueError, TypeError) as e:
        return None, None, f"Data parsing error: {e}"
    except Exception as e:
        return None, None, f"Unexpected error: {e}"


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


def parse_merqury_qv(filepath, num_assemblies):
    """
    Parse Merqury QV file.
    Format: name\tErrors\tBases\tQV\tErrorRate
    Returns list of QV values for each assembly.
    """
    qv_values = []
    
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        for i in range(min(num_assemblies, len(lines))):
            parts = lines[i].strip().split('\t')
            if len(parts) >= 4:
                qv_values.append(float(parts[3]))
            else:
                qv_values.append(None)
        
        # Pad with None if fewer lines than expected
        while len(qv_values) < num_assemblies:
            qv_values.append(None)
            
    except Exception as e:
        print(f"Warning: Could not parse Merqury QV file: {e}")
        qv_values = [None] * num_assemblies
    
    return qv_values


def parse_merqury_completeness(filepath, num_assemblies):
    """
    Parse Merqury completeness file.
    Format: name\tall\tFound\tTotal\tCompleteness
    Returns list of completeness values for each assembly.
    """
    completeness_values = []
    
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        for i in range(min(num_assemblies, len(lines))):
            parts = lines[i].strip().split('\t')
            if len(parts) >= 5:
                completeness_values.append(float(parts[4]))
            else:
                completeness_values.append(None)
        
        while len(completeness_values) < num_assemblies:
            completeness_values.append(None)
            
    except Exception as e:
        print(f"Warning: Could not parse Merqury completeness file: {e}")
        completeness_values = [None] * num_assemblies
    
    return completeness_values


def find_merqury_plots(merqury_dir, asm_id):
    """
    Find Merqury plot files in the output directory.
    
    Haploid (1 asm) produces:
      - {prefix}.spectra-cn.fl.png
      - {prefix}.spectra-asm.fl.png
    
    Diploid (2 asm) produces:
      - {prefix}.{asm1_name}.spectra-cn.fl.png
      - {prefix}.{asm2_name}.spectra-cn.fl.png
      - {prefix}.spectra-asm.fl.png
      - {prefix}.spectra-cn.fl.png (combined)
    """
    plots = {
        'spectra_cn': [],
        'spectra_asm': None,
        'spectra_cn_combined': None
    }
    
    if not merqury_dir or not os.path.isdir(merqury_dir):
        return plots
    
    # Find ALL spectra-cn plots
    cn_pattern = os.path.join(merqury_dir, "*.spectra-cn.fl.png")
    cn_files = sorted(glob.glob(cn_pattern))
    
    for cn_file in cn_files:
        basename = os.path.basename(cn_file)
        
        # FIX: Use Regex instead of dot counting.
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
    asm_pattern = os.path.join(merqury_dir, "*.spectra-asm.fl.png")
    asm_files = glob.glob(asm_pattern)
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


def generate_report(species_name, assembly_id, gfastats_list, compleasm_list, 
                   merqury_qv_values, merqury_completeness_values,
                   haploid_number, haploid_source, genomescope_plot, 
                   merqury_plots, hic_files, output_file):
    """Generate the markdown report supporting 1 or 2 assemblies."""
    
    num_assemblies = len(gfastats_list)
    is_diploid = num_assemblies == 2
    
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
    
    # ---- Merqury metrics ----
    
    if merqury_qv_values and any(v is not None for v in merqury_qv_values):
        ratings = [get_rating(v, 'merqury_qv') for v in merqury_qv_values]
        add_metric("MERQ QV", merqury_qv_values, ratings)
    
    if merqury_completeness_values and any(v is not None for v in merqury_completeness_values):
        ratings = [get_rating(v, 'merqury_completeness') for v in merqury_completeness_values]
        add_metric("MERQ Compl.", merqury_completeness_values, ratings)
    
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
        haploid_info = f"‡ = Haploid number is {haploid_number} ({haploid_source}, [GoaT](https://goat.genomehubs.org))<br>"
    else:
        haploid_info = "‡ = Haploid number not found on [GoaT](https://goat.genomehubs.org)<br>"
    
    # Build the complete report
    report_lines = [
        "GEP2 genome stats report",
        "---",
        f"### {species_name}",
        f"#### {assembly_id}",
        "",
        header,
        separator,
    ]
    report_lines.extend(table_rows)
    
    # Add legend only for metrics that were included
    report_lines.append("")
    
    if has_compleasm:
        report_lines.append("e = eukaryota (odb12)<br>")
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
                "### Merqury Plots",
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
    if hic_files:
        # Filter to ensure files actually exist
        valid_hic_files = [f for f in hic_files if os.path.exists(f)]
        
        if valid_hic_files:
            report_dir = os.path.dirname(os.path.abspath(output_file))
            
            report_lines.extend([
                "",
                "---",
                "### Hi-C Contact Maps",
            ])
            
            for i, hic_plot in enumerate(valid_hic_files):
                rel_path = os.path.relpath(hic_plot, report_dir)
                
                # Determine label: use "asmX" if multiple files, otherwise generic "Assembly"
                if len(valid_hic_files) > 1:
                    label = f"asm{i+1}"
                else:
                    label = "Assembly"
                
                report_lines.extend([
                    "",
                    f"#### {label}",
                    f"![Hi-C Map - {label}]({rel_path})",
                ])



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
    
    parser.add_argument('-s', '--species', required=True, 
                       help='Species name (with underscore, e.g., Homo_sapiens)')
    parser.add_argument('-a', '--assembly', required=True, 
                       help='Assembly ID')
    parser.add_argument('-g', '--gfastats', required=True, nargs='+',
                       help='gfastats output file(s) - one per assembly')
    parser.add_argument('-c', '--compleasm', required=False, nargs='+', default=[],
                       help='compleasm output file(s) - one per assembly')
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
        if args.compleasm and len(args.compleasm) != num_assemblies:
            print(f"Error: Number of compleasm files ({len(args.compleasm)}) must match "
                  f"number of gfastats files ({num_assemblies})", file=sys.stderr)
            sys.exit(1)
        
        # Clean species name
        species_clean = args.species.replace('_', ' ')
        print(f"Processing species: {species_clean}")
        print(f"Number of assemblies: {num_assemblies}")
        
        # Get haploid number from GoaT
        print("Fetching haploid number from GoaT...")
        haploid_number, haploid_source, error = get_species_genomic_data_from_goat(species_clean)
        
        if error:
            print(f"Warning: Could not retrieve haploid number from GoaT: {error}")
            haploid_number, haploid_source = None, None
        else:
            print(f"Retrieved haploid number: {haploid_number} (source: {haploid_source})")
        
        # Parse gfastats files
        print("Parsing gfastats output(s)...")
        gfastats_list = []
        for gf_file in args.gfastats:
            gfastats_list.append(parse_gfastats(gf_file))
        
        # Parse compleasm files (if provided)
        print("Parsing compleasm output(s)...")
        compleasm_list = []
        if args.compleasm:
            for comp_file in args.compleasm:
                compleasm_list.append(parse_compleasm(comp_file))
        else:
            # Empty compleasm results
            compleasm_list = [{'eukaryota_single': None, 'eukaryota_dupl': None,
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
        
        # Create output directory if needed
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        
        # Generate report
        print("Generating report...")
        generate_report(
            species_clean, 
            args.assembly, 
            gfastats_list, 
            compleasm_list,
            merqury_qv_values,
            merqury_completeness_values,
            haploid_number, 
            haploid_source, 
            args.genomescope_plot,
            merqury_plots,
            args.hic,
            args.output
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
