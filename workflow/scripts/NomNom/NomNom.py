#!/usr/bin/env python3
r"""
(∪• ﻌ •∪)
"""
# by Diego De Panis, 2025
# This script is part of the GAME-pipeline
# See https://github.com/diegomics/GAME/tree/main/scripts/NomNom

import pandas as pd
import numpy as np
import re
import os
import sys
import yaml
from typing import Dict, List, Any, Optional
import openpyxl
import csv
import argparse
from numbers_parser import Document
from odf.opendocument import load
from odf.table import Table, TableRow, TableCell
from odf.text import P
import xlrd

import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

def debug_wrapper(func):
    """Decorator to catch and report errors in functions"""
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            print(f"ERROR in {func.__name__}: {str(e)}")
            import traceback
            traceback.print_exc()
            if "bool" in str(e) and "not iterable" in str(e):
                # Find which argument might be a boolean
                for i, arg in enumerate(args):
                    if isinstance(arg, bool):
                        print(f"Argument {i} is bool: {arg}")
                for k, v in kwargs.items():
                    if isinstance(v, bool):
                        print(f"Keyword argument {k} is bool: {v}")
            return None
    return wrapper

class TableAnalyzer:
    """
    A comprehensive tool for analyzing tabular data with specific focus on
    genomic sequencing data tables. Features include:
    - Automatic delimiter detection
    - Data validation and error checking
    - Paired-end read detection for specific read types
    - Hierarchical YAML output generation
    """
    
    def __init__(self):
        self.analysis_results = {"file_info": {}, "warnings": [], "errors": []}
        self.nested_categories = {}
        self.column_data_types = {}
        self.column_values = {}
        self.warnings = []
        self.errors = []
        self.structured_data = {}
        self.current_df = None  # Store the current DataFrame for later access

    def _detect_delimiter(self, file_path: str) -> str:
        """
        Detect the delimiter in a file by examining the first few lines.
        Supports comma, semicolon, tab, pipe, and vertical bar delimiters.
        """
        delimiters = {',': 0, ';': 0, '\t': 0, '|': 0, '\\': 0}
        
        with open(file_path, 'r', errors='ignore') as file:
            # Check first 5 lines or fewer if file is smaller
            lines = []
            for _ in range(5):
                line = file.readline()
                if not line:
                    break
                lines.append(line)
            
            # For each line, count occurrences
            for line in lines:
                # Count tabs more accurately by checking actual tab characters
                delimiters['\t'] += line.count('\t')
                
                # Count other delimiters
                for delimiter in [',', ';', '|', '\\']:
                    delimiters[delimiter] += line.count(delimiter)
        
        # Special check for tab-delimited files that might have spaces
        # If we have consistent tab counts across lines, it's likely tab-delimited
        if lines and delimiters['\t'] > 0:
            tab_counts = [line.count('\t') for line in lines if line.strip()]
            if tab_counts and all(count == tab_counts[0] for count in tab_counts):
                return '\t'
        
        # Return the most common delimiter, or default to comma if none found
        most_common = max(delimiters.items(), key=lambda x: x[1])
        return most_common[0] if most_common[1] > 0 else ','
    
    def _is_boolean(self, value) -> bool:
        """Check if a value represents a boolean"""
        if value is None or pd.isna(value) or value == '':
            return False
            
        # If it's already a boolean type, return True
        if isinstance(value, bool):
            return True
            
        # Convert to string and check against boolean values
        try:
            str_value = str(value).strip().lower()
            boolean_values = {'true', 'false', 'yes', 'no', 'y', 'n', 't', 'f', '1', '0'}
            return str_value in boolean_values
        except:
            return False

    def _to_boolean(self, value) -> bool:
        """Convert a string value to a boolean"""
        if value is None or pd.isna(value) or value == '':
            return False
            
        # If it's already a boolean type, return it directly
        if isinstance(value, bool):
            return value
            
        # Convert to string and check against true values
        try:
            str_value = str(value).strip().lower()
            true_values = {'true', 'yes', 'y', 't', '1'}
            return str_value in true_values
        except:
            return False
    
    def _is_integer(self, value) -> bool:
        """Check if a value represents an integer"""
        if value is None or pd.isna(value) or value == '':
            return False
            
        str_value = str(value).strip()
        try:
            int(str_value)
            return True
        except:
            return False


    def _infer_column_types(self, df: pd.DataFrame) -> Dict[str, str]:
        """
        Infer column types based on content analysis, excluding special placeholder values.
        """
        inferred_types = {}
        
        # Define special values that should be excluded from type inference calculations
        special_values = {'', 'na', 'nan', 'none', 'null', 'auto', 'missing', '-', ' '}
        
        # Define read file extensions to identify path columns
        read_extensions = ['.fastq', '.fq', '.fasta', '.fa', '.fna', '.bam', '.vcf', '.bcf', '.gz']
        
        # Columns that often contain paired reads
        likely_paired_columns = ['hic_reads', 'short_reads', 'illumina_reads', '10x_reads']

        # Columns that are likely paths - specifically include asm_file
        path_column_indicators = ['file', 'files', 'path', 'paths', 'reads', 'assembly', 'asm_file', 'read_files', 'reference', 'ref']
        
        # Force these columns to always be treated as simple paths
        always_simple_path_columns = ['asm_file', 'assembly_file', 'reference_file', 'ref_file', 'genome_file']
        
        for column in df.columns:
            # First check: if column name is in our always_simple_path list, immediately set type
            if column.lower() in always_simple_path_columns or any(indicator in column.lower() for indicator in always_simple_path_columns):
                inferred_types[column] = "simple_path"
                continue
                
            # Initialize counters for different value types
            bool_values = 0
            int_values = 0
            path_values = []
            paired_patterns = []
            multi_paths = []
            categorical_values = set()
            non_special_count = 0  # Count of values that aren't special placeholders
            
            # NEW: Track non-boolean values that disqualify boolean classification
            non_boolean_values = []
            
            # Get column name in lowercase for comparison
            col_lower = column.lower()
            
            # Special case for known column types based on name - as a last resort
            # after content-based analysis doesn't yield clear results
            is_likely_integer = col_lower in ['genome_size', 'kmer_size', 'kmer', 'size', 'length']
            is_likely_boolean = col_lower in ['trim', 'read_qc', 'qc', 'flag']
            is_likely_paired = col_lower in likely_paired_columns
            
            if hasattr(self, 'merged_paired_columns') and column in self.merged_paired_columns:
                # Force it to be complex_path since it contains paired data
                inferred_types[column] = "complex_path"
                continue

            # Analyze each value in the column
            for value in df[column]:
                # Skip None/NaN
                if pd.isna(value):
                    continue
                    
                # Convert to string for easier comparison
                str_value = str(value).strip().lower() if not isinstance(value, bool) else str(value).lower()
                
                # Check if this is a special placeholder value
                if str_value in special_values:
                    continue
                    
                # Increment non-special value counter
                non_special_count += 1
                
                # Track unique values for categorical analysis
                categorical_values.add(str_value)
                
                # Check for boolean values - but exclude numeric strings like "0", "1200000"
                if self._is_boolean(value) and not (str_value.isdigit() and len(str_value) > 1):
                    bool_values += 1
                else:
                    # NEW: Track values that are clearly NOT boolean
                    non_boolean_values.append(str_value)
                    
                # Check for integer values - improved to catch large numbers
                if self._is_integer(value):
                    int_values += 1
                    
                # Skip non-string values for path analysis
                if not isinstance(value, str):
                    continue
                    
                # Check if it's a path
                contains_path_separator = '/' in value or '\\' in value
                
                # CRITICAL FIX: Check if it contains comma-separated paths
                # This is the key check for complex paths
                has_multiple_paths = ',' in value and (value.count('/') > 1 or value.count('\\') > 1)
                
                # Check for read file extensions
                has_read_extension = any(ext in value.lower() for ext in read_extensions)
                
                # Check for paired read patterns
                has_paired_pattern = any(pattern in value for pattern in 
                                        ['_R1', '_R2', '_1.', '_2.', '.1.', '.2.'])
                
                if contains_path_separator:
                    path_values.append(value)
                        
                if has_paired_pattern:
                    paired_patterns.append(value)
                    
                if has_multiple_paths:
                    multi_paths.append(value)
            
            # If we have no non-special values, use column name as a fallback
            if non_special_count == 0:
                if is_likely_integer:
                    inferred_types[column] = "integer"
                elif is_likely_boolean:
                    inferred_types[column] = "boolean"
                elif is_likely_paired:
                    inferred_types[column] = "complex_path"
                elif any(indicator in col_lower for indicator in path_column_indicators):
                    inferred_types[column] = "simple_path"  # Default to simple_path when empty
                else:
                    inferred_types[column] = "string"
            
            # Now make a decision based on the actual content (excluding special values)
            # 1. Integer type - check first to prevent large numbers from being classified as boolean
            elif int_values / non_special_count >= 0.7 or is_likely_integer:  # Relaxed threshold and added name check
                inferred_types[column] = "integer"
            
            # 2. Path detection - check BEFORE boolean to prevent file paths from being treated as boolean
            elif len(path_values) > 0:
                # CRITICAL: Check for multiple paths or paired patterns
                if multi_paths or len(paired_patterns) > 0:
                    inferred_types[column] = "complex_path"
                else:
                    inferred_types[column] = "simple_path"
            
            # 3. Boolean type - ENHANCED: Only if ALL non-special values can be boolean
            # AND we don't have obvious non-boolean content like file paths
            elif (bool_values / non_special_count >= 0.5 and 
                not is_likely_integer and 
                len(non_boolean_values) == 0 and  # NEW: No clearly non-boolean values
                len(path_values) == 0):  # NEW: No file paths
                inferred_types[column] = "boolean"
                    
            # 4. Complex Path - if the column contains comma-separated paths or paired read patterns
            elif multi_paths or (len(paired_patterns) > 0 and len(path_values) > 0):
                inferred_types[column] = "complex_path"
                    
            # 5. Categorical type - if there are few unique values relative to non-special values
            elif len(categorical_values) <= min(10, non_special_count * 0.3):
                inferred_types[column] = "category"

            # 6. Path column by name - if the column name suggests it contains paths
            elif any(indicator in col_lower for indicator in path_column_indicators):
                # IMPORTANT: Check if we have evidence of multiple paths
                if multi_paths or len(paired_patterns) > 0:
                    inferred_types[column] = "complex_path"
                else:
                    inferred_types[column] = "simple_path"
                            
            # Default to string for all other columns
            else:
                inferred_types[column] = "string"
        
        return inferred_types


    def _infer_paired_data_columns(self, df: pd.DataFrame) -> List[str]:
        """
        Identify columns that likely contain paired-end sequencing data.
        
        Args:
            df: DataFrame to analyze
            
        Returns:
            List of column names that likely contain paired-end data
        """
        paired_columns = []
        
        for column in df.columns:
            # Get non-null values for analysis
            values = [str(v) for v in df[column].dropna() if pd.notna(v)]
            
            # Skip empty columns
            if not values:
                continue
                
            # Check for indicators of paired data
            paired_indicators = [
                # R1/R2 file pattern
                sum('_R1' in v and '_R2' in v for v in values),
                # _1/_2 file pattern
                sum('_1.' in v and '_2.' in v for v in values),
                # Paired file paths separated by comma
                sum(v.count(',') == 1 for v in values),
                # Even number of file paths per entry
                sum(v.count(',') % 2 == 1 for v in values),
            ]
            
            # If multiple paired indicators are present, it's likely paired data
            if sum(bool(i) for i in paired_indicators) >= 2:
                paired_columns.append(column)
                continue
                
            # Check for paired data based on column name as a fallback
            col_lower = column.lower()
            if any(term in col_lower for term in ['paired', 'pe', 'hi-c', 'hic', '10x', 'short']):
                paired_columns.append(column)
                
        return paired_columns

    def _fix_malformed_headers(self, lines):
        """
        Fix headers that have multiple spaces acting as delimiters.
        Common issue: "HiC_R2  Long_Reads" should be "HiC_R2\tLong_Reads"
        """
        if not lines:
            return lines
        
        header = lines[0]
        
        # Check if we have multiple consecutive spaces that might be intended as delimiters
        if '  ' in header:  # Two or more spaces
            # Count different delimiter candidates
            tab_count = header.count('\t')
            multi_space_matches = re.findall(r' {2,}', header)
            multi_space_count = len(multi_space_matches)
            
            # If we have multiple spaces and few/no tabs, spaces are likely the delimiter
            if multi_space_count > 0:
                # Get the most common multi-space length
                space_lengths = [len(match) for match in multi_space_matches]
                
                # Replace multiple spaces with tabs
                fixed_header = re.sub(r' {2,}', '\t', header)
                
                # Verify this creates a reasonable number of columns
                potential_cols = fixed_header.split('\t')
                if len(potential_cols) > 1:
                    self.warnings.append(f"Fixed header with {multi_space_count} multiple-space delimiters")
                    
                    # Create new lines list with fixed header
                    fixed_lines = [fixed_header]
                    
                    # Fix data lines too
                    for i in range(1, len(lines)):
                        if '  ' in lines[i]:
                            fixed_line = re.sub(r' {2,}', '\t', lines[i])
                            fixed_lines.append(fixed_line)
                        else:
                            fixed_lines.append(lines[i])
                    
                    return fixed_lines
        
        return lines

    def _detect_read_type_column(self, df: pd.DataFrame) -> Optional[str]:
        """
        Identify the column that contains read type information.
        
        Args:
            df: DataFrame to analyze
            
        Returns:
            Name of column containing read types, or None if not found
        """
        # Check for columns with common read type values
        read_type_values = {
            'hic', 'hi-c', 'hifi', 'clr', 'ont', 'nanopore', 'illumina', 
            'pacbio', '10x', 'short', 'long', 'paired', 'single'
        }
        
        for column in df.columns:
            # Get lowercase string values
            values = [str(v).lower() for v in df[column].dropna() if pd.notna(v)]
            
            # Skip empty columns
            if not values:
                continue
                
            # Count how many values match known read types
            matches = sum(any(rt in v for rt in read_type_values) for v in values)
            
            # If more than 50% of values match known read types, it's likely the read type column
            if matches / len(values) >= 0.5:
                return column
                
        # Check column names as a fallback
        for column in df.columns:
            col_lower = column.lower()
            if 'type' in col_lower or 'platform' in col_lower or 'technology' in col_lower:
                return column
                
        # If no clear read type column is found, return None
        return None

    def _merge_paired_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Detect and merge paired columns (like Library_R1/Library_R2, HiC_R1/HiC_R2)
        into a single column with paired values using _x suffix.
        
        Examples:
        - Library_R1/Library_R2 → Library_Rx
        - HiC_R1/HiC_R2 → HiC_Rx
        - Read_1/Read_2 → Read_x
        - IlluminaF/IlluminaR → Illuminax
        - Forward/Reverse → Paired_x (special case)
        - R1/R2 → Paired_x (special case)
        """
        # Common patterns for paired columns
        paired_patterns = [
            ('_R1', '_R2'),
            ('_1', '_2'),
            ('_F', '_R'),
            ('_fwd', '_rev'),
            ('_forward', '_reverse'),
            ('_Forward', '_Reverse'),
            ('_LEFT', '_RIGHT'),
            ('_left', '_right'),
            ('_L', '_R'),
            ('F', 'R'),  # For cases like IlluminaF/IlluminaR
            ('1', '2'),  # For cases like Read1/Read2
            ('_r1', '_r2'),
            ('.R1', '.R2'),
            ('.1', '.2'),
            ('_1.', '_2.'),
            ('.r1', '.r2'),
            ('_READ1', '_READ2'),
            ('_read1', '_read2'),
        ]
        
        # Special standalone column names that need custom handling
        special_standalone_pairs = [
            ('Forward', 'Reverse'),
            ('forward', 'reverse'),
            ('R1', 'R2'),
            ('r1', 'r2'),
            ('Left', 'Right'),
            ('left', 'right'),
            ('F', 'R'),
            ('1', '2'),
        ]
        
        merged_df = df.copy()
        columns_to_drop = []
        processed_columns = set()  # Track which columns we've already processed
        
        # First, check for special standalone pairs
        for col1_name, col2_name in special_standalone_pairs:
            if col1_name in df.columns and col2_name in df.columns:
                if col1_name not in processed_columns and col2_name not in processed_columns:
                    # Create new column name for special cases
                    new_col_name = "Paired_x"
                    
                    # If Paired_x already exists, add a number
                    if new_col_name in merged_df.columns:
                        counter = 2
                        while f"Paired{counter}_x" in merged_df.columns:
                            counter += 1
                        new_col_name = f"Paired{counter}_x"
                    
                    # Merge the columns
                    merged_values = []
                    for idx, row in df.iterrows():
                        val1 = row[col1_name]
                        val2 = row[col2_name]
                        
                        # Handle None/NaN values
                        if pd.isna(val1) or val1 == 'None':
                            val1 = 'None'
                        if pd.isna(val2) or val2 == 'None':
                            val2 = 'None'
                            
                        # If both are None, keep as None
                        if val1 == 'None' and val2 == 'None':
                            merged_values.append('None')
                        else:
                            # Create paired value
                            merged_values.append(f"{val1}, {val2}")
                    
                    # Add merged column
                    merged_df[new_col_name] = merged_values
                    
                    # Mark original columns for dropping
                    columns_to_drop.extend([col1_name, col2_name])
                    processed_columns.update([col1_name, col2_name])
                    
                    # Mark this as a paired column
                    if not hasattr(self, 'merged_paired_columns'):
                        self.merged_paired_columns = {}
                    self.merged_paired_columns[new_col_name] = True
        
        # Now check for pattern-based pairs
        for col in df.columns:
            # Skip if already processed
            if col in columns_to_drop or col in processed_columns:
                continue
                
            # Check each pattern
            for suffix1, suffix2 in paired_patterns:
                # For single character suffixes (like F/R), ensure they're at the end
                if len(suffix1) == 1:
                    if col.endswith(suffix1) and len(col) > 1:
                        # Find the base name
                        base_name = col[:-len(suffix1)]
                        paired_col = base_name + suffix2
                    else:
                        continue
                else:
                    if suffix1 in col:
                        # Find where the suffix appears
                        suffix_pos = col.rfind(suffix1)
                        if suffix_pos == -1:
                            continue
                        
                        # Get base name and construct paired column name
                        base_name = col[:suffix_pos]
                        remainder = col[suffix_pos + len(suffix1):]
                        paired_col = base_name + suffix2 + remainder
                    else:
                        continue
                
                # Check if paired column exists
                if paired_col in df.columns and paired_col not in processed_columns:
                    # Clean up the base name for the new column
                    # Remove trailing underscores or dots
                    clean_base = base_name.rstrip('_.')
                    
                    # Create new column name with appropriate suffix
                    # Strategy: preserve 'R' in R1/R2 patterns only, use 'x' for everything else
                    
                    # For R1/R2 patterns, preserve the 'R'
                    if suffix1 in ['_R1', '_r1', '.R1', '.r1', '_READ1', '_read1'] and \
                    suffix2 in ['_R2', '_r2', '.R2', '.r2', '_READ2', '_read2']:
                        # Determine separator and case
                        separator = '.' if suffix1.startswith('.') else '_'
                        use_upper = 'R' in suffix1  # Preserve case
                        new_suffix = f"{separator}{'Rx' if use_upper else 'rx'}"
                        new_col_name = f"{clean_base}{new_suffix}" if clean_base else f"Paired{new_suffix}"
                    
                    # For single character F/R suffixes (like IlluminaF/IlluminaR)
                    elif suffix1 in ['F', 'R', 'L'] and len(suffix1) == 1:
                        new_col_name = f"{clean_base}x" if clean_base else "Pairedx"
                    
                    # For numeric patterns and all other patterns, use appropriate separator with 'x'
                    else:
                        # Determine the separator based on the suffix pattern
                        if suffix1.startswith('.'):
                            separator = '.'
                        elif suffix1.startswith('_'):
                            separator = '_'
                        else:
                            separator = '_'  # Default separator
                        
                        new_col_name = f"{clean_base}{separator}x" if clean_base else f"Paired{separator}x"
                    
                    # Merge the columns
                    merged_values = []
                    for idx, row in df.iterrows():
                        val1 = row[col]
                        val2 = row[paired_col]
                        
                        # Handle None/NaN values
                        if pd.isna(val1) or val1 == 'None':
                            val1 = 'None'
                        if pd.isna(val2) or val2 == 'None':
                            val2 = 'None'
                            
                        # If both are None, keep as None
                        if val1 == 'None' and val2 == 'None':
                            merged_values.append('None')
                        else:
                            # Create paired value
                            merged_values.append(f"{val1}, {val2}")
                    
                    # Add merged column
                    merged_df[new_col_name] = merged_values
                    
                    # Mark original columns for dropping
                    columns_to_drop.extend([col, paired_col])
                    processed_columns.update([col, paired_col])
                    
                    # Mark this as a paired column
                    if not hasattr(self, 'merged_paired_columns'):
                        self.merged_paired_columns = {}
                    self.merged_paired_columns[new_col_name] = True
                    
                    # Break after finding a match
                    break
        
        # Drop the original paired columns
        if columns_to_drop:
            merged_df = merged_df.drop(columns=columns_to_drop)
            
            # Log what we did
            self.warnings.append(f"Merged {len(columns_to_drop)//2} paired column sets with _x suffix")
        
        return merged_df


    def _apply_inferred_types(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Apply inferred types to the dataframe for consistent processing.
        RESPECTS attributes file settings when available.
        
        Args:
            df: DataFrame to process
            
        Returns:
            DataFrame with values converted to their inferred types
        """
        # Create a copy to avoid modifying the original
        df_out = df.copy()
        
        # Infer column types
        inferred_types = self._infer_column_types(df)
        
        # Store the inferred types for later use
        self.inferred_column_types = inferred_types
        
        # CRITICAL FIX: Check if we should respect attributes file
        # Get attributes if not overwriting
        attributes = {}
        if not self.analysis_results.get("file_info", {}).get("overwrite_attributes", False):
            file_path = self.analysis_results.get("file_info", {}).get("path", "")
            custom_attrs_path = self.analysis_results.get("file_info", {}).get("attributes_path", None)
            attributes = self.read_attributes_file(file_path, custom_attrs_path) if file_path else {}
        
        # Process each column based on its type
        for column, col_type in inferred_types.items():
            # Check if attributes file has a type for this column
            col_attrs = attributes.get(column, {})
            col_type_from_attrs = col_attrs.get("type", "unknown")
            
            # If attributes file specifies a type, DON'T apply inferred type conversions
            # This prevents the bug where empty values get converted to False before
            # we know the column should be treated as a path
            if col_type_from_attrs != "unknown":
                # Skip type conversion - let _apply_defaults_selectively handle it later
                continue
                
            # Only apply conversions for columns without explicit attributes
            if col_type == "boolean":
                # Convert to boolean
                for i, value in enumerate(df[column]):
                    if pd.isna(value) or value == '' or value == '-':
                        df_out.at[i, column] = False
                    else:
                        df_out.at[i, column] = self._to_boolean(value)
                        
            elif col_type == "integer":
                # Convert to integer - FIXED to preserve original numeric values
                for i, value in enumerate(df[column]):
                    if pd.isna(value) or value == '' or value == '-' or isinstance(value, bool):
                        df_out.at[i, column] = 0
                    else:
                        try:
                            # Preserve the original numeric value
                            numeric_value = float(value)
                            # Store as integer if it's a whole number
                            if numeric_value.is_integer():
                                df_out.at[i, column] = int(numeric_value)
                            else:
                                df_out.at[i, column] = numeric_value
                        except:
                            if isinstance(value, str) and value.lower() == 'auto':
                                df_out.at[i, column] = 'auto'
                            else:
                                df_out.at[i, column] = 0
                                self.warnings.append(f"Row {i+1}: Non-integer value '{value}' in {column} column, defaulting to 0")
            
            # Don't change path or string columns - keep them as is
        
        return df_out


    def _parse_read_paths(self, paths_str, read_type):
        """
        Parse read paths according to the read type.
        For paired-end reads (Hi-C, short, 10x, Illumina), group them into pairs.
        Uses Path1, Path2, etc. for naming.
        
        Special case: Accession numbers (SRR*, ERR*, DRR*, GCA_*, GCF_*) are 
        never paired - each accession gets its own Path entry.
        """
        # Handle None, NaN, empty strings, and booleans
        if paths_str is None or pd.isna(paths_str) or paths_str == '' or paths_str == 'None':
            return []
        
        # Convert boolean values to strings
        if isinstance(paths_str, bool):
            return []  # Boolean values aren't valid paths
                
        # Ensure we're working with a string
        paths_str = str(paths_str).strip()
                
        # Split the string into individual paths
        if ',' in paths_str:
            paths = [p.strip() for p in paths_str.split(',')]
        else:
            paths = [paths_str.strip()]
                
        # Remove any empty paths
        paths = [p for p in paths if p and p != 'None']
            
        # If no valid paths remain, return empty list
        if not paths:
            return []
        
        # SPECIAL CASE: Check for accession numbers
        # Accession numbers should NEVER be paired - each gets its own Path entry
        if self._contains_accession_numbers(paths):
            result = {}
            for i, path in enumerate(paths):
                result[f"Path{i+1}"] = path
            return result
        
        # Expanded list of paired-end sequencing types with variations
        paired_types = [
            'hi-c', 'hic', 'hi_c', 'hi c',  # Hi-C variations
            'short', 'short-read',          # Short read variations
            '10x', '10-x', 'ten-x',         # 10X Genomics variations
            'illumina', 'nextseq', 'miseq', 'hiseq', 'novaseq',  # Illumina platform variations
            'paired', 'paired-end', 'pe'    # Generic paired-end indicators
        ]
        
        # Check if this is a Hi-C or other paired read type based on read_type
        is_paired_by_type = False
        if isinstance(read_type, str):
            read_type_lower = read_type.lower()
            for pt in paired_types:
                if pt in read_type_lower:
                    is_paired_by_type = True
                    break
        
        # Check for R1/R2 pattern in file names to detect paired reads automatically
        r1_files = [p for p in paths if '_R1' in p or '_1.' in p or '.1.' in p or '_r1' in p or '_1_' in p]
        r2_files = [p for p in paths if '_R2' in p or '_2.' in p or '.2.' in p or '_r2' in p or '_2_' in p]
        
        # Check if we have a pattern that suggests pairing
        has_pairing_pattern = False
        if len(r1_files) > 0 and len(r2_files) > 0 and len(r1_files) == len(r2_files):
            has_pairing_pattern = True
        
        # Check for even number of files that might be pairs
        could_be_paired = len(paths) % 2 == 0 and len(paths) >= 2
        
        # Enhanced pair detection: look for files that differ only by 1/2 pattern
        detected_pairs = []
        if could_be_paired and not has_pairing_pattern:
            # Try to match files based on similar names with 1/2 differences
            remaining_files = paths.copy()
            
            for file1 in paths:
                if file1 not in remaining_files:
                    continue
                    
                # Look for a file that matches this one but with 1->2 or similar
                for file2 in remaining_files:
                    if file1 == file2:
                        continue
                        
                    # Check if they form a 1/2 pair
                    if self._are_paired_files(file1, file2):
                        detected_pairs.append((file1, file2))
                        remaining_files.remove(file1)
                        remaining_files.remove(file2)
                        break
        
        # Determine if we should treat this as paired data
        should_treat_as_paired = is_paired_by_type or has_pairing_pattern or len(detected_pairs) > 0
        
        # If we have matched R1/R2 files using the original pattern
        if should_treat_as_paired and has_pairing_pattern:
            # Sort to ensure R1 and R2 from same library are paired
            r1_files.sort()
            r2_files.sort()
            
            # Group into pairs
            result = {}
            for i, (r1, r2) in enumerate(zip(r1_files, r2_files)):
                result[f"Path{i+1}"] = f"{r1}, {r2}"
            return result
        
        # If we detected pairs through name matching
        elif len(detected_pairs) > 0:
            result = {}
            for i, (file1, file2) in enumerate(detected_pairs):
                result[f"Path{i+1}"] = f"{file1}, {file2}"
            return result
        
        # If we have an even number of files and pairing is suggested by read type
        elif should_treat_as_paired and could_be_paired:
            result = {}
            for i in range(0, len(paths), 2):
                result[f"Path{(i//2)+1}"] = f"{paths[i]}, {paths[i+1]}"
            return result
        
        # Otherwise list as individual files
        else:
            result = {}
            for i, path in enumerate(paths):
                result[f"Path{i+1}"] = path
            return result


    def _are_paired_files(self, file1: str, file2: str) -> bool:
        """
        Check if two files are likely paired reads (e.g., lib1_1.fq and lib1_2.fq)
        """
        # Common paired patterns
        pair_patterns = [
            ('_1.', '_2.'), ('.1.', '.2.'), ('_1_', '_2_'),
            ('_R1', '_R2'), ('_r1', '_r2'), ('.R1', '.R2'),
            ('_READ1', '_READ2'), ('_read1', '_read2'),
            ('_F.', '_R.'), ('_forward', '_reverse'),
            ('_LEFT', '_RIGHT'), ('_left', '_right')
        ]
        
        for pattern1, pattern2 in pair_patterns:
            # Check if file1 has pattern1 and file2 has pattern2
            if pattern1 in file1 and pattern2 in file2:
                # Check if the rest of the filename is the same
                base1 = file1.replace(pattern1, '__PAIR__')
                base2 = file2.replace(pattern2, '__PAIR__')
                if base1 == base2:
                    return True
            
            # Also check the reverse (file2 has pattern1, file1 has pattern2)
            if pattern1 in file2 and pattern2 in file1:
                base1 = file1.replace(pattern2, '__PAIR__')
                base2 = file2.replace(pattern1, '__PAIR__')
                if base1 == base2:
                    return True
        
        return False


    def _read_numbers_file(self, file_path: str) -> pd.DataFrame:
        """
        Read Apple Numbers file and convert to pandas DataFrame.
        """
        try:
            # Parse the Numbers document
            doc = Document(file_path)
            
            # Check if document has sheets
            if not doc.sheets:
                self.errors.append("Numbers file contains no sheets")
                return pd.DataFrame()
            
            # If multiple sheets, warn user and use first one
            if len(doc.sheets) > 1:
                sheet_names = [sheet.name for sheet in doc.sheets]
                self.warnings.append(f"Numbers file contains {len(doc.sheets)} sheets: {', '.join(sheet_names)}. Using first sheet: '{sheet_names[0]}'")
            
            sheet = doc.sheets[0]
            
            # Check if sheet has tables
            if not sheet.tables:
                self.errors.append(f"Sheet '{sheet.name}' contains no tables")
                return pd.DataFrame()
            
            # If multiple tables, warn user and use first one
            if len(sheet.tables) > 1:
                self.warnings.append(f"Sheet contains {len(sheet.tables)} tables. Using first table.")
            
            table = sheet.tables[0]
            
            # Handle empty table
            if table.num_rows == 0:
                self.warnings.append("Numbers table is empty")
                return pd.DataFrame()
            
            # Extract headers from first row
            headers = []
            for col_idx in range(table.num_cols):
                try:
                    cell_value = table.cell(0, col_idx).value
                    if cell_value is None or str(cell_value).strip() == '':
                        header = f"Column_{col_idx + 1}"
                    else:
                        header = str(cell_value).strip()
                    headers.append(header)
                except Exception as e:
                    headers.append(f"Column_{col_idx + 1}")
                    self.warnings.append(f"Error reading header for column {col_idx + 1}: {str(e)}")
            
            # Extract data rows
            data_rows = []
            for row_idx in range(1, table.num_rows):
                row_data = []
                for col_idx in range(table.num_cols):
                    try:
                        cell_value = table.cell(row_idx, col_idx).value
                        row_data.append(cell_value)
                    except Exception as e:
                        row_data.append(None)
                        self.warnings.append(f"Error reading cell at row {row_idx + 1}, column {col_idx + 1}: {str(e)}")
                data_rows.append(row_data)
            
            # Create DataFrame
            if data_rows:
                df = pd.DataFrame(data_rows, columns=headers)
            else:
                df = pd.DataFrame(columns=headers)
            
            return df
            
        except Exception as e:
            self.errors.append(f"Failed to parse Numbers file: {str(e)}")
            return pd.DataFrame()


    def _read_libreoffice_file(self, file_path: str) -> pd.DataFrame:
        """
        Read LibreOffice files (.ods, .ots) and convert to pandas DataFrame.
        """
        file_ext = os.path.splitext(file_path)[1].lower()
        
        # Method 1: Try pandas direct reading (works for .ods)
        if file_ext == '.ods':
            try:
                df = pd.read_excel(file_path, engine='odf')
                
                # Clean whitespace
                for column in df.columns:
                    df[column] = df[column].apply(
                        lambda x: x.strip() if isinstance(x, str) else x
                    )
                df.columns = [col.strip() for col in df.columns]
                
                return df
                
            except Exception as e:
                self.warnings.append(f"pandas ODS reading failed: {str(e)}, trying alternative method")
        
        # Method 2: Use odfpy directly for more control
        try:
            # Load the ODF document
            doc = load(file_path)
            
            # Get all tables (sheets)
            tables = doc.getElementsByType(Table)
            
            if not tables:
                self.errors.append("No tables found in LibreOffice file")
                return pd.DataFrame()
            
            # Use first table by default
            if len(tables) > 1:
                self.warnings.append(f"LibreOffice file contains {len(tables)} tables/sheets. Using first one.")
            
            table = tables[0]
            
            # Extract data from the table
            rows_data = []
            table_rows = table.getElementsByType(TableRow)
            
            for row in table_rows:
                row_data = []
                cells = row.getElementsByType(TableCell)
                
                for cell in cells:
                    # Handle repeated columns (LibreOffice optimization)
                    repeat_attr = cell.getAttribute('numbercolumnsrepeated')
                    repeat_count = int(repeat_attr) if repeat_attr else 1
                    
                    # Extract text from cell
                    cell_value = ""
                    paragraphs = cell.getElementsByType(P)
                    if paragraphs:
                        cell_texts = []
                        for p in paragraphs:
                            # Get all text content
                            text_content = ""
                            for node in p.childNodes:
                                if hasattr(node, 'data'):
                                    text_content += node.data
                            cell_texts.append(text_content)
                        cell_value = " ".join(cell_texts).strip()
                    
                    # Handle empty cells and repetition
                    if not cell_value:
                        cell_value = None
                    
                    # Add repeated columns
                    for _ in range(repeat_count):
                        row_data.append(cell_value)
                
                if row_data:  # Only add non-empty rows
                    rows_data.append(row_data)
            
            if not rows_data:
                self.warnings.append("LibreOffice table appears to be empty")
                return pd.DataFrame()
            
            # Determine number of columns (use the longest row)
            max_cols = max(len(row) for row in rows_data) if rows_data else 0
            
            # Pad shorter rows with None values
            for row in rows_data:
                while len(row) < max_cols:
                    row.append(None)
            
            # First row as headers, rest as data
            if len(rows_data) > 0:
                headers = []
                for i, header in enumerate(rows_data[0]):
                    if header is None or str(header).strip() == '':
                        headers.append(f"Column_{i + 1}")
                    else:
                        headers.append(str(header).strip())
                
                # Data rows (skip header)
                data_rows = rows_data[1:] if len(rows_data) > 1 else []
                
                # Create DataFrame
                if data_rows:
                    df = pd.DataFrame(data_rows, columns=headers)
                else:
                    df = pd.DataFrame(columns=headers)
                
                # Clean whitespace
                for column in df.columns:
                    df[column] = df[column].apply(
                        lambda x: x.strip() if isinstance(x, str) else x
                    )
                
                return df
            
            return pd.DataFrame()
            
        except Exception as e:
            self.errors.append(f"Failed to read LibreOffice file: {str(e)}")
            return pd.DataFrame()


    def read_file(self, file_path: str) -> pd.DataFrame:
        """
        Read a table file with automatic format and delimiter detection.
        Supports Excel, Numbers, LibreOffice, CSV, and text files with various delimiters.
        """
        # Store file path in analysis results right away
        self.analysis_results["file_info"]["path"] = file_path
        
        file_ext = os.path.splitext(file_path)[1].lower()
        
        # Excel files
        if file_ext in ['.xlsx', '.xls', '.xlsm', '.xltx', '.xltm']:
            try:
                # Use appropriate engine based on file type
                engine = 'openpyxl' if file_ext in ['.xlsx', '.xlsm', '.xltx', '.xltm'] else 'xlrd'
                df = pd.read_excel(file_path, engine=engine)
            except Exception as e:
                self.errors.append(f"Failed to read Excel file: {str(e)}")
                return pd.DataFrame()
                
            self.analysis_results["file_info"]["format"] = "Microsoft Excel"
            self.analysis_results["file_info"]["delimiter"] = "N/A"
            
            # Clean whitespace from all string columns
            for column in df.columns:
                df[column] = df[column].apply(
                    lambda x: x.strip() if isinstance(x, str) else x
                )
                
            # Also clean column names
            df.columns = [col.strip() for col in df.columns]
            
            return df
        
        # Apple Numbers files
        elif file_ext == '.numbers':
            df = self._read_numbers_file(file_path)
            if not df.empty:
                self.analysis_results["file_info"]["format"] = "Apple Numbers"
                self.analysis_results["file_info"]["delimiter"] = "N/A"
                
                # Clean whitespace from all string columns
                for column in df.columns:
                    df[column] = df[column].apply(
                        lambda x: x.strip() if isinstance(x, str) else x
                    )
                
                # Also clean column names  
                df.columns = [col.strip() for col in df.columns]
            
            return df
        
        # LibreOffice files (.ods = Calc spreadsheet, .ots = Calc template)
        elif file_ext in ['.ods', '.ots', '.fods']:
            df = self._read_libreoffice_file(file_path)
            if not df.empty:
                format_names = {
                    '.ods': "LibreOffice Calc",
                    '.ots': "LibreOffice Calc Template", 
                    '.fods': "Flat OpenDocument Spreadsheet"
                }
                self.analysis_results["file_info"]["format"] = format_names.get(file_ext, "LibreOffice")
                self.analysis_results["file_info"]["delimiter"] = "N/A"
            return df
        
        # For text files, detect and fix malformed headers FIRST
        try:
            with open(file_path, 'r', errors='ignore') as f:
                lines = f.readlines()
                
            if not lines:
                self.errors.append("File is empty")
                return pd.DataFrame()
            
            # Fix malformed headers before any other processing
            lines = self._fix_malformed_headers(lines)
            
            # Write fixed lines to a temporary file
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as temp_file:
                temp_file.writelines(lines)
                temp_path = temp_file.name
            
            # Now use pandas to read the fixed file
            delimiter = self._detect_delimiter(temp_path)
            self.analysis_results["file_info"]["format"] = "Delimited Text"
            self.analysis_results["file_info"]["delimiter"] = delimiter
            
            try:
                df = pd.read_csv(temp_path, sep=delimiter, engine='python')
                os.unlink(temp_path)  # Clean up temp file
                
                # Clean whitespace from all string columns
                for column in df.columns:
                    df[column] = df[column].apply(
                        lambda x: x.strip() if isinstance(x, str) else x
                    )
                    
                # Also clean column names
                df.columns = [col.strip() for col in df.columns]
                
                return df
                
            except Exception as e:
                os.unlink(temp_path)  # Clean up temp file
                # Fall back to original file reading logic
                pass
                
        except Exception:
            # If header fixing fails, continue with original logic
            pass
        
        # CSV and other delimited files - standard approach
        delimiter = self._detect_delimiter(file_path)
        self.analysis_results["file_info"]["format"] = "Delimited Text"
        self.analysis_results["file_info"]["delimiter"] = delimiter
        
        # Try multiple approaches to read the file
        try:
            df = pd.read_csv(file_path, sep=delimiter, engine='python')
        except Exception as e:
            try:
                df = pd.read_csv(file_path, sep=delimiter, quoting=csv.QUOTE_NONE)
            except Exception as e:
                try:
                    # Try common delimiters with skippinitialspace
                    for delim in [',', ';', '\t', '|']:
                        try:
                            df = pd.read_csv(file_path, sep=delim, skipinitialspace=True)
                            if len(df.columns) > 1:  # Successfully parsed into multiple columns
                                self.warnings.append(f"Used delimiter '{delim}' with skipped spaces")
                                self.analysis_results["file_info"]["delimiter"] = delim
                                break
                        except:
                            continue
                except Exception:
                    pass
                
                # Last resort: very permissive reading
                try:
                    df = pd.read_csv(file_path, sep=delimiter, quoting=csv.QUOTE_NONE, 
                                    on_bad_lines='skip')
                except Exception as e:
                    self.errors.append(f"Failed to read file: {str(e)}")
                    # Return empty DataFrame as a fallback
                    return pd.DataFrame()

        # Clean column names more aggressively
        cleaned_columns = []
        for col in df.columns:
            # Replace multiple spaces with single space
            cleaned_col = ' '.join(col.split())
            cleaned_columns.append(cleaned_col.strip())
        df.columns = cleaned_columns

        # Clean whitespace from all string columns
        for column in df.columns:
            df[column] = df[column].apply(
                lambda x: x.strip() if isinstance(x, str) else x
            )
        
        return df

    
    def analyze_string_table(self, table_str: str) -> Dict[str, Any]:
        """
        Analyze a table provided as a string.
        Attempts to detect the delimiter and parse the table.
        """
        # Create a temporary file to leverage existing file reading logic
        import tempfile
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp:
            tmp.write(table_str)
            tmp_path = tmp.name
        
        try:
            results = self.analyze_file(tmp_path)
            # Override the file path to indicate this was from a string
            if "file_info" in results:
                results["file_info"]["path"] = "<string input>"
            return results
        finally:
            os.unlink(tmp_path)  # Clean up the temporary file


    def analyze_file(self, file_path: str, attributes_file_path: str = None, 
                sanitise_fields: bool = False, overwrite_attributes: bool = False,
                validate_paths: bool = False, discover_paths: bool = False) -> Dict[str, Any]:
        """
        Analyze a table file and return detailed information about its structure,
        validate data according to specific requirements, and prepare a hierarchical
        data structure.
        
        Args:
            file_path: Path to the table file
            attributes_file_path: Optional path to a custom attributes file
            sanitise_fields: Whether to sanitize hierarchical field values
            overwrite_attributes: Whether to regenerate attributes file ignoring existing one
            validate_paths: Whether to validate that all file paths in the output exist
            discover_paths: Whether to expand glob patterns (*, ?, []) in path columns to actual file paths
        """
        # Reset state
        self.warnings = []
        self.errors = []
        self.structured_data = {}
        self.analysis_results = {"file_info": {}, "warnings": [], "errors": []}
        self.inferred_column_types = {}
        
        # Store the custom attributes file path if provided
        if attributes_file_path:
            self.analysis_results["file_info"]["attributes_path"] = attributes_file_path
        
        # Check attributes file BEFORE any processing, and store the result
        # If --overwrite-attributes is used, act as if no attributes file exists for detection
        custom_attrs_path = attributes_file_path
        default_attrs_path = os.path.splitext(file_path)[0] + ".attributes"
        
        # Check by trying to open the files - most reliable method
        attributes_exist = False
        if not overwrite_attributes:  # Only check if not overwriting
            try:
                if custom_attrs_path:
                    with open(custom_attrs_path, 'r') as f:
                        attributes_exist = True
                else:
                    with open(default_attrs_path, 'r') as f:
                        attributes_exist = True
            except (IOError, FileNotFoundError):
                attributes_exist = False
        
        # Store the initial detection state for later use
        self.analysis_results["file_info"]["attributes_found_initially"] = attributes_exist
        
        # Store overwrite flag for attributes generation
        self.analysis_results["file_info"]["overwrite_attributes"] = overwrite_attributes
        
        try:
            # Read the file
            df = self.read_file(file_path)

            # Expand glob patterns if requested (BEFORE merging paired columns)
            if discover_paths:
                input_dir = os.path.dirname(os.path.abspath(file_path))
                df = self._expand_glob_patterns(df, input_dir)

            # Merge paired columns before further processing
            df = self._merge_paired_columns(df)
            
            # Store the DataFrame in the object for later access
            self.current_df = df
            
            # Update file info with basic stats
            self.analysis_results["file_info"]["rows"] = len(df)
            self.analysis_results["file_info"]["columns"] = len(df.columns)
            
            # If the DataFrame is empty, we had an error reading the file
            if df.empty and self.errors:
                self.analysis_results["errors"] = self.errors
                return self.analysis_results
            
            # Replace empty strings with None for better handling
            df = df.infer_objects(copy=False).replace(['', ' '], np.nan)
            
            # Infer column types for better data handling
            df = self._apply_inferred_types(df)
            
            # Detect paired data columns
            paired_columns = self._infer_paired_data_columns(df)
            
            # Detect read type column
            read_type_column = self._detect_read_type_column(df)
            
            # Store this information in analysis results
            self.analysis_results["column_info"] = {
                "inferred_types": self.inferred_column_types,
                "paired_columns": paired_columns,
                "read_type_column": read_type_column
            }
            
            # Basic validation checks - use fresh detection if overwriting attributes
            if overwrite_attributes:
                # When overwriting, act as if no attributes file exists for validation
                self._validate_data(df, None)
            else:
                self._validate_data(df, attributes_file_path)
            
            # Update errors in the results
            self.analysis_results["errors"] = self.errors
            
            # If there are errors, return early
            if self.errors:
                return self.analysis_results
            
            # Apply default values for certain fields (using the selective method)
            df_with_defaults = self._apply_defaults_selectively(df, attributes_file_path)
            
            # Create df_clean ALWAYS - this ensures it's always defined
            df_clean = df_with_defaults.copy()
            
            # Convert object columns to strings for consistent grouping
            for col in df_clean.columns:
                if df_clean[col].dtype == 'object':
                    df_clean[col] = df_clean[col].astype(str)
            
            # Identify hierarchy columns based on data patterns
            hierarchy_columns = self._detect_hierarchy_columns(df_clean)

            # Handle field sanitization - always check attributes file for sanitise flags
            # But ignore existing attributes if overwrite flag is set
            if overwrite_attributes:
                # Act as if no attributes file exists - use fresh detection only
                attributes = {}
            else:
                # Read attributes file to check for individual column sanitise flags
                attributes = self.read_attributes_file(file_path, attributes_file_path)
            
            # Apply sanitization based on:
            # 1. Individual column 'sanitise' flags in .attributes file (if not overwriting)
            # 2. --sanitise-fields flag (sanitizes ALL hierarchy columns)
            df_clean = self._sanitize_hierarchical_values(df_clean, hierarchy_columns, 
                                                        attributes, sanitise_fields)

            # Build structured hierarchical data
            self._build_structured_data(df_clean, paired_columns, read_type_column)
            
            # Validate paths if requested
            if validate_paths:
                input_dir = os.path.dirname(os.path.abspath(file_path))
                paths_valid = self._validate_paths_in_structured_data(
                    self.structured_data, 
                    input_dir
                )
                
                # If validation failed, update errors and return early (do not generate YAML)
                if not paths_valid:
                    self.analysis_results["errors"] = self.errors
                    self.analysis_results["warnings"] = self.warnings
                    return self.analysis_results
            
            # Update the results
            self.analysis_results["warnings"] = self.warnings
            self.analysis_results["structured_data"] = self.structured_data
            
            # Generate attributes file - note this creates the file if it doesn't exist
            # or overwrites if --overwrite-attributes flag is used
            try:
                # Use the original df for attributes generation, not df_clean
                self.generate_attributes_file(file_path, df, attributes_file_path,
                                            sanitise_fields, overwrite_attributes)
            except Exception as e:
                self.warnings.append(f"Could not generate attributes file: {str(e)}")
            
            return self.analysis_results
            
        except Exception as e:
            # Catch any unexpected errors and provide useful debugging info
            error_msg = f"Unexpected error during analysis: {str(e)}"
            self.errors.append(error_msg)
            self.analysis_results["errors"] = self.errors
            
            # Print stack trace for debugging
            import traceback
            print(f"DEBUG: {error_msg}")
            traceback.print_exc()
            
            return self.analysis_results


    def _apply_defaults_selectively(self, df: pd.DataFrame, attributes_file_path: str = None, 
                                overwrite_attributes: bool = False) -> pd.DataFrame:
        """
        Apply default values to columns based on their data types, but preserve
        existing numeric values. Respects attributes file settings as primary authority.
        
        Args:
            df: DataFrame to process
            attributes_file_path: Optional path to a custom attributes file
            
        Returns:
            DataFrame with default values applied where appropriate
        """
        # Create a copy to avoid modifying the original
        df_out = df.copy()
        
        # Check if there's an attributes file to use for defaults
        # But ignore if overwrite flag is set
        if overwrite_attributes:
            attributes = {}  # Act as if no attributes file exists
        else:
            attributes = self.read_attributes_file(self.analysis_results["file_info"]["path"], attributes_file_path)
        
        # Process each column
        for column in df.columns:
            # Get column attributes if available, otherwise use built-in logic
            col_attrs = attributes.get(column, {})
            col_type = col_attrs.get("type", "unknown")
            col_default = col_attrs.get("default", None)
            
            # PRIMARY: Handle columns based on attributes file type (if specified)
            if col_type != "unknown":
                if col_type == "boolean":
                    for i, value in enumerate(df[column]):
                        # Skip if it's already a boolean
                        if isinstance(value, bool):
                            df_out.at[i, column] = value
                            continue
                            
                        if pd.isna(value) or value == '' or value == '-' or not self._is_boolean(value):
                            # Use default from attributes file
                            default_val = False
                            if col_default is not None:
                                if col_default.lower() == 'none':
                                    default_val = 'None'  # String 'None', not boolean
                                else:
                                    default_val = str(col_default).lower() in ['true', 'yes', 'y', 't', '1']
                            df_out.at[i, column] = default_val
                        else:
                            # Keep the actual value as-is for boolean columns
                            bool_value = str(value).lower() in ['true', 'yes', 'y', 't', '1', 'on']
                            df_out.at[i, column] = bool_value
                            
                elif col_type == "integer":
                    for i, value in enumerate(df[column]):
                        # Skip if it's already an integer
                        if isinstance(value, int) and not isinstance(value, bool):
                            continue
                        
                        # For existing numeric values, try to convert but preserve the value
                        if not pd.isna(value) and value != '' and value != '-' and not isinstance(value, bool):
                            try:
                                numeric_val = float(value)
                                if numeric_val.is_integer():
                                    df_out.at[i, column] = int(numeric_val)
                                else:
                                    df_out.at[i, column] = numeric_val
                                continue
                            except:
                                pass
                        
                        # For empty or non-numeric values, apply default from attributes
                        default_val = 0
                        if col_default is not None:
                            if col_default.lower() == 'none':
                                default_val = 'None'
                            else:
                                try:
                                    default_val = int(col_default)
                                except:
                                    default_val = 0
                        df_out.at[i, column] = default_val
                        
                elif col_type in ["path", "simple_path", "complex_path"]:
                    # Handle path types - apply defaults for empty fields
                    for i, value in enumerate(df[column]):
                        if pd.isna(value) or value == '' or value == '-':
                            # Use default from attributes file
                            if col_default is not None:
                                df_out.at[i, column] = col_default
                            else:
                                df_out.at[i, column] = "None"
                        elif isinstance(value, str):
                            # Strip whitespace from string values
                            df_out.at[i, column] = value.strip()
                            
                elif col_type in ["string", "category"]:
                    # CRITICAL FIX: Handle string and category types properly
                    for i, value in enumerate(df[column]):
                        if pd.isna(value) or value == '' or value == '-':
                            # Use default from attributes file
                            if col_default is not None:
                                df_out.at[i, column] = col_default
                            else:
                                df_out.at[i, column] = "None"
                        else:
                            # Convert to string and strip whitespace
                            # IMPORTANT: Don't convert pd.NA/np.nan to string ".nan"
                            df_out.at[i, column] = str(value).strip()
                
                # For unknown types with attributes, still handle empty values
                else:
                    for i, value in enumerate(df[column]):
                        if pd.isna(value) or value == '' or value == '-':
                            if col_default is not None:
                                df_out.at[i, column] = col_default
                            else:
                                df_out.at[i, column] = "None"
                        elif isinstance(value, str):
                            df_out.at[i, column] = value.strip()
                
            # FALLBACK: Handle columns with unknown type using built-in logic
            else:
                # ... (keep existing fallback logic)
                # Handle boolean fields based on column name patterns
                if column.lower() in ['trim', 'read_qc', 'trim_qc']:
                    for i, value in enumerate(df[column]):
                        if isinstance(value, bool):
                            df_out.at[i, column] = value
                            continue
                            
                        if pd.isna(value) or value == '' or value == '-' or not self._is_boolean(value):
                            default_val = False
                            if col_default is not None:
                                if col_default.lower() == 'none':
                                    default_val = 'None'
                                else:
                                    default_val = str(col_default).lower() in ['true', 'yes', 'y', 't', '1']
                            self.warnings.append(f"Row {i+1}: Empty or non-boolean value '{value}' in {column} column, defaulting to {default_val}")
                            df_out.at[i, column] = default_val
                        else:
                            bool_value = str(value).lower() in ['true', 'yes', 'y', 't', '1', 'on']
                            df_out.at[i, column] = bool_value
                
                # Handle integer fields based on column name patterns
                elif column.lower() in ['kmer_size', 'genome_size'] or column.lower().endswith('_size'):
                    for i, value in enumerate(df[column]):
                        if isinstance(value, int) and not isinstance(value, bool):
                            continue
                        
                        if not pd.isna(value) and value != '' and value != '-' and not isinstance(value, bool):
                            try:
                                numeric_val = float(value)
                                if numeric_val.is_integer():
                                    df_out.at[i, column] = int(numeric_val)
                                else:
                                    df_out.at[i, column] = numeric_val
                                continue
                            except:
                                pass
                        
                        default_val = 0
                        if col_default is not None:
                            if col_default.lower() == 'none':
                                default_val = 'None'
                            else:
                                try:
                                    default_val = int(col_default)
                                except:
                                    default_val = 0
                        self.warnings.append(f"Row {i+1}: Non-integer value '{value}' in {column} column, defaulting to {default_val}")
                        df_out.at[i, column] = default_val
                
                # Handle string fields
                else:
                    for i, value in enumerate(df[column]):
                        if pd.isna(value) or value == '' or value == '-':
                            if col_default is not None:
                                df_out.at[i, column] = col_default
                            else:
                                if pd.api.types.is_numeric_dtype(df[column].dtype):
                                    df_out.at[i, column] = np.nan
                                else:
                                    df_out.at[i, column] = "None"
                        elif isinstance(value, str):
                            df_out.at[i, column] = value.strip()
            
        return df_out


    def _validate_data(self, df: pd.DataFrame, attributes_file_path: str = None) -> None:
        """
        Validate data according to specific requirements and attribute constraints
        
        Args:
            df: DataFrame to validate
            attributes_file_path: Optional path to a custom attributes file
        """
        
        # Check if there's an attributes file to use for validation
        attributes = self.read_attributes_file(self.analysis_results["file_info"]["path"], attributes_file_path)
        
        # Check for required columns from attributes file
        mandatory_columns = []
        for col, attrs in attributes.items():
            if attrs.get("mandatory", False):
                mandatory_columns.append(col)
                
        # If no mandatory columns are defined in attributes, default to first column only
        if not mandatory_columns and len(df.columns) > 0:
            first_column = df.columns[0]
            mandatory_columns = [first_column]
            
            # Ultra reliable file check by trying to open them
            custom_exists = False
            default_exists = False
            
            if attributes_file_path:
                try:
                    with open(attributes_file_path, 'r') as f:
                        custom_exists = True
                except (IOError, FileNotFoundError):
                    custom_exists = False
                    
            default_attrs_path = os.path.splitext(self.analysis_results["file_info"]["path"])[0] + ".attributes"
            try:
                with open(default_attrs_path, 'r') as f:
                    default_exists = True
            except (IOError, FileNotFoundError):
                default_exists = False
            
            # Only add warning if no attributes file actually exists
            if not (custom_exists or default_exists):
                self.warnings.append(f"Since no attributes file found, assuming only the first column '{first_column}' is mandatory.")
                        
        # Check for missing mandatory columns
        for col in mandatory_columns:
            if col not in df.columns:
                self.errors.append(f"Missing required column: {col}")
        
        # If any mandatory columns are missing, return early
        if any(col not in df.columns for col in mandatory_columns):
            return
                    
        # Check for empty fields in mandatory columns
        for col in mandatory_columns:
            if df[col].isna().any():
                rows = df[df[col].isna()].index.tolist()
                rows_str = ', '.join(str(r+1) for r in rows)  # +1 for human-readable row numbers
                self.errors.append(f"Empty {col} found in rows: {rows_str}")
                    
        # Validate against specific attribute constraints if available
        for col, attrs in attributes.items():
            if col not in df.columns:
                continue
                
            # Check category values if defined
            if attrs["type"] == "category" and attrs["values"]:
                valid_values = set(attrs["values"])
                for i, value in enumerate(df[col]):
                    if pd.notna(value) and str(value).strip() not in valid_values:
                        self.warnings.append(f"Row {i+1}: Value '{value}' in column {col} is not one of the allowed values: {', '.join(valid_values)}")

 

    def _apply_defaults(self, df: pd.DataFrame, attributes_file_path: str = None) -> pd.DataFrame:
        """
        Apply default values to columns based on their data types
        
        Args:
            df: DataFrame to process
            attributes_file_path: Optional path to a custom attributes file
        """
        
        # Create a copy to avoid modifying the original
        df_out = df.copy()
        
        # Check if there's an attributes file to use for defaults
        attributes = self.read_attributes_file(self.analysis_results["file_info"]["path"], attributes_file_path)
        
        # Process each column
        for column in df.columns:
            # Get column attributes if available, otherwise use built-in logic
            col_attrs = attributes.get(column, {})
            col_type = col_attrs.get("type", "unknown")
            col_default = col_attrs.get("default", None)
            
            # Handle boolean fields - special focus on Trim and Read_QC fields
            if col_type == "boolean" or column.lower() in ['trim', 'read_qc']:
                for i, value in enumerate(df[column]):
                    # Skip if it's already a boolean
                    if isinstance(value, bool):
                        df_out.at[i, column] = value
                        continue
                        
                    if pd.isna(value) or value == '' or value == '-' or not self._is_boolean(value):
                        default_val = False
                        if col_default is not None:
                            # Convert string default to boolean
                            default_val = str(col_default).lower() in ['true', 'yes', 'y', 't', '1']
                        self.warnings.append(f"Row {i+1}: Empty or non-boolean value '{value}' in {column} column, defaulting to {default_val}")
                        df_out.at[i, column] = default_val
                    else:
                        # Explicitly convert string boolean values to Python booleans
                        bool_value = str(value).lower() in ['true', 'yes', 'y', 't', '1', 'on']
                        df_out.at[i, column] = bool_value
            
            # Handle integer fields - Apply stricter integer conversion
            elif col_type == "integer" or column.lower() in ['kmer_size', 'genome_size'] or column.lower().endswith('_size'):
                for i, value in enumerate(df[column]):
                    # Skip if it's already an integer
                    if isinstance(value, int) and not isinstance(value, bool):  # Explicit check to exclude boolean
                        df_out.at[i, column] = value
                        continue
                        
                    # For non-integer values, apply default
                    if pd.isna(value) or value == '' or value == '-' or not self._is_integer(value) or isinstance(value, bool):
                        default_val = 0
                        if col_default is not None:
                            try:
                                default_val = int(col_default)
                            except:
                                pass
                        self.warnings.append(f"Row {i+1}: Non-integer value '{value}' in {column} column, defaulting to {default_val}")
                        df_out.at[i, column] = default_val
                    else:
                        try:
                            df_out.at[i, column] = int(value)
                        except:
                            df_out.at[i, column] = 0
                            self.warnings.append(f"Row {i+1}: Failed to convert '{value}' to integer in {column} column, defaulting to 0")
            
            # Handle string fields - replace empty with None
            elif col_type == "string" or col_type == "path" or col_type == "unknown":
                for i, value in enumerate(df[column]):
                    if pd.isna(value) or value == '' or value == '-':
                        # Check column dtype to ensure compatible type assignment
                        if pd.api.types.is_numeric_dtype(df[column].dtype):
                            df_out.at[i, column] = np.nan  # Use NaN for numeric columns
                        else:
                            df_out.at[i, column] = "None"
                    elif isinstance(value, str):
                        # Strip whitespace from string values
                        df_out.at[i, column] = value.strip()
        
        return df_out


    def _build_structured_data(self, df_clean: pd.DataFrame, paired_columns: List[str] = None, 
                            read_type_column: str = None) -> None:
        """
        Build a structured hierarchical data representation based on detected hierarchy.
        Uses a general approach that works for any data structure.
        Ensures consistent structure by including all fields with appropriate default values.
        
        Args:
            df_clean: Pre-cleaned DataFrame with consistent types for grouping
            paired_columns: List of columns containing paired data
            read_type_column: Column containing read type information
        """
        # df_clean is now passed in as a parameter, already cleaned
        # No need to create it here anymore
        
        # Get all columns in left-to-right order
        all_columns = list(df_clean.columns)
        
        # Identify hierarchy columns based on data patterns
        hierarchy_columns = self._detect_hierarchy_columns(df_clean)
        
        # CRITICAL FIX: For simple cases like Species->Read type->properties,
        # ensure we detect the correct hierarchical structure
        if len(all_columns) <= 5:  # Simple tables
            # Check if second column should definitely be hierarchical
            if len(all_columns) >= 2:
                second_col = all_columns[1]
                second_col_lower = second_col.lower()
                
                # Strong indicators that second column should be hierarchical
                if any(indicator in second_col_lower for indicator in ['type', 'platform', 'technology', 'method']):
                    # Check if it actually varies within first column groups
                    first_col = all_columns[0]
                    varies_within_groups = False
                    
                    for _, group in df_clean.groupby(first_col):
                        unique_vals = group[second_col].dropna().unique()
                        if len(unique_vals) > 1:
                            varies_within_groups = True
                            break
                    
                    # If it varies within groups, it should definitely be hierarchical
                    if varies_within_groups and second_col not in hierarchy_columns:
                        hierarchy_columns.insert(1, second_col)
        
        # Initialize structured data
        structured_data = {}
        
        # Special handling for merged paired columns
        merged_paired_cols = []
        if hasattr(self, 'merged_paired_columns'):
            merged_paired_cols = list(self.merged_paired_columns.keys())
        
        # Determine property placement for non-hierarchy columns
        property_columns = [col for col in all_columns if col not in hierarchy_columns]
        property_placement = self._determine_property_placement(df_clean, hierarchy_columns, property_columns)
        
        # Create a template of all properties with their default values for consistency
        property_defaults = self._create_property_defaults(df_clean, property_columns, merged_paired_cols)
        
        # Group by hierarchy columns to process related rows together
        if hierarchy_columns:
            grouped = df_clean.groupby(hierarchy_columns, dropna=False)
            
            for group_keys, group_df in grouped:
                # Ensure group_keys is a tuple
                if not isinstance(group_keys, tuple):
                    group_keys = (group_keys,)
                
                # Build the hierarchy structure dynamically
                self._build_hierarchy_path(structured_data, hierarchy_columns, group_keys, 
                                        group_df, property_columns, property_placement, 
                                        merged_paired_cols, read_type_column, property_defaults)
        
        # Store the result
        self.structured_data = structured_data


    def _create_property_defaults(self, df: pd.DataFrame, property_columns: List[str], 
                                merged_paired_cols: List[str]) -> Dict[str, Any]:
        """
        Create default values for all property columns to ensure consistent structure.
        Each default is a unique object to prevent YAML anchors.
        """
        defaults = {}
        
        # Get attributes if not overwriting
        attributes = {}
        if not self.analysis_results.get("file_info", {}).get("overwrite_attributes", False):
            file_path = self.analysis_results.get("file_info", {}).get("path", "")
            custom_attrs_path = self.analysis_results.get("file_info", {}).get("attributes_path", None)
            attributes = self.read_attributes_file(file_path, custom_attrs_path) if file_path else {}
        
        for col in property_columns:
            # Check attributes file first for type
            col_attrs = attributes.get(col, {})
            col_type_from_attrs = col_attrs.get("type", "unknown")
            col_default_from_attrs = col_attrs.get("default", None)
            
            # Use attributes type if available, otherwise fall back to inferred
            if col_type_from_attrs != "unknown":
                col_type = col_type_from_attrs
            else:
                col_type = self.inferred_column_types.get(col, "string")
            
            # Create appropriate default based on type
            if col in merged_paired_cols:
                # Create a unique dict for each column to avoid anchors
                defaults[col] = {"Path1": "None"}
            elif col_type == "boolean":
                # Use default from attributes if available
                if col_default_from_attrs is not None:
                    if col_default_from_attrs.lower() == 'none':
                        defaults[col] = 'None'
                    else:
                        defaults[col] = col_default_from_attrs.lower() == 'true'
                else:
                    defaults[col] = False
            elif col_type == "integer":
                if col_default_from_attrs is not None:
                    if col_default_from_attrs.lower() == 'none':
                        defaults[col] = 'None'
                    else:
                        try:
                            defaults[col] = int(col_default_from_attrs)
                        except:
                            defaults[col] = 0
                else:
                    defaults[col] = 0
            elif col_type in ["path", "simple_path", "complex_path"]:
                # For path types, check the default from attributes
                if col_default_from_attrs is not None:
                    if col_default_from_attrs.lower() == 'none':
                        defaults[col] = {"Path1": "None"}
                    else:
                        defaults[col] = {"Path1": col_default_from_attrs}
                else:
                    defaults[col] = {"Path1": "None"}
            else:
                # String and other types
                if col_default_from_attrs is not None:
                    defaults[col] = col_default_from_attrs
                else:
                    defaults[col] = "None"
        
        return defaults


    def _build_hierarchy_path(self, structured_data: Dict, hierarchy_columns: List[str], 
                            group_keys: tuple, group_df: pd.DataFrame, 
                            property_columns: List[str], property_placement: Dict[str, int],
                            merged_paired_cols: List[str], read_type_column: str = None,
                            property_defaults: Dict[str, Any] = None) -> None:
        """
        Build the hierarchy path for a specific group and place properties at appropriate levels.
        Creates labeled hierarchy levels for columns that represent categorical groupings.
        Ensures all properties are included with appropriate defaults for consistency.
        """
        # Determine which hierarchy columns should have explicit labels
        # Always label columns beyond the first one that have categorical nature
        labeled_columns = set()
        for i, col in enumerate(hierarchy_columns[1:], 1):  # Skip first column
            col_lower = col.lower()
            # Columns that typically represent categorical groupings that should be labeled
            if any(pattern in col_lower for pattern in ['type', 'id', 'name', 'sample', 'library', 'group', 'class', 'category']):
                labeled_columns.add(i)
        
        # Group properties by their placement level
        properties_by_level = {}
        for col in property_columns:
            target_level = property_placement.get(col, len(hierarchy_columns) - 1)
            if target_level not in properties_by_level:
                properties_by_level[target_level] = []
            properties_by_level[target_level].append(col)
        
        # Build hierarchy structure and place properties at each level
        for level in range(len(hierarchy_columns)):
            # Navigate to the current level
            current = structured_data
            for nav_level in range(level + 1):
                hier_col = hierarchy_columns[nav_level]
                hier_val = group_keys[nav_level]
                
                # For the first level, just use the value directly
                if nav_level == 0:
                    if hier_val not in current:
                        current[hier_val] = {}
                    if nav_level == level:
                        level_dict = current[hier_val]
                    else:
                        current = current[hier_val]
                else:
                    # For subsequent levels, check if this should be a labeled level
                    if nav_level in labeled_columns:
                        # Create labeled level structure: column_name: { value: {...} }
                        if hier_col not in current:
                            current[hier_col] = {}
                        if hier_val not in current[hier_col]:
                            current[hier_col][hier_val] = {}
                        if nav_level == level:
                            level_dict = current[hier_col][hier_val]
                        else:
                            current = current[hier_col][hier_val]
                    else:
                        # Use value directly (unlabeled level)
                        if hier_val not in current:
                            current[hier_val] = {}
                        if nav_level == level:
                            level_dict = current[hier_val]
                        else:
                            current = current[hier_val]
            
            # Add ALL properties that belong to this level (including defaults for missing ones)
            if level in properties_by_level:
                for col in properties_by_level[level]:
                    # CRITICAL FIX: Don't filter out empty values prematurely
                    # Get ALL values from the group
                    all_values = group_df[col].tolist()
                    
                    # For single-row groups, handle the value directly
                    if len(all_values) == 1:
                        value = all_values[0]
                        # Check if it's empty/None
                        if pd.isna(value) or str(value).strip() in ['', 'None']:
                            # Use default value for this column
                            if property_defaults and col in property_defaults:
                                level_dict[col] = property_defaults[col]
                            else:
                                level_dict[col] = "None"
                        else:
                            # Process the non-empty value
                            self._add_property_to_level(level_dict, col, [value], 
                                                    merged_paired_cols, read_type_column, group_df)
                    else:
                        # Multiple values - filter for non-empty ones
                        values = [v for v in all_values 
                                if pd.notna(v) and str(v).strip() != 'None' and str(v).strip() != '']
                        
                        # If we have actual values, use them
                        if values:
                            self._add_property_to_level(level_dict, col, values, 
                                                    merged_paired_cols, read_type_column, group_df)
                        else:
                            # All values are empty - use default
                            if property_defaults and col in property_defaults:
                                level_dict[col] = property_defaults[col]
                            else:
                                level_dict[col] = "None"


    def _add_property_to_level(self, target_dict: Dict, col: str, values: List, 
                            merged_paired_cols: List[str], read_type_column: str = None,
                            group_df: pd.DataFrame = None) -> None:
        """
        Add a property to a specific level in the hierarchy.
        """
        # Special handling for merged paired columns
        if col in merged_paired_cols:
            # Process each row's paired data
            path_counter = 1
            paths_dict = {}
            
            for value in values:
                # This is already a comma-separated pair
                paths_dict[f"Path{path_counter}"] = str(value).strip()
                path_counter += 1
            
            if paths_dict:
                target_dict[col] = paths_dict
            else:
                target_dict[col] = {"Path1": "None"}  # FIX: Use proper default for empty path
        
        # If all values are the same, store as a single value
        elif len(set(str(v) for v in values)) == 1:
            value = values[0]
            processed = self._process_value(value, col, 
                                        group_df[read_type_column].iloc[0] if read_type_column and group_df is not None else None)
            target_dict[col] = processed
        
        # Handle multiple different values
        else:
            # CRITICAL FIX: Don't treat categorical/hierarchy columns as file paths
            # Check if this column contains actual file paths (not just has "read" in the name)
            contains_actual_paths = any('/' in str(v) or '\\' in str(v) for v in values)
            is_path_type = self.inferred_column_types.get(col) in ["simple_path", "complex_path"]
            has_file_indicator = any(indicator in col.lower() for indicator in ['file', 'files', 'path', 'paths'])
            
            # Only treat as paths if it actually contains path separators OR is explicitly a file column
            should_treat_as_paths = contains_actual_paths or (is_path_type and has_file_indicator)
            
            if should_treat_as_paths:
                # Create numbered paths for multiple file values
                paths_dict = {}
                for i, value in enumerate(values, 1):
                    paths_dict[f"Path{i}"] = str(value).strip()
                target_dict[col] = paths_dict
            else:
                # For non-path columns with multiple values, this indicates a structural issue
                # This should not happen if hierarchy detection worked correctly
                self.warnings.append(f"Column '{col}' has multiple different values at the same hierarchy level: {values}")
                self.warnings.append(f"This suggests '{col}' should be a hierarchy level, not a property")
                
                # Take the first value as a fallback, but this indicates the hierarchy detection failed
                processed = self._process_value(values[0], col, 
                                            group_df[read_type_column].iloc[0] if read_type_column and group_df is not None else None)
                target_dict[col] = processed


    def _determine_property_placement(self, df: pd.DataFrame, hierarchy_columns: List[str], 
                                    property_columns: List[str]) -> Dict[str, int]:
        """
        Determine at which hierarchy level each property should be placed.
        Enhanced with domain knowledge about genomic data relationships.
        FIXED: Better placement logic based on where properties actually vary.
        """
        placement = {}
        
        # Domain knowledge: which properties belong to which entity types
        # Map properties to their expected hierarchy level
        property_associations = {
            # Assembly-level properties (usually with asm_id)
            'assembly': ['asm_file', 'assembly_file', 'genome_file', 'sex_chr', 
                        'chromosome', 'masking', 'mask', 'reference', 'related_sp', 'related_species'],
            # Sample-level properties (usually with sample_id)  
            'sample': ['sample_name', 'individual', 'tissue', 'population'],
            # Read-level properties (usually with read_type or deepest level)
            'reads': ['read_files', 'reads_path', 'trim_qc', 'read_qc', 'trim', 
                    'adapter', 'quality', 'coverage', 'read_paths']
        }
        
        # Create a mapping of hierarchy column names to their indices
        hierarchy_map = {col: i for i, col in enumerate(hierarchy_columns)}
        
        for prop_col in property_columns:
            prop_lower = prop_col.lower()
            best_level = len(hierarchy_columns) - 1  # Default to deepest
            placed = False
            
            # First, try domain knowledge placement
            # Check assembly-level properties
            if any(term in prop_lower for term in property_associations['assembly']):
                # Find the asm_id level if it exists
                for hier_col in ['asm_id', 'assembly_id', 'assembly', 'genome_id', 'genome']:
                    if hier_col in hierarchy_map:
                        best_level = hierarchy_map[hier_col]
                        placed = True
                        break
                # If no asm_id level found but it's an assembly property, place after first level
                if not placed and len(hierarchy_columns) > 1:
                    best_level = 1
                    placed = True
            
            # Check sample-level properties  
            elif any(term in prop_lower for term in property_associations['sample']):
                # Find the sample_id level if it exists
                for hier_col in ['sample_id', 'sample', 'library_id', 'library', 'individual_id']:
                    if hier_col in hierarchy_map:
                        best_level = hierarchy_map[hier_col]
                        placed = True
                        break
                # If no sample_id level found, place at second-to-last level
                if not placed and len(hierarchy_columns) > 2:
                    best_level = len(hierarchy_columns) - 2
                    placed = True
                    
            # Check read-level properties
            elif any(term in prop_lower for term in property_associations['reads']):
                # These go at the deepest level or with read_type
                for hier_col in ['read_type', 'seq_type', 'data_type', 'sequencing_type']:
                    if hier_col in hierarchy_map:
                        best_level = hierarchy_map[hier_col]
                        placed = True
                        break
                # Default to deepest level for read properties
                if not placed:
                    best_level = len(hierarchy_columns) - 1
                    placed = True
            
            # If not placed by domain knowledge, use data-driven approach
            if not placed:
                # Find where property becomes constant by checking from deepest to shallowest
                for level in range(len(hierarchy_columns)):
                    group_cols = hierarchy_columns[:level + 1]
                    
                    # Check if property is constant within groups at this level
                    is_constant = True
                    for _, group in df.groupby(group_cols, dropna=False):
                        # Get non-null unique values
                        unique_vals = group[prop_col].dropna().unique()
                        if len(unique_vals) > 1:
                            # Property varies within this group, need deeper level
                            is_constant = False
                            break
                    
                    if is_constant:
                        best_level = level
                        break
            
            placement[prop_col] = best_level
        
        return placement


    def _determine_property_level(self, df: pd.DataFrame, property_col: str, hierarchy_columns: List[str]) -> int:
        """
        Determine the level at which a property column should be placed in the hierarchy.
        The property is placed at the level where it first remains constant when grouped by
        hierarchy columns up to that level.
        
        Args:
            df: DataFrame with the data
            property_col: Column name of the property to place
            hierarchy_columns: List of column names that form the hierarchy
            
        Returns:
            The index in hierarchy_columns where this property should be placed
        """
        # Start checking from the first hierarchy level
        for i in range(len(hierarchy_columns)):
            # Group by hierarchy columns up to and including this level
            group_cols = hierarchy_columns[:i+1]
            
            # Check if the property is constant within each group
            constants_in_all_groups = True
            for _, group in df.groupby(group_cols):
                if len(group[property_col].unique()) > 1:
                    constants_in_all_groups = False
                    break
            
            # If constant within all groups at this level, place it here
            if constants_in_all_groups:
                return i
        
        # If not constant at any level, place at the deepest level
        return len(hierarchy_columns) - 1


    def _process_value(self, value, column, read_type=None):
        """
        Process a value based on its column type and content.
        Respects attributes file settings as the primary authority.
        """
        # Check for null-like values
        is_null = (pd.isna(value) or value == '' or value == 'None' or 
                value == 'NA' or value == '-' or value == ' ')
        
        # Get column attributes from the attributes file if available
        # But ignore if --overwrite-attributes flag is set
        overwrite_attrs = self.analysis_results.get("file_info", {}).get("overwrite_attributes", False)
        
        if overwrite_attrs:
            # Act as if no attributes file exists - use smart detection only
            attributes = {}
            col_attrs = {}
            col_type_from_attrs = "unknown"
            col_default_from_attrs = None
        else:
            # Read attributes file normally
            file_path = self.analysis_results.get("file_info", {}).get("path", "")
            custom_attrs_path = self.analysis_results.get("file_info", {}).get("attributes_path", None)
            attributes = self.read_attributes_file(file_path, custom_attrs_path) if file_path else {}
            
            col_attrs = attributes.get(column, {})
            col_type_from_attrs = col_attrs.get("type", "unknown")
            col_default_from_attrs = col_attrs.get("default", None)
        
        # Handle null-like values based on attributes file or inferred type
        if is_null:
            # Use default from attributes file if available
            if col_default_from_attrs is not None:
                if col_default_from_attrs.lower() in ['true', 'false']:
                    return col_default_from_attrs.lower() == 'true'
                elif col_default_from_attrs == 'None':
                    return 'None'
                else:
                    try:
                        if col_type_from_attrs == "integer":
                            return int(col_default_from_attrs)
                        else:
                            return col_default_from_attrs
                    except:
                        return col_default_from_attrs
            
            # Fallback to inferred type defaults
            elif self.inferred_column_types.get(column) == "integer":
                return 0
            else:
                return 'None'
        
        # Check if this is a merged paired column
        if hasattr(self, 'merged_paired_columns') and column in self.merged_paired_columns:
            parsed_paths = self._parse_read_paths(value, read_type)
            if not parsed_paths:
                return {"Path1": "None"}
            return parsed_paths

        # CRITICAL: Special handling for assembly file columns (only for single files, not comma-separated lists)
        if any(name in column.lower() for name in ['asm_file', 'assembly_file', 'reference', 'ref_file', 'genome_file']):
            str_value = str(value).strip()
            # If it's comma-separated, let the normal path logic handle it
            if ',' not in str_value:
                return {"Path1": str_value}
        
        # PRIMARY: Process based on attributes file type (if specified and not "unknown")
        # Attributes file is now AUTHORITATIVE - no smart overrides
        if col_type_from_attrs and col_type_from_attrs != "unknown":
            if col_type_from_attrs == "boolean":
                return self._to_boolean(value)
            elif col_type_from_attrs == "integer":
                try:
                    numeric_val = float(value)
                    if numeric_val.is_integer():
                        return int(numeric_val)
                    else:
                        return numeric_val
                except:
                    if isinstance(value, str) and value.lower() == 'auto':
                        return 'auto'
                    else:
                        return 0
            elif col_type_from_attrs in ["path", "simple_path", "complex_path"]:
                # Handle path types from attributes
                str_value = str(value).strip()
                
                # CRITICAL FIX: Even for "path" or "simple_path", check if it contains multiple paths OR accessions
                # Split by comma to check for multiple values
                if ',' in str_value:
                    values_list = [v.strip() for v in str_value.split(',')]
                    # Check if these are accession numbers OR multiple paths
                    if self._contains_accession_numbers(values_list) or str_value.count('/') > 1 or str_value.count('\\') > 1:
                        # This looks like multiple paths/accessions - parse them
                        parsed_paths = self._parse_read_paths(value, read_type)
                        if not parsed_paths:
                            return {"Path1": "None"}
                        return parsed_paths
                
                # Single value or not accessions/multiple paths
                if col_type_from_attrs == "simple_path" or col_type_from_attrs == "path":
                    return {"Path1": str_value}
                else:  # complex_path
                    parsed_paths = self._parse_read_paths(value, read_type)
                    if not parsed_paths:
                        return {"Path1": "None"}
                    return parsed_paths
            else:
                # For string, category, or other types, return as string
                return str(value).strip()
        
        # FALLBACK: Smart value-level processing for unknown types only
        # Check if this individual value is clearly a boolean
        if self._is_boolean(value) and not ('/' in str(value) or '\\' in str(value)):
            return self._to_boolean(value)
        
        # Check if this individual value is clearly an integer
        if self._is_integer(value) and self.inferred_column_types.get(column) != "simple_path":
            try:
                numeric_val = float(value)
                if numeric_val.is_integer():
                    return int(numeric_val)
                else:
                    return numeric_val
            except:
                pass
        
        # Process based on inferred column type (fallback when no attributes)
        if self.inferred_column_types.get(column) == "boolean":
            return self._to_boolean(value)
        elif self.inferred_column_types.get(column) == "integer":
            try:
                numeric_val = float(value)
                if numeric_val.is_integer():
                    return int(numeric_val)
                else:
                    return numeric_val
            except:
                if isinstance(value, str) and value.lower() == 'auto':
                    return 'auto'
                else:
                    return 0
        elif self.inferred_column_types.get(column) in ["simple_path", "complex_path"]:
            str_value = str(value).strip()
            looks_like_path = (
                '/' in str_value or '\\' in str_value or
                '.' in str_value and len(str_value) > 3 or
                any(ext in str_value.lower() for ext in ['.gz', '.fa', '.fasta', '.fq', '.fastq', '.bam', '.vcf'])
            )
            
            if looks_like_path:
                if self.inferred_column_types.get(column) == "simple_path":
                    return {"Path1": str_value}
                else:
                    parsed_paths = self._parse_read_paths(value, read_type)
                    if not parsed_paths:
                        return {"Path1": "None"}
                    return parsed_paths
            else:
                return str_value
        elif isinstance(value, str) and ('/' in value or '\\' in value) and ('file' in column.lower() or 'path' in column.lower()):
            return {"Path1": str(value).strip()}
        else:
            return str(value).strip()


    def _build_correct_hierarchy(self, df: pd.DataFrame, hierarchy_columns: List[str], 
                                result_dict: Dict, paired_columns: List[str] = None,
                                read_type_column: str = None) -> None:
        """
        Build a hierarchical structure where:
        - The first column is the top level
        - Certain columns are treated as properties rather than nested levels
        - Other columns create true nested levels
        
        Args:
            df: DataFrame with data to process
            hierarchy_columns: List of column names that form the hierarchy
            result_dict: Dictionary to populate with the hierarchical structure
            paired_columns: List of columns containing paired data
            read_type_column: Column containing read type information
        """
        # Get all columns
        all_columns = list(df.columns)
        
        # Determine which columns are actual hierarchy levels vs. properties
        # This is determined by analyzing what values remain constant when grouped by preceding columns
        level_columns = []
        property_columns_by_level = {}
        
        # The first column is always a level
        level_columns.append(hierarchy_columns[0])
        property_columns_by_level[0] = []
        
        # Analyze each hierarchy column to decide if it's a level or a property
        for i in range(1, len(hierarchy_columns)):
            col = hierarchy_columns[i]
            
            # Group by all preceding hierarchy columns
            prev_columns = hierarchy_columns[:i]
            
            # Check if this column has multiple distinct values for any group
            has_multiple_values = False
            for _, group in df.groupby(prev_columns):
                if len(group[col].unique()) > 1:
                    has_multiple_values = True
                    break
            
            # If this column has multiple values within a group, it's a level
            # Otherwise, it's a property of the preceding level
            if has_multiple_values:
                level_columns.append(col)
                property_columns_by_level[i] = []
            else:
                # Determine which preceding level this property belongs to
                for j in range(i-1, -1, -1):
                    if j < len(level_columns):
                        property_columns_by_level[j].append(col)
                        break
        
        # Add non-hierarchy columns as properties at the appropriate level
        for col in all_columns:
            if col not in hierarchy_columns:
                # Determine which level this property belongs to
                added = False
                for i in range(len(level_columns)-1, -1, -1):
                    # Group by hierarchy columns up to this level
                    level_cols = level_columns[:i+1]
                    
                    # Check if this column has multiple values for any group
                    has_multiple_values = False
                    for _, group in df.groupby(level_cols):
                        if len(group[col].unique()) > 1:
                            has_multiple_values = True
                            break
                    
                    if not has_multiple_values:
                        # This property belongs to this level
                        property_columns_by_level[i].append(col)
                        added = True
                        break
                
                # If not added yet, it belongs to the deepest level
                if not added:
                    property_columns_by_level[len(level_columns)-1].append(col)
        
        # Group by all level columns to process each unique combination
        grouped = df.groupby(level_columns)
        
        # Process each group
        for group_keys, group_df in grouped:
            # Ensure group_keys is a tuple
            if not isinstance(group_keys, tuple):
                group_keys = (group_keys,)
                
            # Get the first row for this group (for property values)
            row = group_df.iloc[0]
            
            # Build the hierarchy
            current_dict = result_dict
            for i, (level_col, level_val) in enumerate(zip(level_columns, group_keys)):
                # Create the current level if needed
                if level_val not in current_dict:
                    current_dict[level_val] = {}
                
                # Add properties for this level
                for prop_col in property_columns_by_level[i]:
                    prop_val = row[prop_col]
                    
                    # Skip empty values
                    if pd.isna(prop_val) or prop_val == '' or prop_val == 'None' or prop_val == '-':
                        continue
                    
                    # Process based on column type
                    col_type = self.inferred_column_types.get(prop_col, "unknown")
                    
                    if col_type == "boolean":
                        current_dict[level_val][prop_col] = self._to_boolean(prop_val)
                    elif col_type == "integer":
                        try:
                            numeric_val = float(prop_val)
                            if numeric_val.is_integer():
                                current_dict[level_val][prop_col] = int(numeric_val)
                            else:
                                current_dict[level_val][prop_col] = numeric_val
                        except:
                            if isinstance(prop_val, str) and prop_val.lower() == 'auto':
                                current_dict[level_val][prop_col] = 'auto'
                            else:
                                current_dict[level_val][prop_col] = 0
                    elif col_type == "complex_path":
                        read_type_value = row.get(read_type_column) if read_type_column else None
                        is_paired = prop_col in (paired_columns or [])
                        paths = self._parse_read_paths(prop_val, read_type_value)
                        if paths:
                            current_dict[level_val][prop_col] = paths
                    else:
                        current_dict[level_val][prop_col] = str(prop_val).strip()
                
                # Move to the next level for the next iteration (if not the last level)
                if i < len(level_columns) - 1:
                    current_dict = current_dict[level_val]


    def _detect_hierarchy_columns(self, df: pd.DataFrame) -> List[str]:
        """
        Detect columns that form the hierarchy based on data patterns.
        Enhanced for genomic data with predefined patterns.
        FIXED: Trust domain knowledge for genomic hierarchy patterns.
        """
        columns = list(df.columns)
        
        # Read attributes for hierarchical flags
        attributes = {}
        if not self.analysis_results.get("file_info", {}).get("overwrite_attributes", False):
            file_path = self.analysis_results.get("file_info", {}).get("path", "")
            custom_attrs_path = self.analysis_results.get("file_info", {}).get("attributes_path", None)
            attributes = self.read_attributes_file(file_path, custom_attrs_path) if file_path else {}
        
        # Columns that should never be hierarchy levels
        non_hierarchy_patterns = ['_x', '_Rx', '_rx', '.Rx', '.rx', 'Paired_x', 'Paired_Rx', 'Pairedx']
        
        # Start with the first column as root
        hierarchy_columns = [columns[0]]
        
        # CRITICAL: Define genomic hierarchy patterns with TRUST LEVELS
        # These patterns represent domain knowledge about genomic data structure
        genomic_hierarchy_patterns = {
            # High trust - almost always hierarchical in genomic workflows
            'high_trust': {
                'asm_id': ['asm_id', 'assembly_id', 'assembly', 'genome_id', 'genome'],
                'sample_id': ['sample_id', 'sample', 'library_id', 'library', 'individual_id', 'individual'],
                'read_type': ['read_type', 'seq_type', 'sequencing_type', 'platform']
                # NOTE: Removed 'data_type' - it's often 1:1 with sample
            },
            # Medium trust - often hierarchical but check data
            'medium_trust': {
                'data_type': ['data_type'],  # Moved here to check with data
                'method': ['method', 'approach', 'technique'],
                'replicate': ['replicate', 'rep', 'replicate_id']
            }
        }
        
        # Check each remaining column
        for i in range(1, len(columns)):
            current_col = columns[i]
            
            # Skip merged paired columns
            if any(pattern in current_col for pattern in non_hierarchy_patterns):
                continue
            
            # Check if forced to be property by attributes
            if current_col in attributes and attributes[current_col].get("property", False):
                continue  # Skip adding to hierarchy
            
            # Check if forced hierarchical by attributes
            if current_col in attributes and attributes[current_col].get("hierarchical", False):
                hierarchy_columns.append(current_col)
                continue
            
            # Check HIGH TRUST genomic patterns first - these ALWAYS become hierarchical
            is_high_trust_genomic = False
            for pattern_type, patterns in genomic_hierarchy_patterns['high_trust'].items():
                if current_col.lower() in [p.lower() for p in patterns]:
                    is_high_trust_genomic = True
                    # For high-trust patterns, ALWAYS add them - no 1:1 check
                    # asm_id should come before sample_id, sample_id before read_type
                    if pattern_type == 'asm_id' and current_col not in hierarchy_columns:
                        hierarchy_columns.append(current_col)
                    elif pattern_type == 'sample_id' and current_col not in hierarchy_columns:
                        hierarchy_columns.append(current_col)
                    elif pattern_type == 'read_type' and current_col not in hierarchy_columns:
                        hierarchy_columns.append(current_col)
                    break
            
            if is_high_trust_genomic:
                continue  # Already added, move to next column
            
            # Check MEDIUM TRUST patterns - verify with data including 1:1 check
            is_medium_trust_genomic = False
            for pattern_type, patterns in genomic_hierarchy_patterns.get('medium_trust', {}).items():
                if current_col.lower() in [p.lower() for p in patterns]:
                    is_medium_trust_genomic = True
                    
                    # For medium trust, check if it's 1:1 with existing hierarchy
                    if len(hierarchy_columns) > 0:
                        is_one_to_one = True
                        for _, group in df.groupby(hierarchy_columns, dropna=False):
                            unique_vals = group[current_col].dropna().unique()
                            if len(unique_vals) > 1:
                                is_one_to_one = False
                                break
                        
                        # If it's 1:1, don't make it hierarchical
                        if is_one_to_one:
                            break
                    
                    # Also check if it creates subdivisions
                    is_hierarchical = self._is_hierarchical_column(df, hierarchy_columns, current_col)
                    if is_hierarchical:
                        hierarchy_columns.append(current_col)
                    break
            
            if is_medium_trust_genomic:
                continue
            
            # SPECIAL CHECK for columns in the second position (only for non-genomic columns)
            if i == 1 and not is_high_trust_genomic and not is_medium_trust_genomic:
                # Check if this is a 1:1 mapping with the first column
                first_col = columns[0]
                
                # Group by first column and check if each group has only one unique value in current column
                is_one_to_one = True
                for _, group in df.groupby(first_col, dropna=False):
                    unique_vals = group[current_col].dropna().unique()
                    if len(unique_vals) > 1:
                        is_one_to_one = False
                        break
                
                # If it's 1:1, it should be a property, not hierarchical
                if is_one_to_one:
                    continue  # Skip adding to hierarchy
            
            # Otherwise use smart detection for non-genomic columns
            is_hierarchical = self._is_hierarchical_column(df, hierarchy_columns, current_col)
            
            if is_hierarchical:
                hierarchy_columns.append(current_col)
        
        return hierarchy_columns


    def _should_be_hierarchical(self, df: pd.DataFrame, existing_hierarchy: List[str], 
                            test_col: str, remaining_cols: List[str]) -> bool:
        """
        Improved logic to determine if a column should be part of the hierarchy.
        More conservative approach - only hierarchical if it actually creates subdivisions.
        """
        # Get column content analysis
        unique_values = df[test_col].dropna().unique()
        total_rows = len(df)
        non_empty_rows = len(df[test_col].dropna())
        unique_count = len(unique_values)
        
        # Skip if all values are paths or files
        if any('/' in str(v) or '\\' in str(v) for v in unique_values):
            return False
        
        # Skip boolean-like columns
        str_values = [str(v).lower().strip() for v in unique_values if pd.notna(v)]
        boolean_like = {'true', 'false', 'yes', 'no', 'y', 'n', '1', '0'}
        if all(v in boolean_like for v in str_values) and unique_count <= 3:
            return False
        
        col_lower = test_col.lower()
        
        # CRITICAL: First check if this column actually creates subdivisions
        # This is the MOST IMPORTANT check
        if len(existing_hierarchy) > 0:
            creates_any_subdivision = False
            
            # Group by ALL existing hierarchy levels
            grouped = df.groupby(existing_hierarchy, dropna=False)
            
            for group_key, group in grouped:
                # Get unique values in this group for the test column
                group_values = group[test_col].dropna()
                if len(group_values) > 0:
                    unique_in_group = group_values.nunique()
                    # If this group has more than one unique value, it creates subdivisions
                    if unique_in_group > 1:
                        creates_any_subdivision = True
                        break
            
            # If NO group has multiple values, this column should NOT be hierarchical
            if not creates_any_subdivision:
                return False
        
        # Special handling for columns with "type" in the name
        # They're ONLY hierarchical if they create actual subdivisions
        if 'type' in col_lower:
            # We already checked for subdivisions above
            # If we're here, it means it DOES create subdivisions
            # But let's be extra careful - require strong evidence
            if len(existing_hierarchy) > 0:
                # Count how many groups have subdivisions
                groups_with_subdivisions = 0
                total_groups = 0
                
                for _, group in df.groupby(existing_hierarchy, dropna=False):
                    total_groups += 1
                    if group[test_col].dropna().nunique() > 1:
                        groups_with_subdivisions += 1
                
                # Require at least 20% of groups to have subdivisions for "type" columns
                if total_groups > 0 and groups_with_subdivisions / total_groups < 0.2:
                    return False
                    
            return True  # Has subdivisions and enough evidence
        
        # For columns matching genomic hierarchy patterns
        # These are more likely to be hierarchical even with limited data
        genomic_patterns = ['sample_id', 'library_id', 'asm_id', 'assembly_id']
        if any(pattern in col_lower for pattern in genomic_patterns):
            # For these special columns, we're more permissive
            # They can be hierarchical even if current data doesn't show subdivisions
            # (because we know from domain knowledge they often are)
            if non_empty_rows > 0:
                return True
        
        # For other columns, only hierarchical if they have reasonable cardinality
        # AND create subdivisions (already checked above)
        if 2 <= unique_count <= min(20, total_rows * 0.5):
            return True
        
        return False


    def _is_hierarchical_column(self, df: pd.DataFrame, existing_hierarchy: List[str], test_col: str) -> bool:
        """
        Determine if a column should be part of the hierarchy based on data patterns.
        FIXED: More strict about requiring actual subdivisions.
        """
        # Get unique values and analyze content
        unique_values = df[test_col].dropna().unique()
        total_rows = len(df[test_col].dropna())
        unique_count = len(unique_values)
        
        # Rule 1: Skip columns that look like file paths or URLs
        if any('/' in str(v) or '\\' in str(v) for v in unique_values):
            return False
        
        # Rule 2: Skip columns that look like boolean flags
        str_values = [str(v).lower().strip() for v in unique_values if pd.notna(v)]
        boolean_like_values = {'true', 'false', 'yes', 'no', 'y', 'n', '1', '0', 'none', 'na', 'false', 'true'}
        if all(v in boolean_like_values for v in str_values) and unique_count <= 3:
            return False
        
        # Rule 3: Skip columns with very high cardinality
        if total_rows <= 5:
            max_cardinality_threshold = total_rows
        else:
            max_cardinality_threshold = total_rows * 0.7

        if unique_count > max_cardinality_threshold:
            return False
        
        col_lower = test_col.lower()
        
        # CRITICAL FIX: For ANY column (including those with "type" in the name),
        # REQUIRE that it creates actual subdivisions
        if len(existing_hierarchy) > 0:
            has_any_subdivision = False
            groups_with_multiple_rows = 0
            total_groups = 0
            
            # Group by existing hierarchy and check if this column varies within ANY group
            for group_key, group in df.groupby(existing_hierarchy, dropna=False):
                total_groups += 1
                
                # IMPORTANT: Only consider groups that have more than one row
                # A column can only create subdivisions if the group has multiple rows to subdivide!
                if len(group) > 1:
                    groups_with_multiple_rows += 1
                    group_values = group[test_col].dropna().unique()
                    if len(group_values) > 1:
                        has_any_subdivision = True
                        # Don't break - we want to count all groups for better decisions
            
            # NEW LOGIC: If NO groups have multiple rows, this column CAN'T be hierarchical
            # because there's nothing to subdivide!
            if groups_with_multiple_rows == 0:
                # Special case: each parent value has exactly one child value
                # This means it's a 1:1 mapping and should be a property, not hierarchy
                return False
            
            # If NO group has multiple values, it's NOT hierarchical
            if not has_any_subdivision:
                return False
            
            # Additional check for "type" columns - require significant evidence
            if 'type' in col_lower:
                # Require at least 20% of multi-row groups to have subdivisions
                if groups_with_multiple_rows > 0:
                    subdivision_ratio = 1.0 if has_any_subdivision else 0.0
                    # For type columns, we want clear evidence of subdivision
                    # If only a few groups have subdivisions, it's likely not hierarchical
                    if subdivision_ratio < 0.2:
                        return False
        
        # Only after confirming it creates subdivisions, check other criteria
        
        # Special handling for "type" columns - now they need subdivisions too
        if 'type' in col_lower:
            # We already verified it has subdivisions above
            return True
        
        # Special handling for "level" columns
        if 'level' in col_lower:
            # Already verified subdivisions
            return False  # Usually level describes a property
        
        # Strong hierarchy indicators
        strong_hierarchical_hints = ['id', 'name', 'sample', 'library']
        if any(hint in col_lower for hint in strong_hierarchical_hints):
            # Already verified it creates subdivisions
            return True
        
        # Property column patterns
        property_hints = [
            'file', 'files', 'path', 'paths', 'reads', 'qc', 'trim', 'mask', 'masking',
            'size', 'length', 'count', 'score', 'rate', 'lineage', 'asm', 'db',
            'chr', 'chromosome', 'scaffold', 'contig'
        ]
        
        if any(hint in col_lower for hint in property_hints):
            return False
        
        # If we get here, it passed the subdivision test and other checks
        return True


    def _sanitize_hierarchical_values(self, df: pd.DataFrame, hierarchy_columns: List[str], 
                                    attributes: Dict[str, Dict[str, Any]], 
                                    sanitise_fields_flag: bool = False) -> pd.DataFrame:
        """
        Sanitize hierarchical field values by replacing problematic characters with underscores.
        Only affects columns that will become YAML keys (hierarchy columns).
        
        Args:
            df: DataFrame to sanitize
            hierarchy_columns: List of columns that form the hierarchy
            attributes: Attributes dictionary to check for 'sanitise' flag
            sanitise_fields_flag: If True, sanitize ALL hierarchy columns (from --sanitise-fields flag)
            
        Returns:
            DataFrame with sanitized hierarchical values
        """
        df_sanitized = df.copy()
        sanitized_columns = []
        
        for col in hierarchy_columns:
            should_sanitize = False
            
            # Check if --sanitise-fields flag was used (sanitize all hierarchy columns)
            if sanitise_fields_flag:
                should_sanitize = True
            
            # Always check attributes file for individual column sanitise flags
            elif col in attributes and attributes[col].get('sanitise', False):
                should_sanitize = True
            
            if should_sanitize:
                # Sanitize values in this column
                df_sanitized[col] = df_sanitized[col].apply(self._sanitize_value)
                sanitized_columns.append(col)
        
        if sanitized_columns:
            if sanitise_fields_flag:
                self.warnings.append(f"--sanitise-fields flag: Sanitized all hierarchical values in columns: {', '.join(sanitized_columns)}")
            else:
                self.warnings.append(f"Sanitized hierarchical values based on .attributes file in columns: {', '.join(sanitized_columns)}")
        
        return df_sanitized


    def _sanitize_value(self, value) -> str:
        """
        Sanitize a single value by replacing problematic characters with underscores.
        
        Args:
            value: Value to sanitize
            
        Returns:
            Sanitized string value
        """
        if pd.isna(value) or value is None:
            return value
        
        # Convert to string
        str_value = str(value)
        
        # Replace problematic characters with underscores
        # Problematic chars: spaces, quotes, brackets, special chars that cause YAML issues
        import re
        sanitized = re.sub(r'[^\w\-.]', '_', str_value)
        
        # Clean up multiple consecutive underscores
        sanitized = re.sub(r'_+', '_', sanitized)
        
        # Remove leading/trailing underscores
        sanitized = sanitized.strip('_')
        
        return sanitized


    def _is_path_column(self, column_name: str) -> bool:
        """Check if column name suggests it contains file paths"""
        indicators = ['path', 'file', 'read', 'asm', 'reference', 'ref']
        return any(ind in column_name.lower() for ind in indicators)


    def _contains_accession_numbers(self, paths_list: list) -> bool:
        """
        Check if a list of paths contains accession numbers rather than actual file paths.
        
        Accession patterns:
        - SRR* (SRA/ENA run accessions)
        - ERR* (ENA run accessions)
        - DRR* (DDBJ run accessions)
        - GCA_* (GenBank assembly accessions)
        - GCF_* (RefSeq assembly accessions)
        
        Args:
            paths_list: List of path strings to check
            
        Returns:
            True if list contains accession numbers, False otherwise
        """
        import re
        
        # Accession patterns
        accession_patterns = [
            r'^SRR\d+',     # SRA run: SRR123456
            r'^ERR\d+',     # ENA run: ERR123456
            r'^DRR\d+',     # DDBJ run: DRR123456
            r'^GCA_\d+',    # GenBank assembly: GCA_123456789.1
            r'^GCF_\d+',    # RefSeq assembly: GCF_123456789.1
        ]
        
        # Check if any path matches accession patterns
        for path in paths_list:
            path = path.strip()
            for pattern in accession_patterns:
                if re.match(pattern, path):
                    return True
        
        return False


    def _expand_glob_patterns(self, df: pd.DataFrame, base_dir: str) -> pd.DataFrame:
        """
        Expand glob patterns in path columns to actual file paths.
        
        Args:
            df: DataFrame with potential glob patterns
            base_dir: Base directory for resolving relative paths
        
        Returns:
            DataFrame with expanded file paths
        """
        import glob
        
        df_expanded = df.copy()
        expanded_count = 0
        patterns_expanded = []
        
        for column in df.columns:
            # Only process columns that likely contain paths
            if not self._is_path_column(column):
                continue
            
            for idx in df.index:
                value = df.loc[idx, column]
                
                # Skip empty or non-string values
                if pd.isna(value) or not isinstance(value, str):
                    continue
                
                # Check if contains glob patterns
                if any(char in value for char in ['*', '?', '[']):
                    expanded = self._expand_single_glob_value(value, base_dir)
                    if expanded != value:  # Something changed
                        df_expanded.loc[idx, column] = expanded
                        expanded_count += 1
                        patterns_expanded.append(f"{column}[row {idx}]: {value[:60]}...")
        
        if expanded_count > 0:
            self.warnings.append(f"✓ Expanded {expanded_count} glob pattern(s) to file paths")
        
        return df_expanded


    def _expand_single_glob_value(self, value: str, base_dir: str) -> str:
        """
        Expand a single cell value that may contain one or more glob patterns.
        
        Args:
            value: String potentially containing glob patterns (e.g., "*.fq" or "path1/*.fq,path2/*.fq")
            base_dir: Base directory for resolving relative paths
        
        Returns:
            Comma-separated string of expanded file paths
        """
        import glob
        
        # Handle quoted values (CSV may have quotes around comma-separated paths)
        value = value.strip('"').strip("'")
        
        # Split by comma - might have multiple patterns
        patterns = [p.strip() for p in value.split(',')]
        
        all_expanded = []
        
        for pattern in patterns:
            # Skip if this doesn't look like a glob pattern
            if not any(char in pattern for char in ['*', '?', '[']):
                # Not a glob, keep as-is
                all_expanded.append(pattern)
                continue
            
            # Resolve relative paths
            if not os.path.isabs(pattern):
                full_pattern = os.path.join(base_dir, pattern)
            else:
                full_pattern = pattern
            
            # Expand glob
            matches = sorted(glob.glob(full_pattern))  # Sort for consistency!
            
            if not matches:
                # No matches - warn and keep original
                self.warnings.append(f"⚠ Glob pattern matched no files: {pattern}")
                all_expanded.append(pattern)
            else:
                # Add all matches
                all_expanded.extend(matches)
        
        # Return comma-separated (no space after comma to match NomNom's format)
        return ','.join(all_expanded)


    def _validate_paths_in_structured_data(self, structured_data: dict, base_dir: str) -> bool:
        """
        Validate that all file paths in the structured data exist.
        
        Args:
            structured_data: The hierarchical structured data to validate
            base_dir: Base directory for resolving relative paths (directory of input file)
        
        Returns:
            bool: True if all paths exist, False if any are missing
        """
        missing_paths = []
        checked_paths = 0
        
        def resolve_path(path: str) -> str:
            """Resolve a path (absolute or relative to base_dir)"""
            path = path.strip()
            if os.path.isabs(path):
                return path
            else:
                return os.path.join(base_dir, path)
        
        def check_single_path(path: str, context: str) -> bool:
            """Check if a single path exists and record if missing"""
            if not path or path.strip() == '':
                return True
            
            nonlocal checked_paths
            checked_paths += 1
            
            resolved_path = resolve_path(path)
            
            if not os.path.exists(resolved_path):
                missing_paths.append({
                    'path': path,
                    'resolved': resolved_path,
                    'context': context
                })
                return False
            return True
        
        def traverse_and_validate(data, context_path=''):
            """Recursively traverse the structure and validate paths"""
            if isinstance(data, dict):
                for key, value in data.items():
                    new_context = f"{context_path}/{key}" if context_path else key
                    
                    # Check if this looks like a path field
                    if any(indicator in key.lower() for indicator in ['path', 'file', 'asm_', 'read_']):
                        if isinstance(value, str):
                            # Handle different path formats
                            if ',' in value:
                                # Comma-separated paired reads: "file1.fq, file2.fq"
                                for path in value.split(','):
                                    check_single_path(path.strip(), new_context)
                            
                            elif ';' in value:
                                # Semicolon-separated multi-library: "Path1: lib1.fq; Path2: lib2.fq"
                                for lib_entry in value.split(';'):
                                    if ':' in lib_entry:
                                        _, paths = lib_entry.split(':', 1)
                                        for path in paths.split(','):
                                            check_single_path(path.strip(), new_context)
                                    else:
                                        check_single_path(lib_entry.strip(), new_context)
                            else:
                                # Single path
                                check_single_path(value, new_context)
                        else:
                            # Recurse into nested structure
                            traverse_and_validate(value, new_context)
                    else:
                        # Not a path field, continue traversing
                        traverse_and_validate(value, new_context)
            
            elif isinstance(data, list):
                for item in data:
                    traverse_and_validate(item, context_path)
        
        # Start validation
        traverse_and_validate(structured_data)
        
        # Report results
        if missing_paths:
            self.errors.append(f"Path validation failed: {len(missing_paths)} of {checked_paths} paths not found")
            
            # Add detailed information about missing paths
            for item in missing_paths[:15]:  # Show first 15 missing paths
                self.errors.append(f"  ✗ Missing: {item['path']}")
                self.errors.append(f"    Context: {item['context']}")
                # Only show "Resolved to:" if it's different from the original path (i.e., relative path)
                if item['resolved'] != item['path']:
                    self.errors.append(f"    Resolved to: {item['resolved']}")
            
            if len(missing_paths) > 15:
                self.errors.append(f"  ... and {len(missing_paths) - 15} more missing paths")
            
            return False
        
        # All paths validated successfully
        self.warnings.append(f"✓ Path validation passed: all {checked_paths} file paths exist")
        return True


    def _parse_paired_reads(self, paths_str, is_paired_data=False, column_name=""):
        """
        Generalized parser for path values that handles both paired and unpaired paths.
        
        Args:
            paths_str: String containing one or more file paths
            is_paired_data: Flag indicating if data should be treated as paired-end
            column_name: Name of the column containing these paths
                
        Returns:
            Dictionary with path entries using consistent Path1, Path2, etc. naming
        """
        # Handle None, NaN, empty strings
        if paths_str is None or pd.isna(paths_str) or paths_str == '' or paths_str == 'None':
            return {}
        
        # Ensure we're working with a string
        paths_str = str(paths_str).strip()
        
        # Split the string into individual paths
        if ',' in paths_str:
            paths = [p.strip() for p in paths_str.split(',')]
        else:
            paths = [paths_str.strip()]
        
        # Remove any empty paths
        paths = [p for p in paths if p and p != 'None']
        
        # If no paths, return empty dict
        if not paths:
            return {}
        
        # Always use "Path" as the prefix for consistency
        prefix = "Path"
        
        # SPECIAL CASE: Check if these are accession numbers (ERR, SRR, DRR, GCA_, GCF_)
        # If they are, treat each one as a separate path (no pairing)
        if self._contains_accession_numbers(paths):
            # Each accession gets its own Path entry
            return {f"{prefix}{i+1}": path for i, path in enumerate(paths)}
        
        # If odd number of paths and we're expecting pairs, warn but proceed
        if is_paired_data and len(paths) % 2 != 0:
            self.warnings.append(f"Found odd number of paths ({len(paths)}), unable to treat as paired-end data: {paths_str[:50]}...")
        
        # Check if files contain standard paired-end patterns
        has_r1_r2 = any('_R1' in p or '_R2' in p for p in paths)
        has_1_2_pattern = any(('_1.' in p or '_2.' in p or '.1.' in p or '.2.' in p) for p in paths)
        
        # If we have standard patterns or we're forced to treat as paired
        if (has_r1_r2 or has_1_2_pattern or is_paired_data) and len(paths) % 2 == 0:
            # Get matched pairs
            pairs, remaining = self._match_read_pairs(paths)
            
            # Create output dictionary
            result = {}
            
            # Add all pairs
            for i, (file1, file2) in enumerate(pairs):
                result[f"{prefix}{i+1}"] = f"{file1}, {file2}"
            
            # Add any remaining unpaired files
            for i, file in enumerate(remaining):
                result[f"{prefix}{len(pairs)+i+1}"] = file
                
            return result
        
        # Return as individual entries if we can't confidently pair them
        return {f"{prefix}{i+1}": path for i, path in enumerate(paths)}


    def _match_read_pairs(self, paths):
        """
        Match read pairs using common patterns in genomic data.
        
        Args:
            paths: List of file paths to analyze
            
        Returns:
            Tuple of (matched_pairs, remaining_unpaired_files)
            - matched_pairs: List of tuples, each containing a pair of files (file1, file2)
            - remaining_unpaired_files: List of files that couldn't be paired
        """
        # Common patterns for paired reads
        patterns = [
            ('_R1', '_R2'),
            ('_1.', '_2.'),
            ('.1.', '.2.'),
            ('_1_', '_2_'),
            ('_r1', '_r2'),
            ('_r1.', '_r2.'),
            ('_READ1', '_READ2'),
            ('_read1', '_read2'),
            ('.R1.', '.R2.'),
            ('_1r', '_2r'),
        ]
        
        # Make a copy to avoid modifying the original
        remaining = list(paths)
        pairs = []
        
        # Try exact patterns first
        for pattern1, pattern2 in patterns:
            if not remaining:
                break
                
            # Find all files with pattern1
            p1_files = [p for p in remaining if pattern1 in p]
            
            for p1_file in p1_files:
                if p1_file not in remaining:
                    continue
                    
                # Construct the expected pattern2 filename
                p2_file = p1_file.replace(pattern1, pattern2)
                
                # If the pattern2 file exists, we have a pair
                if p2_file in remaining:
                    pairs.append((p1_file, p2_file))
                    remaining.remove(p1_file)
                    remaining.remove(p2_file)
        
        # If we still have unpaired files, try positional matching
        if remaining and len(remaining) % 2 == 0:
            # Sort to help match similar filenames
            remaining.sort()
            
            # Group by prefix similarities
            i = 0
            while i < len(remaining) - 1:
                file1 = remaining[i]
                file2 = remaining[i+1]
                
                # Check if they differ only by 1→2 in some position
                has_1_2_diff = False
                for pos, (c1, c2) in enumerate(zip(file1, file2)):
                    if c1 == '1' and c2 == '2':
                        # Look at surrounding characters to confirm
                        context_size = 3
                        start = max(0, pos - context_size)
                        end = min(min(len(file1), len(file2)), pos + context_size + 1)
                        
                        # If surrounding context matches, this is likely a pair
                        if file1[start:pos] == file2[start:pos] and file1[pos+1:end] == file2[pos+1:end]:
                            has_1_2_diff = True
                            break
                
                if has_1_2_diff:
                    pairs.append((file1, file2))
                    remaining.remove(file1)
                    remaining.remove(file2)
                    # Don't increment i since we've modified the list
                else:
                    i += 1  # Check next file
        
        # Return the pairs and remaining unpaired files
        return pairs, remaining


    def to_yaml(self) -> str:
        """Convert analysis results to YAML format using the correct top-level key"""
        # Clean up the structured data to ensure proper boolean representation
        cleaned_data = self._clean_for_yaml(self.structured_data)
        
        # Get the correct top-level key from the first column name if available
        top_level_key = "Sample_ID"  # Default fallback
        
        # If we have the current DataFrame available, use its first column name
        if self.current_df is not None and not self.current_df.empty:
            top_level_key = self.current_df.columns[0]
        
        # Only include structured data in the YAML output
        output_data = {
            top_level_key: cleaned_data
        }
        
        # Create a custom YAML dumper that doesn't use anchors/aliases
        class NoAliasDumper(yaml.SafeDumper):
            def ignore_aliases(self, data):
                return True
        
        # Use the custom dumper to prevent anchors/references
        return yaml.dump(output_data, 
                        Dumper=NoAliasDumper,
                        sort_keys=False, 
                        default_flow_style=False, 
                        allow_unicode=True)


    def _clean_for_yaml(self, data):
        """
        Clean the data structure for YAML serialization.
        This version is completely general and doesn't assume specific column names.
        """
        if isinstance(data, dict):
            result = {}
            
            # Get attributes to check column types
            attributes = {}
            if not self.analysis_results.get("file_info", {}).get("overwrite_attributes", False):
                file_path = self.analysis_results.get("file_info", {}).get("path", "")
                custom_attrs_path = self.analysis_results.get("file_info", {}).get("attributes_path", None)
                attributes = self.read_attributes_file(file_path, custom_attrs_path) if file_path else {}
            
            # Get hierarchy columns to avoid converting them to strings
            # Reconstruct hierarchy columns from the structured data
            hierarchy_columns = set()
            if hasattr(self, 'current_df') and self.current_df is not None:
                df_clean = self.current_df.copy()
                for col in df_clean.columns:
                    if df_clean[col].dtype == 'object':
                        df_clean[col] = df_clean[col].astype(str)
                hierarchy_columns = set(self._detect_hierarchy_columns(df_clean))
            
            for key, value in data.items():
                # Clean up whitespace in keys and ensure string keys
                clean_key = str(key).strip() if isinstance(key, str) else str(key)
                
                # CRITICAL: Never convert hierarchy columns to strings
                if clean_key in hierarchy_columns:
                    # This is a hierarchy level, process its nested structure
                    result[clean_key] = self._clean_for_yaml(value)
                    continue
                
                # Check if attributes file specifies this column type
                col_attrs = attributes.get(clean_key, {})
                col_type_from_attrs = col_attrs.get("type", "unknown")
                
                # For non-hierarchy columns, respect attributes file type settings
                if col_type_from_attrs == "string" and not isinstance(value, dict):
                    # Simple value that should be a string
                    result[clean_key] = str(value).strip() if value is not None else 'None'
                    continue
                elif col_type_from_attrs == "string" and isinstance(value, dict) and 'Path1' in value:
                    # This was incorrectly wrapped as a path, unwrap it
                    result[clean_key] = str(value['Path1'])
                    continue
                
                # Handle different value types
                if isinstance(value, (list, dict)):
                    # Special case for simple-path columns incorrectly parsed as complex
                    if isinstance(value, dict) and len(value) == 1 and 'Path1' in value and ',' not in str(value['Path1']):
                        # Check if this should be kept as path or converted to string
                        if col_type_from_attrs == "string":
                            result[clean_key] = str(value['Path1'])
                        else:
                            # Create a completely new dict to avoid references
                            result[clean_key] = {"Path1": str(value['Path1'])}
                    # Fix for integer/float values incorrectly nested
                    elif isinstance(value, dict) and 'Path1' in value:
                        # Check if this could be a numeric value stored as a path
                        try:
                            # Try to convert to numeric
                            numeric_val = float(value['Path1'])
                            # Store as integer if it's a whole number
                            if numeric_val.is_integer():
                                result[clean_key] = int(numeric_val)
                            else:
                                result[clean_key] = float(numeric_val)
                        except:
                            # Not numeric, use type-specific cleanup
                            if str(value['Path1']).strip() == 'auto':
                                result[clean_key] = 'auto'
                            elif col_type_from_attrs == "string":
                                result[clean_key] = str(value['Path1'])
                            else:
                                # Create a deep copy to avoid YAML references
                                result[clean_key] = self._deep_copy_dict(value)
                    else:
                        result[clean_key] = self._clean_for_yaml(value)
                
                # Type-specific value handling
                elif isinstance(value, (int, float)):
                    # Ensure integers are stored as integers, not floats
                    if isinstance(value, float) and value.is_integer():
                        result[clean_key] = int(value)
                    else:
                        result[clean_key] = value
                elif isinstance(value, str):
                    # Clean up string values and create new string objects
                    clean_value = str(value).strip()
                    if clean_value == 'None' or clean_value == '-':
                        result[clean_key] = 'None'
                    # Special handling for path-like strings ONLY if not explicitly string type
                    elif ('/' in clean_value or '\\' in clean_value) and 'file' in clean_key.lower():
                        if col_type_from_attrs != "string":
                            result[clean_key] = {"Path1": str(clean_value)}
                        else:
                            result[clean_key] = str(clean_value)
                    else:
                        # Check if this is a numeric string that should be a number
                        try:
                            numeric_val = float(clean_value)
                            # If it's a whole number, store as integer, otherwise as float
                            if numeric_val.is_integer():
                                result[clean_key] = int(numeric_val)
                            else:
                                result[clean_key] = float(numeric_val)
                        except:
                            # Not numeric, keep as string
                            result[clean_key] = str(clean_value)
                else:
                    result[clean_key] = value
                    
            return result
        elif isinstance(data, list):
            return [self._clean_for_yaml(item) for item in data]
        else:
            return data


    def _deep_copy_dict(self, data):
        """
        Create a deep copy of dictionary to avoid YAML references/anchors.
        """
        if isinstance(data, dict):
            # Create a completely new dict with copied values
            return {str(k): self._deep_copy_dict(v) for k, v in data.items()}
        elif isinstance(data, list):
            return [self._deep_copy_dict(item) for item in data]
        else:
            # For primitive values, return a copy
            if isinstance(data, str):
                return str(data)  # Create new string object
            return data
        

    def save_yaml(self, output_path: str) -> None:
        """Save analysis results to a YAML file"""
        with open(output_path, 'w') as f:
            f.write(self.to_yaml())
    

    def read_attributes_file(self, file_path: str, custom_attributes_path: str = None) -> Dict[str, Dict[str, Any]]:
        """
        Read an attributes file using the simplified format:
        {column name}: {type of data} (missing = {default value}); {mandatory}
        
        Args:
            file_path: Path to the table file (used to derive default attributes path)
            custom_attributes_path: Optional custom path to an attributes file
        """
        # Determine which attributes file to use
        attributes_file = custom_attributes_path
        if not attributes_file:
            # Generate the default attributes file path
            base_path = os.path.splitext(file_path)[0]
            attributes_file = f"{base_path}.attributes"
        
        # Ultra reliable file check - try to actually open it
        try:
            with open(attributes_file, 'r') as test_file:
                pass  # Just testing if we can open it
        except (IOError, FileNotFoundError):
            return {}
        
        attributes = {}
        
        # Read the attributes file
        with open(attributes_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Parse the line using the simplified format
                column_attrs = self._parse_attribute_line(line)
                if column_attrs:
                    attributes.update(column_attrs)
        
        return attributes


    def _parse_attribute_line(self, line: str) -> Dict[str, Dict[str, Any]]:
        """Parse an attribute line in the simplified format."""
        # Split column name and attributes
        parts = line.split(':', 1)
        if len(parts) != 2:
            return {}
            
        column = parts[0].strip()
        attrs_str = parts[1].strip()
        
        # Initialize column attributes
        attrs = {
            "type": "unknown",
            "mandatory": False,
            "default": None,
            "values": [],
            "sanitise": False,
            "hierarchical": False  # NEW: Add hierarchical flag
        }
        
        # Check for flags (mandatory, sanitise, hierarchical)
        flags = []
        if ';' in attrs_str:
            main_part, *flag_parts = attrs_str.split(';')
            attrs_str = main_part.strip()
            flags = [flag.strip() for flag in flag_parts]

        # Process flags
        for flag in flags:
            if flag == 'mandatory':
                attrs["mandatory"] = True
            elif flag == 'sanitise':
                attrs["sanitise"] = True
            elif flag == 'hierarchical':  # NEW
                attrs["hierarchical"] = True

        # Parse type (everything before the first parenthesis)
        if '(' in attrs_str:
            type_name = attrs_str.split('(', 1)[0].strip()
            attrs["type"] = type_name
        else:
            attrs["type"] = attrs_str
        
        # Parse default value if present
        if '(missing = ' in attrs_str:
            default_part = attrs_str.split('(missing = ', 1)[1]
            if ')' in default_part:
                default_value = default_part.split(')', 1)[0].strip()
                attrs["default"] = default_value
        
        return {column: attrs}
 

    def generate_attributes_file(self, file_path: str, df: pd.DataFrame, 
                        custom_attributes_path: str = None,
                        sanitise_requested: bool = False,
                        overwrite_attributes: bool = False) -> str:  # ← NEW

        """
        Generate an attributes file based on the dataframe using the simplified format.
        If an existing attributes file exists, preserve its settings while updating with new information.
        
        Args:
            file_path: Path to the table file
            df: DataFrame to analyze
            custom_attributes_path: Optional path to a custom attributes file
            sanitise_requested: Whether --sanitise-fields was requested
            
        Returns:
            Path to the created attributes file
        """
        # Determine the output attributes file path
        output_attributes_path = custom_attributes_path
        if not output_attributes_path:
            # Default: use input filename with .attributes extension
            base_path = os.path.splitext(file_path)[0]
            output_attributes_path = f"{base_path}.attributes"
        
        # Read existing attributes if available
        existing_attributes = {}
        if os.path.exists(output_attributes_path):
            existing_attributes = self.read_attributes_file(file_path, output_attributes_path)
        
        # Detect attributes for columns not in existing attributes
        detected_attributes = self._detect_column_attributes(df)
        
        # Get hierarchy columns for sanitization flag placement
        # Create a temporary df_clean to detect hierarchy columns (same logic as in analyze_file)
        df_clean = df.copy()
        for col in df_clean.columns:
            if df_clean[col].dtype == 'object':
                df_clean[col] = df_clean[col].astype(str)
        hierarchy_columns = self._detect_hierarchy_columns(df_clean)
        
        # Merge existing and detected attributes, giving priority to existing values
        merged_attributes = {}
        
        # First, process all columns in the dataframe
        for column in df.columns:
            if column in existing_attributes:
                # Use existing attribute settings
                merged_attributes[column] = existing_attributes[column].copy()
            else:
                # Use newly detected attributes
                merged_attributes[column] = detected_attributes[column].copy()
            
            # Add sanitise flag to hierarchy columns if sanitization was requested
            if sanitise_requested and column in hierarchy_columns:
                merged_attributes[column]['sanitise'] = True
        
        # Write the attributes file using the simplified format
        with open(output_attributes_path, 'w') as f:
            for column, attrs in merged_attributes.items():
                line = f"{column}: {attrs['type']}"
                
                # Add default value as a missing parameter
                if attrs['default'] is not None:
                    line += f" (missing = {attrs['default']})"
                else:
                    # For boolean values without explicit default, use False
                    if attrs['type'] == 'boolean':
                        line += " (missing = False)"
                    else:
                        line += " (missing = None)"
                
                # Add flags
                flags = []
                if attrs['mandatory']:
                    flags.append("mandatory")
                if attrs.get('sanitise', False):
                    flags.append("sanitise")
                
                if flags:
                    line += "; " + "; ".join(flags)
            
                f.write(line + "\n")
        
        return output_attributes_path


    def _should_add_sanitise_flag(self, column: str, hierarchy_columns: List[str], 
                                sanitise_requested: bool) -> bool:
        """
        Determine if a column should have the sanitise flag in the attributes file.
        """
        # Only add sanitise flag to hierarchy columns when --sanitise-fields was used
        return sanitise_requested and column in hierarchy_columns


    def _detect_column_attributes(self, df: pd.DataFrame) -> Dict[str, Dict[str, Any]]:
        """
        Detect attributes for each column in the dataframe.
        Returns a dictionary of column attributes including type, mandatory status, etc.
        """
        attributes = {}
        
        # Get the first column name - this will be considered mandatory by default
        first_column = df.columns[0] if len(df.columns) > 0 else None
        
        for column in df.columns:
            # Initialize column attributes
            col_attrs = {
                "type": "unknown",
                "mandatory": False,
                "default": None,
                "values": []
            }
            
            # Only the first column is mandatory by default if no attributes file exists
            if column == first_column:
                col_attrs["mandatory"] = True
            
            # Get values for analysis
            values = df[column].dropna().tolist()
            unique_values = df[column].dropna().unique().tolist()
            
            # Calculate some statistics to help determine if this is a categorical field
            total_values = len(values)
            unique_count = len(unique_values)
            
            # If there are no values, we can't determine the type
            if not values:
                col_attrs["type"] = "unknown"
            
            # Check for boolean values
            elif all(self._is_boolean(v) for v in values):
                col_attrs["type"] = "boolean"
                col_attrs["default"] = "False"
            
            # Check for integer values
            elif all(self._is_integer(v) for v in values):
                col_attrs["type"] = "integer"
                col_attrs["default"] = "0"
            
            # Check for paths
            elif column.endswith('_file') or any('/' in str(v) or '\\' in str(v) for v in values):
                col_attrs["type"] = "path"
            
            # Enhanced categorical detection
            elif any([
                # Few unique values compared to total values
                unique_count <= min(10, total_values * 0.25),
                
                # More repetition than expected in random data
                unique_count < total_values * 0.5 and unique_count < 20,
                
                # Special column names that are typically categorical
                column.lower() in ['type', 'category', 'class', 'group', 'status', 'state', 'level'],
                column.lower().endswith('_type') or column.lower().endswith('_class') or 
                column.lower().endswith('_category') or column.lower().endswith('_level'),
                
                # Check for repeating values in string columns
                all(isinstance(v, str) for v in values) and len(set(values)) < min(15, len(values) * 0.5),
                
                # Any string column with at least one repeated value and few unique values
                all(isinstance(v, str) for v in values) and unique_count <= 5 and unique_count < total_values,
                
                # For columns with fewer values, be more aggressive in category detection
                all(isinstance(v, str) for v in values) and total_values <= 10 and unique_count < total_values
            ]):
                col_attrs["type"] = "category"
                col_attrs["values"] = [str(v) for v in unique_values]
            
            # Default to string type
            else:
                col_attrs["type"] = "string"
            
            # Set special defaults for certain columns based on naming conventions
            if column.lower() in ['trim', 'read_qc'] or column.lower().endswith('_qc') or column.lower().endswith('_flag'):
                col_attrs["type"] = "boolean"
                col_attrs["default"] = "False"
            elif column.lower() in ['kmer_size', 'genome_size'] or column.lower().endswith('_size') or column.lower().endswith('_length'):
                col_attrs["type"] = "integer"
                col_attrs["default"] = "0"
            
            attributes[column] = col_attrs
        
        return attributes


    def print_analysis_summary(self) -> None:
        """Print a summary of the analysis results with friendly dog emoticons"""
        # Define version and GitHub link
        version = "v.2025.10.01"
        github_link = "https://github.com/diegomics/GAME/tree/main/scripts/NomNom"
        
        if "file_info" not in self.analysis_results:
            print("(U ಥ ﻌ ಥ)")
            print("Error: No file information available in analysis results")
            return
                
        file_info = self.analysis_results["file_info"]
        file_name = os.path.basename(file_info.get('path', 'Unknown'))
        file_path = file_info.get('path', 'Unknown')
        
        print("°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°")
        print("(∪• ﻌ •∪)")
        print(f"Welcome to NomNom {version}")
        print(f"GitHub: {github_link}")
        
        print("")
        print("°°°°°°°°°°")
        print(f"(∪≖ﻌ≖∪)")
        print(f"I see {file_name} has {file_info.get('rows', 0)} rows and {file_info.get('columns', 0)} columns")
        
        # Handle attributes file messaging with --overwrite-attributes consideration
        attributes_found = file_info.get("attributes_found_initially", False)
        overwrite_requested = file_info.get("overwrite_attributes", False)
        custom_attrs_path = file_info.get("attributes_path", None)
        default_attrs_path = os.path.splitext(file_path)[0] + ".attributes"
        
        if overwrite_requested:
            if attributes_found:
                print("")
                print("°°°°°°°°°°")
                print("(∪⚡ﻌ⚡∪)")
                print("--overwrite-attributes flag detected, regenerating attributes file with fresh detection")
            else:
                print("")
                print("°°°°°°°°°°")
                print("(∪oﻌo∪)")
                print("--overwrite-attributes flag detected, creating new attributes file")
        elif custom_attrs_path and attributes_found:
            print("")
            print("°°°°°°°°°°")
            print("(∪＾ﻌ＾∪)")
            print(f"Using custom attributes file: {custom_attrs_path}")
        elif attributes_found:
            print("")
            print("°°°°°°°°°°")
            print("(∪＾ﻌ＾∪)")
            print(f"Using attributes file found at: {default_attrs_path}")
        else:
            print("")
            print("°°°°°°°°°°")
            print("(∪oﻌo∪)")
            print("No attributes file found, will create one")
        
        print("")
        print("°°°°°°°°°°")
        print("(U ᵔﻌᵔ)っ")
        print(f"Ok, will try to nom nom {file_name}...")
        
        # Print errors if any
        if self.errors:
            print("")
            print("°°°°°°°°°°")
            print("(∪☉ﻌ☉∪)")
            print("ERRORS:")
            
            for error in self.errors:
                print(f"  ❌ {error}")

            print("")    
            print("°°°°°°°°°°")
            print("(U ಥ ﻌ ಥ)")
            print(f"{file_name} couldn't be nommed, exiting!")
            return
                    
        # Print warnings if any
        if self.warnings:
            print("")
            print("°°°°°°°°°°")
            print("(U º ﻌ º)")
            print("WARNINGS:")
            
            for warning in self.warnings:
                print(f"  ⚠️  {warning}")
            print("")
        
        # Success message
        print("")
        print("°°°°°°°°°°")
        print("(∪◕ﻌ◕∪)")
        print(f"Table was nommed! YAML file created.")
        
        # YAML file path info
        yaml_path = f"{os.path.splitext(file_path)[0]}.yaml"
        print("")
        print("°°°°°°°°°°")
        print("(∪¯ﻌ¯∪)")
        print(f"YAML results saved to: {yaml_path}")
        
        # Additional message for overwrite-attributes
        if overwrite_requested:
            attrs_path = os.path.splitext(file_path)[0] + ".attributes"
            print("")
            print("°°°°°°°°°°")
            print("(∪🔄ﻌ🔄∪)")
            print(f"Attributes file regenerated: {attrs_path}")
        
        print("")


class YAMLToTable:
    """
    A tool for converting YAML files generated by TableAnalyzer back to tabular format.
    Features include:
    - Automatic flattening of hierarchical YAML structure
    - Support for multiple output formats (CSV, TSV, Excel)
    - Field type preservation using attributes files
    """
    
    def __init__(self):
        self.warnings = []
        self.errors = []
        self.df = None  # Will store the generated DataFrame
        
    def read_yaml(self, yaml_path: str) -> Dict[str, Any]:
        """
        Read a YAML file and return its contents as a dictionary
        """
        try:
            with open(yaml_path, 'r') as f:
                yaml_data = yaml.safe_load(f)
                return yaml_data
        except Exception as e:
            self.errors.append(f"Failed to read YAML file: {str(e)}")
            return {}
    
    def yaml_to_dataframe(self, yaml_data: Dict[str, Any]) -> pd.DataFrame:
        """
        Convert a hierarchical YAML structure to a flat pandas DataFrame
        """
        # Check if this is a NomNom-generated YAML with Sample_ID root
        if 'Sample_ID' not in yaml_data:
            self.errors.append("Invalid YAML format: 'Sample_ID' root key not found")
            return pd.DataFrame()
        
        sample_data = yaml_data['Sample_ID']
        
        # Prepare data for DataFrame
        rows = []
        
        # Process each sample ID
        for sample_id, sample_info in sample_data.items():
            # Process each data type/subsection within the sample
            for data_type, type_data in sample_info.items():
                # Create a row with sample ID and data type
                row = {'Sample_ID': sample_id, 'Data_Type': data_type}
                
                # Add all other fields
                self._flatten_dict(type_data, row)
                
                # Add the row to our collection
                rows.append(row)
        
        # Convert rows to DataFrame
        if not rows:
            self.warnings.append("No data found in YAML file")
            return pd.DataFrame()
            
        df = pd.DataFrame(rows)
        
        # Ensure Sample_ID is first column, followed by Data_Type
        cols = df.columns.tolist()
        if 'Sample_ID' in cols and 'Data_Type' in cols:
            cols.remove('Sample_ID')
            cols.remove('Data_Type')
            cols = ['Sample_ID', 'Data_Type'] + cols
            df = df[cols]
            
        return df
    
    def _flatten_dict(self, nested_dict: Dict[str, Any], result: Dict[str, Any], prefix: str = '') -> None:
        """
        Recursively flatten a nested dictionary into key-value pairs
        
        Args:
            nested_dict: Nested dictionary to flatten
            result: Dictionary to store the flattened key-value pairs
            prefix: Prefix for keys in the current level of nesting
        """
        for key, value in nested_dict.items():
            new_key = f"{prefix}{key}" if prefix else key
            
            # Handle nested dictionaries with special case for read paths
            if isinstance(value, dict):
                # Special case for read paths dictionary
                if all(k.startswith('Read') for k in value.keys()):
                    # Convert read paths dict to comma-separated string
                    read_paths = []
                    for read_num in sorted(value.keys(), key=lambda x: int(x.replace('Read', ''))):
                        read_paths.append(value[read_num])
                    result[new_key] = ', '.join(read_paths)
                else:
                    # Regular nested dict - recurse with updated prefix
                    self._flatten_dict(value, result, f"{new_key}.")
            else:
                # Convert boolean values to strings that match NomNom's expected format
                if isinstance(value, bool):
                    result[new_key] = str(value).lower()
                # Ensure None values are represented as "None" string
                elif value is None:
                    result[new_key] = "None"
                else:
                    result[new_key] = value
    
    def convert_to_table(self, yaml_path: str, output_path: str = None, 
                         attributes_path: str = None, output_format: str = 'csv',
                         delimiter: str = ',') -> str:
        """
        Convert a YAML file to a tabular format
        
        Args:
            yaml_path: Path to the YAML file to convert
            output_path: Path for the output file (if None, derives from yaml_path)
            attributes_path: Optional path to attributes file for type information
            output_format: Format for the output file ('csv', 'tsv', 'excel')
            delimiter: Delimiter to use for CSV/TSV output
            
        Returns:
            Path to the created file
        """
        # Read YAML data
        yaml_data = self.read_yaml(yaml_path)
        if not yaml_data or self.errors:
            return None
        
        # Convert to DataFrame
        self.df = self.yaml_to_dataframe(yaml_data)
        if self.df.empty:
            return None
        
        # Determine output path if not provided
        if not output_path:
            base_path = os.path.splitext(yaml_path)[0]
            if output_format.lower() == 'excel':
                output_path = f"{base_path}.xlsx"
            elif output_format.lower() == 'tsv':
                output_path = f"{base_path}.tsv"
            else:  # Default to CSV
                output_path = f"{base_path}.csv"
        
        # Apply type conversions based on attributes file if provided
        if attributes_path:
            self._apply_attributes(attributes_path)
        
        # Save the table in requested format
        try:
            if output_format.lower() == 'excel':
                self.df.to_excel(output_path, index=False)
            elif output_format.lower() == 'tsv':
                self.df.to_csv(output_path, sep='\t', index=False)
            else:  # Default to CSV
                self.df.to_csv(output_path, sep=delimiter, index=False)
                
            return output_path
        except Exception as e:
            self.errors.append(f"Failed to write output file: {str(e)}")
            return None
    
    def _apply_attributes(self, attributes_path: str) -> None:
        """
        Apply type conversions to DataFrame based on attributes file
        
        Args:
            attributes_path: Path to the attributes file
        """
        if not os.path.exists(attributes_path):
            self.warnings.append(f"Attributes file not found: {attributes_path}")
            return
            
        try:
            attributes = {}
            
            # Read the attributes file
            with open(attributes_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    # Split the line into column name and attributes
                    parts = line.split(':', 1)
                    if len(parts) != 2:
                        continue
                    
                    column = parts[0].strip()
                    attrs_str = parts[1].strip()
                    
                    # Parse attributes - we're mainly interested in type
                    attrs_parts = [p.strip() for p in attrs_str.split(',')]
                    type_part = attrs_parts[0]
                    
                    # Strip out any category values
                    if '(' in type_part:
                        col_type = type_part.split('(', 1)[0].strip()
                    else:
                        col_type = type_part.strip()
                    
                    attributes[column] = col_type
                    
            # Apply type conversions for each column
            for column, col_type in attributes.items():
                if column in self.df.columns:
                    if col_type == 'boolean':
                        # Convert string booleans to actual booleans
                        self.df[column] = self.df[column].apply(
                            lambda x: str(x).lower() in ['true', 'yes', 'y', 't', '1'])
                    elif col_type == 'integer':
                        # Convert to integers, handling errors
                        try:
                            self.df[column] = pd.to_numeric(self.df[column], errors='coerce').fillna(0).astype(int)
                        except Exception:
                            self.warnings.append(f"Could not convert column {column} to integer")
                    
        except Exception as e:
            self.warnings.append(f"Error processing attributes file: {str(e)}")


def is_yaml_file(file_path: str) -> bool:
    """
    Determine if a file is a YAML file based on extension and content
    """
    # Check file extension first
    file_ext = os.path.splitext(file_path)[1].lower()
    
    # Definitive YAML extensions
    if file_ext in ['.yaml', '.yml']:
        return True
    
    # Definitive non-YAML extensions (all supported formats)
    non_yaml_extensions = [
        # Microsoft Office
        '.xlsx', '.xls', '.xlsm', '.xltx', '.xltm',
        # LibreOffice/OpenOffice
        '.ods', '.ots', '.fods',
        # Apple
        '.numbers',
        # Text formats
        '.csv', '.tsv', '.txt', '.dat',
        # Other data formats
        '.json', '.xml', '.parquet', '.feather', '.h5', '.hdf5'
    ]
    
    if file_ext in non_yaml_extensions:
        return False
    
    # If extension isn't definitive, check the content
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            # Read first few lines
            lines = []
            for _ in range(5):
                line = f.readline()
                if not line:
                    break
                lines.append(line.strip())
            
            if not lines:
                return False
            
            # Check for YAML indicators
            yaml_indicators = []
            
            # YAML documents often start with ---
            yaml_indicators.append(any(line.startswith('---') for line in lines))
            
            # YAML has key: value pairs (but so do other formats)
            colon_lines = [line for line in lines if ':' in line and not line.strip().startswith('#')]
            yaml_indicators.append(len(colon_lines) >= 2)
            
            # Specific NomNom YAML indicator
            yaml_indicators.append(any('Sample_ID:' in line for line in lines))
            
            # Strong YAML structure indicators
            yaml_indicators.append(any(line.strip().endswith(':') for line in lines))  # YAML sections
            yaml_indicators.append(any(line.startswith('  ') or line.startswith('\t') for line in lines))  # Indentation
            
            # If multiple YAML indicators are found, it's likely a YAML file
            return sum(yaml_indicators) >= 2
            
    except UnicodeDecodeError:
        # If we can't decode as text, it's definitely not YAML
        return False
    except Exception:
        return False
        
    return False


def main():
    """Command-line interface for NomNom with bi-directional conversion"""
    parser = argparse.ArgumentParser(
        description="NomNom: A tool for dealing with genomic data tables"
    )
    
    parser.add_argument("file", help="Path to the input file (table or YAML)")
    
    # Output options
    parser.add_argument(
        "--output", "-o", 
        help="Path to save output (default: <input_file>.yaml or <input_file>.csv)",
        default=None
    )
    
    # Format options
    parser.add_argument(
        "--format", "-f",
        help="Output format when converting YAML to table (csv, tsv, excel)",
        choices=['csv', 'tsv', 'excel'],
        default='csv'
    )
    
    # Attributes file option
    parser.add_argument(
        "--attributes", "-a",
        help="Path to custom attributes file (default: <input_file>.attributes)",
        default=None
    )
    
    # Processing options
    parser.add_argument(
        "--sanitise-fields", "-s",
        help="Replace spaces and special characters with underscores in hierarchical field values for YAML compatibility",
        action="store_true"
    )
    parser.add_argument(
        "--overwrite-attributes",
        help="Regenerate .attributes file using current smart detection, ignoring existing file",
        action="store_true"
    )
    parser.add_argument(
        "--validate-paths", "-v",
        help="Validate that all file paths in the YAML exist (fails if any path is missing)",
        action="store_true"
    )
    parser.add_argument(
        "--discover-paths", "-D",
        help="Expand glob patterns (*, ?, []) in file path columns to actual file paths",
        action="store_true"
    )
    parser.add_argument(
        "--overwrite-yaml",
        help="Allow overwriting existing YAML output file if it exists",
        action="store_true"
    )
    
    # Other options
    parser.add_argument(
        "--quiet", "-q",
        help="Suppress console output, only save files",
        action="store_true"
    )
    parser.add_argument(
        "--force-delimiter", "-d",
        help="Force a specific delimiter for table files (e.g., ';', ',', '|', 'tab')",
        default=None
    )
    parser.add_argument(
        "--debug", 
        help="Show detailed debugging information",
        action="store_true"
    )
    
    args = parser.parse_args()
    
    # Determine if the input is a YAML or table file
    input_is_yaml = is_yaml_file(args.file)
    
    if input_is_yaml:
        # Convert YAML to table
        if not args.quiet:
            print("°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°")
            print("(∪• ﻌ •∪)")
            print("NomNom detected YAML input, converting to table")
        
        converter = YAMLToTable()
        output_path = converter.convert_to_table(
            args.file, 
            args.output, 
            args.attributes, 
            args.format,
            ',' if args.force_delimiter is None else args.force_delimiter
        )
        
        if converter.errors or output_path is None:
            if not args.quiet:
                print("\n°°°°°°°°°°")
                print("(U ಥ ﻌ ಥ)")
                print("ERRORS:")
                for error in converter.errors:
                    print(f"  ❌ {error}")
                print("Conversion failed!")
            return 1
            
        if not args.quiet:
            print("\n°°°°°°°°°°")
            print("(∪◕ﻌ◕∪)")
            print(f"YAML successfully converted to {args.format.upper()}!")
            print(f"Output saved to: {output_path}")
            
            if converter.warnings:
                print("\n°°°°°°°°°°")
                print("(U º ﻌ º)")
                print("WARNINGS:")
                for warning in converter.warnings:
                    print(f"  ⚠️  {warning}")
        
        return 0
        
    else:
        # Process table file
        analyzer = TableAnalyzer()
        
        try:
            # Process with specified options - pass the new flag
            results = analyzer.analyze_file(args.file, args.attributes, 
                                          args.sanitise_fields, args.overwrite_attributes,
                                          args.validate_paths, args.discover_paths)
            
            if results is None or analyzer.errors:
                if not args.quiet:
                    analyzer.print_analysis_summary()
                return 1
            
            # Determine output YAML path
            if args.output:
                yaml_path = args.output
            else:
                file_base = os.path.splitext(args.file)[0]
                yaml_path = f"{file_base}.yaml"
            
            # Check if output YAML already exists
            if os.path.exists(yaml_path):
                if not args.overwrite_yaml:
                    # File exists and overwrite not allowed - throw error
                    if not args.quiet:
                        print("\n°°°°°°°°°°")
                        print("(U ಥ ﻌ ಥ)")
                        print("ERROR:")
                        print(f"  ❌ Output YAML file already exists: {yaml_path}")
                        print(f"  ❌ To prevent accidental overwrites, NomNom will not continue.")
                        print(f"  ❌ Solutions:")
                        print(f"      1. Remove the existing file: rm {yaml_path}")
                        print(f"      2. Use a different output name: -o another_name.yaml")
                        print(f"      3. Allow overwriting: --overwrite-yaml")
                    return 1
                else:
                    # Overwrite allowed - add warning
                    if not args.quiet:
                        print(f"⚠️  Warning: Overwriting existing file: {yaml_path}")
            
            # Save YAML
            analyzer.save_yaml(yaml_path)
            
            # Print analysis summary unless quiet mode is enabled
            if not args.quiet:
                analyzer.print_analysis_summary()
                
            return 0
            
        except Exception as e:
            if not args.quiet:
                print(f"Error analyzing file: {str(e)}")
                if args.debug:
                    import traceback
                    traceback.print_exc()
            return 1


# Add debug wrappers to relevant methods
TableAnalyzer._is_boolean = debug_wrapper(TableAnalyzer._is_boolean)
TableAnalyzer._to_boolean = debug_wrapper(TableAnalyzer._to_boolean)
TableAnalyzer._apply_defaults = debug_wrapper(TableAnalyzer._apply_defaults)
TableAnalyzer._build_structured_data = debug_wrapper(TableAnalyzer._build_structured_data)
TableAnalyzer._parse_read_paths = debug_wrapper(TableAnalyzer._parse_read_paths)
TableAnalyzer.read_file = debug_wrapper(TableAnalyzer.read_file)
TableAnalyzer._detect_delimiter = debug_wrapper(TableAnalyzer._detect_delimiter)


if __name__ == "__main__":
    sys.exit(main())
