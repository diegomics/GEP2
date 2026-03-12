#!/usr/bin/env python3

# Little resource calculator AKA "Loris appeacer attempt"
# by Diego De Panis, 2025
# This script is part of the GAME-pipeline
# See https://github.com/diegomics/GAME/tree/main/scripts/misc


# This script takes:
# 1) user-defined limits in config/control_panel.yaml (myCPU, myRAM, myRUN)
# 2) tools settings in config/resources.yaml like:
#    resources:
#      <tool_name_1>:
#        <cpu_key>: <value>
#        <ram_key>: <value>
#        <run_key>: <value>                  
#      <tool_name_1>:
#        <cpu_key>: <value>
#        <ram_key>: <value>
#        <run_key>: <value>   

# Please see the README.md for more details


import sys
import yaml
from pathlib import Path
from typing import List, Union

class ResourceCalculator:
    def __init__(self, resource_config_path: str, user_config: dict):
        with open(resource_config_path) as f:
            self.resource_config = yaml.safe_load(f)
        
        self.user_config = user_config
        self.compression_ratios = self.resource_config.get('compression_ratios', {})
        
        # User limits
        self.user_cpu = user_config.get('myCPU', 1)
        self.user_ram = user_config.get('myRAM', 4)
        self.user_runtime = user_config.get('myRUN', 1)
        
        self._validate_user_constraints()
    
    def _validate_user_constraints(self):
        """Validate user resource limits"""
        if self.user_ram < 4:
            raise ValueError(f"myRAM must be at least 4GB (got {self.user_ram}GB)")
        if self.user_cpu < 1:
            raise ValueError(f"myCPU must be at least 1 (got {self.user_cpu})")
        if self.user_runtime < 1:
            raise ValueError(f"myRUN must be at least 1 hour (got {self.user_runtime})")
    
    def _get_file_size_gb(self, filepath: Union[str, Path]) -> float:
        """Get file size in GB"""
        try:
            if isinstance(filepath, str):
                filepath = Path(filepath)
            
            if filepath.exists():
                return filepath.stat().st_size / (1024**3)
            else:
                # Estimate for non-existent files
                return self._estimate_file_size(filepath)
        except:
            return 1.0  # Fallback
    
    def _estimate_file_size(self, filepath: Path) -> float:
        """Estimate file size based on filename patterns"""
        filename_lower = str(filepath).lower()
        
        # Conservative estimates for common genomic files
        if any(x in filename_lower for x in ['.fastq', '.fq']):
            return 5.0  # 5GB for read files
        elif any(x in filename_lower for x in ['.bam', '.cram']):
            return 10.0  # 10GB for alignment files
        elif any(x in filename_lower for x in ['.fa', '.fasta', '.fna']):
            return 3.0  # 3GB for genome files
        else:
            return 1.0
    
    def _get_compression_ratio(self, filepath: Union[str, Path]) -> float:
        """Get compression ratio to estimate uncompressed size"""
        filename = str(filepath).lower()
        
        for ext, ratio in self.compression_ratios.items():
            if ext != 'default' and ext in filename:
                return ratio
        
        return self.compression_ratios.get('default', 1.0)
    
    def _estimate_uncompressed_size(self, input_files: List) -> float:
        """Estimate total uncompressed size of input files"""
        if not input_files:
            return 0.0
        
        total_uncompressed = 0.0
        for filepath in input_files:
            compressed_size = self._get_file_size_gb(filepath)
            compression_ratio = self._get_compression_ratio(filepath)
            uncompressed_size = compressed_size * compression_ratio
            total_uncompressed += uncompressed_size
        
        return total_uncompressed
    
    def calculate_cpu(self, rule_name: str, input_files: List = None, attempt: int = 1) -> int:
        """Calculate CPU requirements with adaptive bounds"""
        rule_config = self.resource_config['resources'].get(rule_name, {})
        max_cpu = rule_config.get('max_cpu', 1)
        
        # Calculate optimal CPU: use all available cores up to tool's maximum
        cpu = min(self.user_cpu, max_cpu)
        
        # Apply retry scaling
        if attempt > 1:
            cpu = cpu * (1.5 ** (attempt - 1))
        
        # Apply final user limit (safety cap)
        cpu = min(int(cpu), self.user_cpu)
        
        return cpu
    
    def calculate_memory(self, rule_name: str, input_files: List = None, attempt: int = 1) -> int:
        """Calculate memory requirements in GB"""
        rule_config = self.resource_config['resources'].get(rule_name, {})
        
        # Handle simple static memory specification
        if 'mem_gb' in rule_config:
            memory = rule_config['mem_gb']
        else:
            # Calculate based on uncompressed file size
            uncompressed_size = self._estimate_uncompressed_size(input_files or [])
            
            base_mem = rule_config.get('base_mem_gb', 4)
            mem_per_gb = rule_config.get('mem_gb_per_uncompressed_gb', 1.0)
            
            memory = base_mem + (uncompressed_size * mem_per_gb)
        
        # Apply retry scaling
        if attempt > 1:
            memory = memory * (1.5 ** (attempt - 1))
        
        # Check if we need more than available (before capping)
        if memory > self.user_ram:
            print(f"⚠️  Warning: {rule_name} ideally needs {memory:.1f}GB but only {self.user_ram}GB available")
            print(f"   Job may run slower or fail due to insufficient memory")
        
        # Apply user limit (cap the allocation)
        memory = min(memory, self.user_ram)
        
        return int(memory)
    
    def calculate_runtime(self, rule_name: str, input_files: List = None, attempt: int = 1) -> int:
        """Calculate runtime requirements in minutes"""
        rule_config = self.resource_config['resources'].get(rule_name, {})
        
        # Handle simple static runtime specification
        if 'runtime_min' in rule_config:
            runtime = rule_config['runtime_min']
        else:
            # Calculate based on uncompressed file size
            uncompressed_size = self._estimate_uncompressed_size(input_files or [])
            
            base_runtime = rule_config.get('base_runtime_min', 30)
            runtime_per_gb = rule_config.get('runtime_min_per_uncompressed_gb', 10)
            
            runtime = base_runtime + (uncompressed_size * runtime_per_gb)
        
        # Apply retry scaling
        if attempt > 1:
            runtime = runtime * (1.5 ** (attempt - 1))
        
        # Apply user limit (convert hours to minutes)
        max_runtime = self.user_runtime * 60
        runtime = min(runtime, max_runtime)
        
        return int(runtime)

# Convenience functions for Snakemake integration
def create_resource_functions(calculator: ResourceCalculator):
    """Create lambda functions for Snakemake resource specifications"""
    
    def cpu_func(rule_name):
        return lambda wildcards, input=None, attempt=1: calculator.calculate_cpu(
            rule_name, list(input) if input else [], attempt
        )
    
    def mem_func(rule_name):
        return lambda wildcards, input=None, attempt=1: calculator.calculate_memory(
            rule_name, list(input) if input else [], attempt
        ) * 1024  # Convert to MB for Snakemake
    
    def time_func(rule_name):
        return lambda wildcards, input=None, attempt=1: calculator.calculate_runtime(
            rule_name, list(input) if input else [], attempt
        )
    
    return cpu_func, mem_func, time_func
