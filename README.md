# GEP2

**Genome Evaluation Pipeline (v2)**

This repository contains a significantly updated version of [GEP](https://git.imp.fu-berlin.de/begendiv/gep) that builds upon lessons learned from [ERGA](https://www.erga-biodiversity.eu/team-1/sac---sequencing-and-assembly-committee) and [GAME](https://github.com/diegomics/GAME).

Data is entered via a simple table, and configuration is managed through a tidy control panel. GEP2 uses a modern [Snakemake](https://snakemake.readthedocs.io) version with [containers](https://apptainer.org) and can run on a server/cluster (SLURM) or a local computer.

**Please cite:** Genome Evaluation Pipeline (GEP): A fully-automated quality control tool for parallel evaluation of genome assemblies. https://doi.org/10.1093/bioadv/vbaf147

---

## Requirements

- `conda`
- `apptainer`

---

## How to Get and Set Up GEP2

### 1) Get the latest version

GEP2 is adding features rapidly, so please download the latest [release](https://github.com/diegomics/GEP2/releases) (or clone the repo for getting hot fixes faster!)

### 2) Create the GEP2 Conda Environment

The environment contains Snakemake packages and [NomNom](https://github.com/diegomics/NomNom). Enter the GEP folder and:

```bash
conda env create -f install.yml
```

### 3) Enter Your Data in a Table

You can use Google Drive, Excel, LibreOffice, Numbers, CSV, TSV, etc.

The table should contain these columns:

| sp_name | asm_id | skip | asm_files | read_type | read_files |
|---------|--------|------|-----------|-----------|------------|
|         |        |      |           |           |            |

**Please see the [example table](https://docs.google.com/spreadsheets/d/1xmsstJGBo45SEQgCPncE76u51_VN5IDG9sFFafYkGfI/edit?gid=1029606022#gid=1029606022)**.Please see the example table. **The easiest is to make a copy of that Google table** (`File`->`Make a copy`) and replace the fields with your data. Remember to change permissions (`Share`-> change `General access` to "Anyone with the link" `viewer`)

#### Column Descriptions:

- **sp_name**: Species name in binomial nomenclature (e.g., `Vultur gryphus`)
- **asm_id**: Assembly identifier (e.g., `hifiasm_l2`, `yahs_test`, or `ASM2260516v1`)
- **skip**: Flag assemblies for selective analysis skipping. Leave empty (or `-` or `off`) to run all analyses. Set to `on` to flag this assembly, then control which analyses to skip in `control_panel.yaml` using `SKIP_KMER`, `SKIP_INSP`, `SKIP_HIC`, etc. Useful for running quick QC on draft assemblies while running full analysis on final assemblies.
- **asm_files**: Path to assembly file, URL, or accession number (e.g., `GCA_022605165.1`). If it's a link or accession, the pipeline will download the data automatically. If Pri/Alt or (Hap1/Hap2) assemblies available, add as comma-separated, like: `GCA_963854735.1, GCA_963694935.1`
- **read_type**: Can be `illumina`, `10x`, `hifi`, or `ont` (variations like `PacBio`, `paired-end`, `linked-read`, `arima`, `promethion` and others should also work fine)
- **read_files**: Comma-separated list of paths to read files. Can also be accession numbers (e.g., `ERR12205285,ERR12205286`). For paired-end reads, list as: `forward1,reverse1,forward2,reverse2`. Also can use pattern expansion in paths, like `/readsA/*.fq.gz, /readsB/*.fq.gz` 


### 4) Configure the Control Panel

Add the table path/address and select different options in:

```
config/control_panel.yaml
```

### 5) Configure Cluster or Computer Parameters

```
GEP2/execution/
├── local/
│   └── config.yaml
└── slurm/
    └── config.yaml
```

**IMPORTANT:** You can tweak per-tool resources boundaries in `GEP2/config/resources.yaml`

### 6) Run!

**First run takes longer** as containers need to be built.

Load the conda environment like `conda activate GEP2_env` and in the GEP2 folder run:

#### On HPC/Server/Cluster using [Slurm](https://slurm.schedmd.com/):
```bash
nohup snakemake --profile execution/slurm &
```

#### On Local Computer:
```bash
nohup snakemake --profile execution/local &
```

#### About the Command:
- `nohup` runs Snakemake in a way that won't be interrupted if you lose connection to the server/cluster
- The trailing `&` runs the command in the background, allowing you to continue using the terminal

#### Dry Run (Recommended):
Before running the full pipeline, perform a dry run to check what will execute and catch any errors:

```bash
snakemake --profile execution/slurm --dry-run
```

You can also inspect:
- `GEP2_results/data_config.yaml`
- `GEP2_results/download_manifest.json`

---

## Results Structure

Open the report with a markdown renderer (VS Code works well):

```
GEP2_results/{sp_name}/{asm_id}/{asm_id}_report.md
```

### Directory Structure:

```
GEP2_results/
├── data/
│   └── {sp_name}/
│       └── reads/
│           ├── {read_type}/
│           │   ├── {read_symlink}
│           │   ├── kmer_db_k{k-mer_length}/
│           │   │   └── {read_name}.meryl
│           │   ├── logs/
│           │   └── processed/
│           │       ├── {read_type}_Path{number}_{read_name}_{process}.fq.gz
│           │       └── reports/
│           │           └── multiqc_report.html
│           └── ...
├── data_config.yaml
├── data_table_{hash}.csv
├── downloaded_data/
│   └── {sp_name}/
│       ├── assemblies/
│       │   └── {asm_file}
│       └── reads/
│           └── {read_type}/
│               └── {read_file}
├── download_manifest.json
└── {sp_name}/
    └── {asm_id}/
        ├── {asm_id}_report.md
        ├── compleasm/
        │   └── {asm_file_name}/
        │       ├── {asm_file_name}_results.tar.gz
        │       └── {asm_file_name}_summary.txt
        ├── gfastats/
        │   └── {asm_file_name}_stats.txt
        ├── hic/
        │   └── {asm_file_name}/
        │       ├── {asm_file_name}.cool
        │       ├── {asm_file_name}.mcool
        │       ├── {asm_file_name}.pairs.gz
        │       ├── {asm_file_name}.pairtools_stats.txt
        │       ├── {asm_file_name}.pretext
        │       ├── {asm_file_name}_tracks.pretext
        │       ├── {asm_file_name}_snapshots
        │       │   └── {asm_file_name}_FullMap.png
        │       └── tracks
        │           └── ...bedgraph
        ├── inspector/
        │   └── {asm_file_name}/
        │       ├── ..
        │       └── summary_statistics
        ├── k{k-mer_length}/
        │   ├── {asm_id}.hist
        │   ├── {asm_id}.meryl
        │   └── genomescope2/
        │       └── {asm_id}_linear_plot.png
        ├── logs/
        └── merqury/
            ├── ..
            ├── {asm_file_name}.completeness.stats
            ├── {asm_file_name}.qv
            └── ...png
```

### Main tools:

| tool | doi |
| :--- | :--- |
|[blobtools](https://github.com/genomehubs/blobtoolkit) | - |
|[chromap](https://github.com/haowenz/chromap) | 10.1038/s41467-021-26865-w |
|[cooler](https://github.com/open2c/cooler) | 10.1093/bioinformatics/btz540 |
|[compleasm](https://github.com/huangnengCSU/compleasm) | 10.1093/bioinformatics/btad595 |
|[diamond](https://github.com/bbuchfink/diamond) | 10.1038/s41592-021-01101-x |
|[enabrowsertools](https://github.com/enasequence/enaBrowserTools) | - |
|[fastp](https://github.com/OpenGene/fastp) | 10.1093/bioinformatics/bty560 |
|[fastqc](https://github.com/s-andrews/FastQC) | - |
|[fcs-gx](https://github.com/ncbi/fcs-gx) | 10.1186/s13059-024-03198-7 |
|[genomescope2](https://github.com/tbenavi1/genomescope2.0) | 10.1038/s41467-020-14998-3 |
|[gfastats](https://github.com/vgl-hub/gfastats) | 10.1093/bioinformatics/btac460 |
|[hifiadapterfilt](https://github.com/sheinasim-USDA/HiFiAdapterFilt) | 10.1186/s12864-022-08375-1 |
|[inspector](https://github.com/Maggi-Chen/Inspector) | 10.1186/s13059-021-02527-4 |
|[longdust](https://github.com/lh3/longdust) | - |
|[merqury](https://github.com/marbl/merqury) | 10.1186/s13059-020-02134-9 |
|[minimap](https://github.com/lh3/minimap2) | 10.1093/bioinformatics/bty191 |
|[multiqc](https://github.com/MultiQC/MultiQC) | 10.1093/bioinformatics/btw354 |
|[nanoplot](https://github.com/wdecoster/NanoPlot) | 10.1093/bioinformatics/btad311 |
|[pairtools](https://github.com/open2c/pairtools) | 10.1101/2023.02.13.528389 |
|[pretextmap](https://github.com/sanger-tol/PretextMap) | - |
|[sambamba](https://github.com/biod/sambamba) | 10.1093/bioinformatics/btv098 |
|[samtools](https://github.com/samtools/samtools) | 10.1093/gigascience/giab008 |
|[sdust](https://github.com/lh3/sdust) | - |
|[tidk](https://github.com/tolkit/telomeric-identifier) | 10.1093/bioinformatics/btaf049 |

