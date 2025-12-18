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

### 1) Get the Latest Release

```bash
git clone --depth=1 --branch=main https://github.com/diegomics/GEP2.git
```

### 2) Create the GEP2 Conda Environment

The environment contains Snakemake and other required packages:

```bash
cd GEP2
conda env create -f install.yml
```

### 3) Enter Your Data in a Table

You can use Google Drive, Excel, LibreOffice, Numbers, CSV, TSV, etc.

The table should contain these columns:

| sp_name | asm_id | asm_files | read_type | read_files |
|---------|--------|-----------|-----------|------------|
|         |        |           |           |            |

**Please see the [example table](https://docs.google.com/spreadsheets/d/1xmsstJGBo45SEQgCPncE76u51_VN5IDG9sFFafYkGfI/edit?gid=1029606022#gid=1029606022)**

#### Column Descriptions:

- **sp_name**: Species name in binomial nomenclature (e.g., `Vultur gryphus`)
- **asm_id**: Assembly identifier (e.g., `hifiasm_l2`, `yahs_test`, or `ASM2260516v1`)
- **asm_files**: Path to assembly file, URL, or accession number (e.g., `GCA_022605165.1`). If it's a link or accession, the pipeline will download the data automatically. If Pri/Alt or (Hap1/Hap2) assemblies available, add as comma-separated, like: `GCA_963854735.1, GCA_963694935.1`
- **read_type**: Can be `illumina`, `10x`, `hifi`, or `ont` (variations like `PacBio`, `paired-end`, `linked-read`, `arima`, `promethion` are also accepted)
- **read_files**: Comma-separated list of paths to read files. Can also be accession numbers (e.g., `ERR12205285,ERR12205286`). For paired-end reads, list as: `forward1,reverse1,forward2,reverse2`

### 4) Configure the Control Panel

Add the table path/address and select different options in:

```
config/control_panel.yaml
```

### 5) Configure Cluster or Computer Parameters

**Note:** Local mode has not been fully tested yet.

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
        ├── k{k-mer_length}/
        │   ├── {asm_id}.hist
        │   ├── {asm_id}.meryl
        │   └── genomescope2/
        │       └── {asm_id}_linear_plot.png
        ├── logs/
        └── merqury/
            └── ...
```

