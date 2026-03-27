# GEP2

**Genome Evaluation Pipeline (v2)**

This repository contains a significantly updated version of [GEP](https://git.imp.fu-berlin.de/begendiv/gep) that builds upon lessons learned from [ERGA](https://www.erga-biodiversity.eu/team-1/sac---sequencing-and-assembly-committee) and [GAME](https://github.com/diegomics/GAME).

Data is entered via a simple table, and configuration is managed through a tidy control panel. GEP2 uses a modern [Snakemake](https://snakemake.readthedocs.io) version with [containers](https://apptainer.org) and can run on a server/cluster (SLURM) or a local computer.

**Please cite:** Genome Evaluation Pipeline (GEP): A fully-automated quality control tool for parallel evaluation of genome assemblies. https://doi.org/10.1093/bioadv/vbaf147

---

## Requirements

- `Conda` (for Snakemake and NomNom)
- `Apptainer`

---

## GEP2 can:
```
• download assemblies & reads (or use the ones in your local storage)
• trim/filter/qc reads (paired-end, 10x, HiFi, ONT)
• classic contiguity metrics and other assembly stats (e.g., N50...)
• completeness metrics based on single-copy orthologs
• kmer-based analyses
• long-reads-based analysis
• Hi-C analysis and curation files production
• contamination screening
• aggregate results and show EBP-metrics guidance
```

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

**Please see the [example table](https://docs.google.com/spreadsheets/d/1xmsstJGBo45SEQgCPncE76u51_VN5IDG9sFFafYkGfI/edit?gid=1029606022#gid=1029606022)**. **The easiest is to make a copy of that Google table** (`File`->`Make a copy`) and replace the fields with your data. Remember to change permissions (`Share`-> change `General access` to "Anyone with the link" `viewer`)

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

Depending on which mode you will run, configure the respective parameters in the config file. **Important:** Don't forget to bind the folders in the _apptainer-args_ field.
```
GEP2/execution/
├── local/
│   └── config.yaml
└── slurm/
    └── config.yaml
```

**Note:** You can tweak per-tool resources boundaries in `GEP2/config/resources.yaml`

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

#### Common Troubleshooting:
- If your process was killed or stopped abruptly, Snakemake might complain about incomplete files when you try to run it again. We can tell Snakemake to identify and rerun those incomplete parts (and always a good idea to try a dry-run first):
```bash
snakemake --profile execution/local --rerun-incomplete —dry-run
```
- If Snakemake was suddenly killed, it might leave a hidden lock on your working directory to prevent other processes from overwriting files, and will tell you the directory is locked. You need to unlock it first before proceed with the run command: `snakemake --unlock`

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
        ├── decontamination/
        │   ├── fcs-gx/
        │   │   └── {asm_file_name}/
        │   │       ├── {asm_file_name}.fcs_gx_report.txt
        │   │       └── {asm_file_name}.taxonomy.rpt
        │   └── blobtools/
        │       └── {asm_file_name}/
        │           ├── Blobdir/
        │           ├── ..
        │           ├── ..blob.circle.png
        │           ├── ..cumulative_plot.png
        │           └── ..snail_plot.png
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

| tool | doi | version | container |
| :--- | :--- | :--- | :--- |
|[bedtools](https://github.com/arq5x/bedtools2) | 10.1093/bioinformatics/btq033 | 2.31.1 | docker://diegomics/hic_analysis:0.2 |
|[blobtools](https://github.com/genomehubs/blobtoolkit) | - | 4.5.1 | docker://genomehubs/blobtoolkit:4.5.1 |
|[busco](https://gitlab.com/ezlab/busco) | 10.1093/molbev/msab199 | 6.0.0 | docker://diegomics/gep2_base:0.2 |
|[bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) | 10.1109/IPDPS.2019.00041 | 2.3 | docker://diegomics/hic_analysis:0.2 |
|[chromap](https://github.com/haowenz/chromap) | 10.1038/s41467-021-26865-w | 0.3.2 | docker://diegomics/hic_analysis:0.2 |
|[cooler](https://github.com/open2c/cooler) | 10.1093/bioinformatics/btz540 | 0.10.4 | docker://diegomics/hic_analysis:0.2 |
|[compleasm](https://github.com/huangnengCSU/compleasm) | 10.1093/bioinformatics/btad595 | 0.2.7 | docker://quay.io/biocontainers/compleasm:0.2.7--pyh7e72e81_0 |
|[diamond](https://github.com/bbuchfink/diamond) | 10.1038/s41592-021-01101-x | 2.1.24 | docker://diegomics/hic_analysis:0.2 |
|[enabrowsertools](https://github.com/enasequence/enaBrowserTools) | - | 1.7.2 | docker://diegomics/gep2_base:0.2 |
|[fastp](https://github.com/OpenGene/fastp) | 10.1093/bioinformatics/bty560 | 1.1.0 | docker://diegomics/gep2_base:0.2 |
|[fastqc](https://github.com/s-andrews/FastQC) | - | 0.12.1 | docker://diegomics/gep2_base:0.2 |
|[fcs-gx](https://github.com/ncbi/fcs-gx) | 10.1186/s13059-024-03198-7 | 0.5.5 | databases/fcs-gx.sif |
|[genomescope2](https://github.com/tbenavi1/genomescope2.0) | 10.1038/s41467-020-14998-3 | 2.0.1 | docker://diegomics/gep2_base:0.2 |
|[gfastats](https://github.com/vgl-hub/gfastats) | 10.1093/bioinformatics/btac460 | 1.3.11 | docker://diegomics/gep2_base:0.2 |
|[hicexplorer](https://github.com/deeptools/HiCExplorer) | 10.1093/gigascience/giac061 | 3.7.6 | docker://diegomics/hic_analysis:0.2 |
|[hifiadapterfilt](https://github.com/sheinasim-USDA/HiFiAdapterFilt) | 10.1186/s12864-022-08375-1 | 3.0.0 | docker://diegomics/gep2_base:0.2 |
|[hifasm](https://github.com/chhylp123/hifiasm) | 10.1038/s41592-024-02269-8 | 0.25.0 | docker://diegomics/gep2_base:0.2 |
|[inspector](https://github.com/Maggi-Chen/Inspector) | 10.1186/s13059-021-02527-4 | 1.3.1 | docker://diegomics/inspector:1.3.1 |
|[longdust](https://github.com/lh3/longdust) | - | 1.4 | docker://diegomics/hic_analysis:0.2 |
|[merqury](https://github.com/marbl/merqury) | 10.1186/s13059-020-02134-9 | 1.3 | docker://diegomics/gep2_base:0.2 |
|[merquryfk](https://github.com/thegenemyers/MERQURY.FK) | - | 1.2 | docker://diegomics/gep2_base:0.2 |
|[minimap2](https://github.com/lh3/minimap2) | 10.1093/bioinformatics/bty191 | 2.30 | docker://diegomics/hic_analysis:0.2 |
|[multiqc](https://github.com/MultiQC/MultiQC) | 10.1093/bioinformatics/btw354 | 1.33 | docker://diegomics/gep2_base:0.2 |
|[nanoplot](https://github.com/wdecoster/NanoPlot) | 10.1093/bioinformatics/btad311 | 1.46.2 | docker://diegomics/gep2_base:0.2 |
|[pairtools](https://github.com/open2c/pairtools) | 10.1101/2023.02.13.528389 | 1.1.3 | docker://diegomics/hic_analysis:0.2 |
|[pretextgraph](https://github.com/sanger-tol/PretextGraph) | - | 0.0.9 | docker://diegomics/hic_analysis:0.2 |
|[pretextmap](https://github.com/sanger-tol/PretextMap) | - | 0.2.4 | docker://diegomics/hic_analysis:0.2 |
|[pretextsnapshot](https://github.com/sanger-tol/PretextSnapshot) | - | 0.0.6 | docker://diegomics/hic_analysis:0.2 |
|[sambamba](https://github.com/biod/sambamba) | 10.1093/bioinformatics/btv098 | 1.0.1 | docker://diegomics/hic_analysis:0.2 |
|[samtools](https://github.com/samtools/samtools) | 10.1093/gigascience/giab008 | 1.22.1 | docker://diegomics/hic_analysis:0.2 |
|[sdust](https://github.com/lh3/sdust) | - | 0.1 | docker://diegomics/hic_analysis:0.2 |
|[seqkit](https://github.com/shenwei356/seqkit) | 10.1002/imt2.191 | 2.13.0 | docker://diegomics/gep2_base:0.2 |
|[seqtk](https://github.com/lh3/seqtk) | - | 1.5 | docker://diegomics/gep2_base:0.2 |
|[tidk](https://github.com/tolkit/telomeric-identifier) | 10.1093/bioinformatics/btaf049 | 0.2.65 | docker://diegomics/hic_analysis:0.2 |


### Outside the pipeline, you can run all of these programs just using the containers!
Once pulled, Snakemake saves the container images in a hidden folder using its md5 name. You can get the image name and save it in a variable, like this:
```
HIC_VAR=$(echo -n "docker://diegomics/hic_analysis:0.2" | md5sum | awk '{print $1}')
```

Then, you can get the GEP2 installation folder in another variable, like:
```
GEP2_FOLDER="/srv/public/users/ddepanis/Software/GEP2"
```

Next, you can combine both to get the full image path:
```
HIC_CONTAINER=${GEP2_FOLDER}/.snakemake/singularity/${HIC_VAR}.simg
```

Finally, you can run a given program included in the container like:
```
apptainer exec -B {path/to/data}:{path/to/data} $HIC_CONTAINER tidk -h

