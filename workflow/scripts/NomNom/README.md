# <span dir="ltr">&nbsp;(∪• ﻌ •∪)</span><br>NomNom

_NomNom_ is a little Python tool that can ingest a genomic data table and convert it into tidy, hierarchical YAML file. It uses data-driven methods to infer table structure and produces organised, consistent outputs for genomic analysis pipelines.

**Index:**<br>
1.&nbsp;Features<br>
2.&nbsp;Usage<br>
3.&nbsp;Attributes file<br>

---

### 1. Features

#### File Reading & Parsing
- **Multi-format support**: Excel (.xlsx, .xls), LibreOffice/OpenOffice (.ods), Apple Numbers, Google Sheets exports, CSV, TSV, and other delimited formats
- **Auto-delimiter detection**: comma, semicolon, tab, pipe, backslash
- **Header fixing**: converts multiple spaces to proper delimiters
- **Whitespace cleaning**: removes extra spaces from headers and data
- **Robust error handling**: graceful fallbacks for problematic file formats
- **Binary format protection**: prevents incorrect YAML detection of Excel/binary files

#### Type Inference
- **Boolean detection**: true/false, yes/no, y/n, 1/0 patterns
- **Integer detection**: numeric values with special handling for large numbers and 'auto' values
- **Path detection**: distinguishes file paths (simple vs. complex paired paths)
- **Category detection**: columns with limited unique values
- **Column name hints**: uses naming patterns to guide type inference
- **Priority logic**: integers checked before booleans to prevent misclassification
- **Context-aware classification**: distinguishes between data types and structural categories

#### Paired Column Detection & Merging
- **Automatic pairing**: detects R1/R2, Forward/Reverse, _1/_2 patterns
- **Pattern varieties**: _R1/_R2, _1/_2, _F/_R, .R1/.R2, Forward/Reverse, Left/Right
- **Special cases**: standalone pairs like "R1"/"R2" become "Paired_x"
- **Suffix preservation**: creates appropriate merged column names (Library_Rx, HiC_x, etc.)
- **Comma-separated storage**: paired values stored as "file1.fq, file2.fq"
- **Enhanced pattern matching**: supports complex genomic naming conventions

#### Data-Driven Hierarchy Detection
- **Pattern analysis**: analyzes how values repeat and vary across rows
- **Content-based rules**: excludes file paths, boolean flags, high-cardinality columns
- **Grouping effectiveness**: tests if columns create meaningful sub-groups
- **Name-based hints**: uses column names as secondary indicators
- **Cardinality analysis**: identifies good grouping variables (2-20 unique values typically)
- **Property vs. hierarchy distinction**: separates data columns from structural columns
- **Simple table optimization**: special handling for tables with ≤5 columns
- **Hierarchical exceptions**: recognizes when columns like "Read type" should be hierarchical despite containing property keywords
- **Enhanced hierarchy detection**: optimized for simple genomic tables (Species → Read type → Properties)

#### Path Handling
- **Single paths**: `Path1: /path/to/file`
- **Paired reads**: `Path1: /path/R1.fq, /path/R2.fq`
- **Multiple libraries**: `Path1: libA_1.fq, libA_2.fq; Path2: libB_1.fq, libB_2.fq`
- **Accession number support**: recognizes and properly handles SRR*, ERR*, DRR*, GCA_*, GCF_* accessions as individual paths
- **Auto-pairing detection**: recognizes paired patterns even without read type info
- **Flexible pairing**: handles 2, 4, 6+ files with smart matching
- **Consistent formatting**: all file references use Path1, Path2, etc. structure
- **Context-aware path detection**: distinguishes actual file paths from categorical data
- **Glob pattern expansion**: `--discover-paths` flag automatically expands wildcards (*, ?, []) to actual file paths
- **Path validation**: `--validate-paths` flag ensures all referenced files exist before generating YAML

#### Hierarchical Structure Building
- **Left-to-right order**: respects column order for hierarchy levels
- **Labeled levels**: creates explicit labels for categorical columns (sample_id:, read_type:)
- **Semantic property placement**: places data at appropriate hierarchy levels based on biological meaning
- **Multi-row handling**: correctly groups multiple rows with same hierarchy values
- **Consistent structure**: ensures all entries have identical field structure
- **Assembly-aware organization**: recognizes assembly-level vs. sample-level vs. read-level properties

#### Attributes System
- **Auto-generation**: creates .attributes files with column metadata
- **Simple format**: e.g. `sp_name: category (missing = None); mandatory; sanitise`
- **Type validation**: validates data against defined constraints
- **Default values**: provides appropriate defaults for missing data
- **Preservation**: maintains existing attributes while adding new ones
- **Mandatory field detection**: automatically identifies required columns
- **Flexible overriding**: `--overwrite-attributes` flag to regenerate using current smart detection

#### Bi-Directional Conversion
- **Table → YAML**: creates hierarchical structure from flat tables
- **YAML → Table**: flattens hierarchical YAML back to tabular format
- **Type preservation**: maintains data types in both directions
- **Structure validation**: ensures consistent schema throughout conversion
- **Format auto-detection**: automatically determines input type (table vs. YAML)

#### Data Processing
- **Consistent defaults**: all samples get same fields with appropriate default values
- **Error handling**: comprehensive validation with helpful error messages
- **Warning system**: alerts for data issues without stopping processing
- **Clean output**: eliminates YAML anchors/references for clean, readable files
- **Reference-free YAML**: produces clean YAML without anchors or aliases for maximum compatibility
- **Overwrite protection**: prevents accidental file overwrites unless `--overwrite-yaml` flag is used

---

### 2. Usage

#### Dependencies
The dependencies for running NomNom are already included in the `game_base` container. To run it in the container:
```bash
apptainer exec docker://diegomics/game_base:0.1 python NomNom.py -h
```

Alternatively, create a standalone conda environment:
```bash
conda env create -f NOM_env.yml
conda activate NOM_env
```

#### Basic Usage

**Convert table to YAML:**
```bash
python NomNom.py my_table.csv
```

**Convert YAML back to table:** #better testing needed
```bash
python NomNom.py my_data.yaml
```

#### Common Options

**Use custom attributes file:**
```bash
python NomNom.py -a my_custom.attributes my_table.csv
```

**Specify output path:**
```bash
python NomNom.py -o output.yaml my_table.csv
```

**Validate that all file paths exist:**
```bash
python NomNom.py --validate-paths my_table.csv
```

**Expand glob patterns in path columns:**
```bash
python NomNom.py --discover-paths my_table.csv
```

**Sanitize hierarchical field values (replace spaces/special chars with underscores):**
```bash
python NomNom.py --sanitise-fields my_table.csv
```

**Allow overwriting existing YAML file:**
```bash
python NomNom.py --overwrite-yaml my_table.csv
```

**Regenerate attributes file ignoring existing one:**
```bash
python NomNom.py --overwrite-attributes my_table.csv
```

**Quiet mode (suppress console output, only save files):**
```bash
python NomNom.py --quiet my_table.csv
```

#### Combined Example
```bash
python NomNom.py \
  -a GAME_table.attributes \
  -o samples.yaml \
  --discover-paths \
  --validate-paths \
  my_samples.csv
```

This will:
1. Use the custom attributes file to guide type inference
2. Expand any glob patterns (like `*.fastq`) in path columns
3. Validate that all referenced files actually exist
4. Save the output to `samples.yaml`

---

### 3. Attributes File

The attributes file controls how NomNom interprets your data. It's a simple text file with one line per column:

#### Format
```
column_name: type (missing = default); mandatory
```

#### Components
- **column_name**: exact name from your table header
- **type**: one of `string`, `integer`, `boolean`, `path`, `simple_path`, `complex_path`, `category`
- **default**: value to use when field is empty (optional)
- **mandatory**: marks column as required (optional)

#### Type Descriptions
- `string`: text data
- `integer`: whole numbers
- `boolean`: true/false values
- `path`: file paths (auto-detects single vs. multiple, handles accessions)
- `simple_path`: single file path per entry
- `complex_path`: multiple or paired file paths
- `category`: categorical data with limited unique values

#### Example Attributes File
```
sp_name: string; mandatory
asm_id: string; mandatory
asm_file: simple_path (missing = None); mandatory
read_type: category; mandatory
read_files: path (missing = None)
trim: boolean (missing = False)
genome_size: integer (missing = auto)
kmer_size: integer (missing = auto)
```

#### Special Flags
- `hierarchical`: force column to be treated as hierarchy level
- `property`: force column to be treated as property (not hierarchy)
- `sanitise`: replace spaces and special characters with underscores in this column's values

#### Auto-Generation
If no attributes file exists, NomNom will automatically create one (`{input_file}.attributes`) using smart detection. You can then edit it to refine the behavior.

To regenerate an attributes file using current smart detection:
```bash
python NomNom.py --overwrite-attributes my_table.csv
```

---

### Algorithm Flow

1. **File Reading**: Auto-detect format, fix headers, clean whitespace
2. **Column Merging**: Detect and merge paired columns with _x suffix
3. **Type Inference**: Analyze content to determine appropriate data types
4. **Hierarchy Detection**: Identify which columns form grouping structure vs. properties
5. **Structure Building**: Create hierarchical organization with semantic property placement
6. **Consistency Enforcement**: Ensure all entries have identical field structure
7. **Output Generation**: Clean data and export to YAML with proper formatting

