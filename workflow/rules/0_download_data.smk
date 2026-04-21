# -------------------------------------------------------------------------------
# GEP2 - Download Data Rules
# -------------------------------------------------------------------------------


# -------------------------------------------------------------------------------
# LOAD DOWNLOAD MANIFEST
# -------------------------------------------------------------------------------

manifest_path = os.path.join(config["OUT_FOLDER"], "GEP2_results", "download_manifest.json")

if os.path.exists(manifest_path):
    with open(manifest_path) as f:
        DOWNLOAD_MANIFEST = json.load(f)
else:
    DOWNLOAD_MANIFEST = []


# -------------------------------------------------------------------------------
# HELPER FUNCTIONS
# -------------------------------------------------------------------------------

def get_manifest_entry(destination):
    """Get download manifest entry for a given destination path"""
    for entry in DOWNLOAD_MANIFEST:
        if entry["destination"] == destination:
            return entry
    return None


# -------------------------------------------------------------------------------
# WATCHDOG SETTINGS (hardcoded with generous headroom)
# -------------------------------------------------------------------------------

WATCHDOG_SCRIPT = str(BASEDIR / "scripts" / "watchdog.sh")


# -------------------------------------------------------------------------------
# RULE ORDER
# -------------------------------------------------------------------------------
# Prefer paired-end download over single-end when both patterns match
# (e.g., ERR123_1.fastq.gz can match {acc}_1.fastq.gz OR {acc}.fastq.gz)

ruleorder: _00_download_reads_sra > _00_download_reads_sra_single
ruleorder: _00_download_reads_sra_single > _00_download_reads_url
ruleorder: _00_download_reads_sra > _00_download_reads_url


# -------------------------------------------------------------------------------
# WILDCARD CONSTRAINTS (global for this module)
# -------------------------------------------------------------------------------
# SRA accessions are typically 3 letters + digits (ERR123456, SRR789, DRR001)
# The constraint prevents acc from containing underscores, so ERR123_1 won't match

wildcard_constraints:
    acc = r"[A-Z]{2,3}[0-9]+",
    species = r"[^/]+",
    read_type = r"[^/]+"


# -------------------------------------------------------------------------------
# RULES
# -------------------------------------------------------------------------------

rule _00_download_assembly:
    """Download assemblies from URLs or NCBI accessions"""
    output:
        asm = "{outdir}/downloaded_data/{species}/assemblies/{filename}"
    params:
        manifest = manifest_path
    threads: cpu_func("download_data")
    resources:
        mem_mb = mem_func("download_data"),
        runtime = time_func("download_data")
    container: CONTAINERS["gep2_base"]
    shell:
        """
        # Verify entry exists in manifest
        MANIFEST_INFO=$(python3 -c "
import json, sys
with open('{params.manifest}') as f:
    manifest = json.load(f)
for item in manifest:
    if item.get('type') == 'assembly' and item['destination'] == '{output.asm}':
        print(item['source'] + '|||' + item['method'])
        sys.exit(0)
print('')
sys.exit(0)
")
        
        if [ -z "$MANIFEST_INFO" ]; then
            echo "[GEP2] Error: No manifest entry found for {output.asm}"
            exit 1
        fi
        
        SOURCE=$(echo "$MANIFEST_INFO" | cut -d'|' -f1)
        METHOD=$(echo "$MANIFEST_INFO" | cut -d'|' -f4)
        
        mkdir -p $(dirname {output.asm})
        
        if [ "$METHOD" = "curl" ]; then
            echo "[GEP2] Downloading assembly from URL: $SOURCE"
            curl -L -C - --retry 3 --retry-delay 5 -o {output.asm}.tmp "$SOURCE"
            
            # Validate download
            if [ ! -s {output.asm}.tmp ]; then
                echo "[GEP2] Error: Downloaded file is empty"
                exit 1
            fi
            
            # Check minimum file size (10KB for assemblies)
            FILE_SIZE=$(stat -c%s "{output.asm}.tmp" 2>/dev/null || echo "0")
            if [ "$FILE_SIZE" -lt 10240 ]; then
                echo "[GEP2] Error: Downloaded file is suspiciously small ($FILE_SIZE bytes)"
                rm -f {output.asm}.tmp
                exit 1
            fi
            
            # Validate gzip integrity if compressed
            if [[ "{output.asm}" == *.gz ]]; then
                echo "[GEP2] Validating gzip integrity..."
                if ! gzip -t {output.asm}.tmp 2>/dev/null; then
                    echo "[GEP2] Error: Downloaded file failed gzip integrity check"
                    rm -f {output.asm}.tmp
                    exit 1
                fi
            fi
            
            mv {output.asm}.tmp {output.asm}
            echo "[GEP2] Downloaded: {output.asm}"
            
        elif [ "$METHOD" = "ncbi_assembly" ]; then
            echo "[GEP2] Downloading NCBI assembly: $SOURCE"
            
            # Parse accession using bash regex (e.g., GCA_963854735.1)
            if [[ $SOURCE =~ ^(GC[AF])_([0-9]{{3}})([0-9]{{3}})([0-9]{{3}})\\.([0-9]+)$ ]]; then
                PREFIX=${{BASH_REMATCH[1]}}
                P1=${{BASH_REMATCH[2]}}
                P2=${{BASH_REMATCH[3]}}
                P3=${{BASH_REMATCH[4]}}
                VERSION=${{BASH_REMATCH[5]}}
            else
                echo "[GEP2] Error: Invalid NCBI accession format: $SOURCE"
                exit 1
            fi
            
            # Build base FTP directory URL
            BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/${{PREFIX}}/${{P1}}/${{P2}}/${{P3}}"
            echo "[GEP2] Looking in: $BASE_URL"
            
            # Find the assembly directory (contains accession + assembly name)
            ASM_DIR=$(curl -sL "$BASE_URL/" | grep -oP "href=\\"${{SOURCE}}_[^/\\"]+" | head -1 | sed 's/href="//')
            
            if [ -z "$ASM_DIR" ]; then
                echo "[GEP2] Error: Could not find assembly directory for $SOURCE"
                exit 1
            fi
            
            # Construct full URL to genomic.fna.gz
            FULL_URL="${{BASE_URL}}/${{ASM_DIR}}/${{ASM_DIR}}_genomic.fna.gz"
            echo "[GEP2] Downloading from: $FULL_URL"
            
            curl -L -C - --retry 5 --retry-delay 10 -o {output.asm}.tmp "$FULL_URL"
            
            # Validate download
            if [ ! -s {output.asm}.tmp ]; then
                echo "[GEP2] Error: Downloaded file is empty"
                exit 1
            fi
            
            # Check minimum file size
            FILE_SIZE=$(stat -c%s "{output.asm}.tmp" 2>/dev/null || echo "0")
            if [ "$FILE_SIZE" -lt 10240 ]; then
                echo "[GEP2] Error: Downloaded file is suspiciously small ($FILE_SIZE bytes)"
                rm -f {output.asm}.tmp
                exit 1
            fi
            
            # Validate gzip file
            echo "[GEP2] Validating gzip integrity..."
            if ! gzip -t {output.asm}.tmp 2>/dev/null; then
                echo "[GEP2] Error: Downloaded file is not a valid gzip file"
                rm -f {output.asm}.tmp
                exit 1
            fi
            
            mv {output.asm}.tmp {output.asm}
            echo "[GEP2] Downloaded NCBI assembly: {output.asm}"
        else
            echo "[GEP2] Error: Unknown download method: $METHOD"
            exit 1
        fi
        """


rule _00_download_reads_sra_single:
    """Download single-end/long reads from SRA/ENA"""
    output:
        reads = "{outdir}/downloaded_data/{species}/reads/{read_type}/{acc}.fastq.gz"
    params:
        outdir = lambda w: os.path.join(w.outdir, "downloaded_data", w.species, "reads", w.read_type),
        manifest = manifest_path,
        watchdog = WATCHDOG_SCRIPT
    threads: cpu_func("download_data")
    resources:
        mem_mb = mem_func("download_data"),
        runtime = time_func("download_data")
    container: CONTAINERS["gep2_base"]
    shell:
        """
        # Verify accession exists in manifest as single-end/long reads
        python3 -c "
import json, sys
with open('{params.manifest}') as f:
    manifest = json.load(f)
found = False
for item in manifest:
    if (item.get('type') == 'reads' and 
        item.get('method') == 'enaDataGet' and 
        item.get('source') == '{wildcards.acc}'):
        if not item.get('paired', False):
            found = True
            break
if not found:
    print('[GEP2] Error: Accession {wildcards.acc} not found in manifest as single-end reads')
    sys.exit(1)
"
        
        source {params.watchdog}
        
        # CHECK IF DOWNLOAD PRODUCED FILES
        # -------------------------------------------------------------------
        check_single_files() {{
            local ACC=$1
            # Check for files in subdirectory or current directory
            # Single-end can come as ACC.fastq.gz, ACC.fastq, ACC_1.fastq.gz, or ACC_1.fastq
            if [ -d "$ACC" ]; then
                [ -f "$ACC/${{ACC}}.fastq.gz" ] || [ -f "$ACC/${{ACC}}.fastq" ] || \
                [ -f "$ACC/${{ACC}}_1.fastq.gz" ] || [ -f "$ACC/${{ACC}}_1.fastq" ]
            else
                [ -f "${{ACC}}.fastq.gz" ] || [ -f "${{ACC}}.fastq" ] || \
                [ -f "${{ACC}}_1.fastq.gz" ] || [ -f "${{ACC}}_1.fastq" ]
            fi
        }}
        
        # MAIN DOWNLOAD LOGIC
        # -------------------------------------------------------------------
        echo "[GEP2] Downloading single-end/long reads: {wildcards.acc}"
        
        mkdir -p {params.outdir}
        cd {params.outdir}
        
        MAX_RETRIES=3
        RETRY_DELAY=60
        USE_ASPERA=true  # Start with Aspera, switch to HTTP if it fails
        
        for ATTEMPT in $(seq 1 $MAX_RETRIES); do
            echo ""
            echo "[GEP2] ------------------------------------------------------------"
            echo "[GEP2] Download attempt $ATTEMPT of $MAX_RETRIES"
            echo "[GEP2] ------------------------------------------------------------"
            
            # Clean slate for this attempt
            rm -rf {wildcards.acc}/ {wildcards.acc}.fastq* {wildcards.acc}_*.fastq* 2>/dev/null || true
            
            # TRY ASPERA FIRST (if not already known to fail)
            # -----------------------------------------------------------------
            if [ "$USE_ASPERA" = "true" ]; then
                echo "[GEP2] Trying Aspera (fast) download..."
                
                # NOTE: `|| ASPERA_EXIT=$?` keeps the command off set -e's kill list.
                # Failure here is expected (fallback drives the retry loop).
                ASPERA_EXIT=0
                gep2_download_with_timeout 14400 enaDataGet.py -a -f fastq -d . {wildcards.acc} || ASPERA_EXIT=$?
                
                # Check if Aspera actually produced files (not just exit code!)
                if check_single_files {wildcards.acc}; then
                    echo "[GEP2] Aspera download produced files"
                else
                    echo "[GEP2] Aspera failed or produced no files (exit code: $ASPERA_EXIT)"
                    echo "[GEP2] Falling back to HTTP..."
                    USE_ASPERA=false
                    
                    # Clean up any partial Aspera artifacts
                    rm -rf {wildcards.acc}/ {wildcards.acc}.fastq* {wildcards.acc}_*.fastq* 2>/dev/null || true
                fi
            fi
            
            # TRY HTTP (if Aspera failed or was skipped)
            # -----------------------------------------------------------------
            if [ "$USE_ASPERA" = "false" ]; then
                echo "[GEP2] Using HTTP download..."
                
                HTTP_EXIT=0
                gep2_download_with_timeout 14400 enaDataGet.py -f fastq -d . {wildcards.acc} || HTTP_EXIT=$?
                
                echo "[GEP2] HTTP download finished (exit code: $HTTP_EXIT)"
            fi
            
            # TRY SUBMITTED FORMAT (if fastq format produced nothing)
            # -----------------------------------------------------------------
            if ! check_single_files {wildcards.acc}; then
                echo "[GEP2] No FASTQ files found, trying submitted format..."
                rm -rf {wildcards.acc}/ {wildcards.acc}.fastq* {wildcards.acc}_*.fastq* 2>/dev/null || true
                
                SUBMITTED_EXIT=0
                gep2_download_with_timeout 14400 enaDataGet.py -f submitted -d . {wildcards.acc} || SUBMITTED_EXIT=$?
                
                echo "[GEP2] Submitted files download finished (exit code: $SUBMITTED_EXIT)"
                
                # Rename unpredictable submitted file to standard naming
                DIR="."
                [ -d "{wildcards.acc}" ] && DIR="{wildcards.acc}"
                
                SUBMITTED=$(find "$DIR" -maxdepth 1 '(' -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" ')' 2>/dev/null | sort | head -1)
                
                if [ -n "$SUBMITTED" ]; then
                    if [[ "$SUBMITTED" == *.gz ]]; then
                        TARGET="$DIR/{wildcards.acc}.fastq.gz"
                    else
                        TARGET="$DIR/{wildcards.acc}.fastq"
                    fi
                    
                    echo "[GEP2] Renaming submitted file to standard naming..."
                    echo "[GEP2]   $SUBMITTED -> $TARGET"
                    mv "$SUBMITTED" "$TARGET"
                else
                    echo "[GEP2] Warning: No submitted fastq files found"
                fi
            fi

            # TRY ENA PORTAL API + HTTPS (last-resort fallback)
            # -----------------------------------------------------------------
            # Handles submissions where enaDataGet.py returns "no files" or
            # crashes, e.g., submissions with FASTQ only in submitted_ftp, not fastq_ftp.
            if ! check_single_files {wildcards.acc}; then
                echo "[GEP2] Trying ENA portal API + HTTPS as last resort..."
                rm -rf {wildcards.acc}/ {wildcards.acc}.fastq* {wildcards.acc}_*.fastq* 2>/dev/null || true
                
                API_EXIT=0
                gep2_ena_download_single {wildcards.acc} || API_EXIT=$?
                
                if [ "$API_EXIT" -eq 0 ] && check_single_files {wildcards.acc}; then
                    echo "[GEP2] Portal-API download produced files"
                else
                    echo "[GEP2] Portal-API fallback failed (exit $API_EXIT)"
                fi
            fi

            # PROCESS AND VALIDATE FILES
            # -----------------------------------------------------------------
            
            # Move files from subdirectory if created
            if [ -d "{wildcards.acc}" ]; then
                echo "[GEP2] Moving files from subdirectory..."
                mv {wildcards.acc}/* . 2>/dev/null || true
                rmdir {wildcards.acc} 2>/dev/null || true
            fi
            
            # Handle different naming conventions from ENA
            # Single-end sometimes comes as _1.fastq even without a _2
            if [ -f "{wildcards.acc}_1.fastq" ] && [ ! -f "{wildcards.acc}_2.fastq" ]; then
                echo "[GEP2] Renaming _1.fastq to .fastq (single-end file)..."
                mv "{wildcards.acc}_1.fastq" "{wildcards.acc}.fastq"
            elif [ -f "{wildcards.acc}_1.fastq.gz" ] && [ ! -f "{wildcards.acc}_2.fastq.gz" ]; then
                echo "[GEP2] Renaming _1.fastq.gz to .fastq.gz (single-end file)..."
                mv "{wildcards.acc}_1.fastq.gz" "{wildcards.acc}.fastq.gz"
            fi
            
            # Compress if needed
            if [ -f "{wildcards.acc}.fastq" ]; then
                echo "[GEP2] Compressing..."
                pigz -p {threads} "{wildcards.acc}.fastq"
            fi
            
            # Clean up any unexpected _2 file (shouldn't exist for single-end)
            rm -f "{wildcards.acc}_2.fastq.gz" "{wildcards.acc}_2.fastq" 2>/dev/null || true
            
            # Check if we got the file
            if [ -f "{output.reads}" ]; then
                echo "[GEP2] Validating downloaded file..."
                
                FILE_SIZE=$(stat -c%s "{output.reads}" 2>/dev/null || echo "0")
                
                echo "[GEP2] File size: $FILE_SIZE bytes"
                
                if [ "$FILE_SIZE" -lt 1024 ]; then
                    echo "[GEP2] Downloaded file is suspiciously small"
                else
                    if gzip -t "{output.reads}" 2>/dev/null; then
                        echo "[GEP2] Downloaded and validated: {wildcards.acc}"
                        exit 0
                    else
                        echo "[GEP2] Downloaded file failed gzip integrity check"
                    fi
                fi
            else
                echo "[GEP2] Expected file not found after download"
                echo "[GEP2] Looking for: {output.reads}"
                echo "[GEP2] Directory contents:"
                ls -la {params.outdir}/ 2>/dev/null || echo "(empty)"
            fi
            
            # RETRY LOGIC
            # -----------------------------------------------------------------
            if [ $ATTEMPT -lt $MAX_RETRIES ]; then
                echo ""
                echo "[GEP2] Download attempt $ATTEMPT failed"
                echo "[GEP2] Retrying in $RETRY_DELAY seconds..."
                sleep $RETRY_DELAY
                rm -rf {wildcards.acc}/ {wildcards.acc}.fastq* {wildcards.acc}_*.fastq* 2>/dev/null || true
            fi
        done
        
        # All retries exhausted
        echo ""
        echo "[GEP2] ------------------------------------------------------------"
        echo "[GEP2] ENA DOWNLOAD FAILED AFTER $MAX_RETRIES ATTEMPTS"
        echo "[GEP2] ------------------------------------------------------------"
        echo "[GEP2]"
        echo "[GEP2] This could be due to:"
        echo "[GEP2]   - ENA/EBI server issues"
        echo "[GEP2]   - Network problems"
        echo "[GEP2]   - Download repeatedly stalling (killed by watchdog)"
        echo "[GEP2]"
        echo "[GEP2] Please verify the accession at:"
        echo "[GEP2]   https://www.ebi.ac.uk/ena/browser/view/{wildcards.acc}"
        echo "[GEP2]"
        echo "[GEP2] Contents of output directory:"
        ls -lh {params.outdir}/ 2>/dev/null || echo "[GEP2] (directory empty or not found)"
        echo ""
        exit 1
        """


rule _00_download_reads_sra:
    """Download paired-end reads from SRA/ENA"""
    output:
        r1 = "{outdir}/downloaded_data/{species}/reads/{read_type}/{acc}_1.fastq.gz",
        r2 = "{outdir}/downloaded_data/{species}/reads/{read_type}/{acc}_2.fastq.gz"
    params:
        outdir = lambda w: os.path.join(w.outdir, "downloaded_data", w.species, "reads", w.read_type),
        manifest = manifest_path,
        watchdog = WATCHDOG_SCRIPT
    threads: cpu_func("download_data")
    resources:
        mem_mb = mem_func("download_data"),
        runtime = time_func("download_data")
    container: CONTAINERS["gep2_base"]
    shell:
        """
        # Verify accession exists in manifest as paired-end reads
        python3 -c "
import json, sys
with open('{params.manifest}') as f:
    manifest = json.load(f)
found = False
for item in manifest:
    if (item.get('type') == 'reads' and 
        item.get('method') == 'enaDataGet' and 
        item.get('paired') == True and
        item.get('source') == '{wildcards.acc}'):
        found = True
        break
if not found:
    print('[GEP2] Error: Accession {wildcards.acc} not found in manifest as paired-end reads')
    sys.exit(1)
"
        
        source {params.watchdog}

        # CHECK IF DOWNLOAD PRODUCED PAIRED FILES
        # -------------------------------------------------------------------
        check_paired_files() {{
            local ACC=$1
            if [ -d "$ACC" ]; then
                {{ [ -f "$ACC/${{ACC}}_1.fastq.gz" ] || [ -f "$ACC/${{ACC}}_1.fastq" ]; }} && \
                {{ [ -f "$ACC/${{ACC}}_2.fastq.gz" ] || [ -f "$ACC/${{ACC}}_2.fastq" ]; }}
            else
                {{ [ -f "${{ACC}}_1.fastq.gz" ] || [ -f "${{ACC}}_1.fastq" ]; }} && \
                {{ [ -f "${{ACC}}_2.fastq.gz" ] || [ -f "${{ACC}}_2.fastq" ]; }}
            fi
        }}

        
        # MAIN DOWNLOAD LOGIC
        # -------------------------------------------------------------------
        echo "[GEP2] Downloading paired-end reads: {wildcards.acc}"
        
        mkdir -p {params.outdir}
        cd {params.outdir}
        
        MAX_RETRIES=3
        RETRY_DELAY=60
        USE_ASPERA=true  # Start with Aspera, switch to HTTP if it fails
        
        for ATTEMPT in $(seq 1 $MAX_RETRIES); do
            echo ""
            echo "[GEP2] ------------------------------------------------------------"
            echo "[GEP2] Download attempt $ATTEMPT of $MAX_RETRIES"
            echo "[GEP2] ------------------------------------------------------------"
            
            # Clean slate for this attempt
            rm -rf {wildcards.acc}/ {wildcards.acc}_*.fastq* 2>/dev/null || true
            
            # TRY ASPERA FIRST (if not already known to fail)
            # -----------------------------------------------------------------
            if [ "$USE_ASPERA" = "true" ]; then
                echo "[GEP2] Trying Aspera (fast) download..."
                
                # NOTE: `|| ASPERA_EXIT=$?` keeps the command off set -e's kill list.
                # Failure here is expected (fallback drives the retry loop).
                ASPERA_EXIT=0
                gep2_download_with_timeout 14400 enaDataGet.py -a -f fastq -d . {wildcards.acc} || ASPERA_EXIT=$?
                
                # Check if Aspera actually produced files (not just exit code!)
                if check_paired_files {wildcards.acc}; then
                    echo "[GEP2] Aspera download produced files"
                else
                    echo "[GEP2] Aspera failed or produced no files (exit code: $ASPERA_EXIT)"
                    echo "[GEP2] Falling back to HTTP..."
                    USE_ASPERA=false
                    
                    # Clean up any partial Aspera artifacts
                    rm -rf {wildcards.acc}/ {wildcards.acc}_*.fastq* 2>/dev/null || true
                fi
            fi
            
            # TRY HTTP (if Aspera failed or was skipped)
            # -----------------------------------------------------------------
            if [ "$USE_ASPERA" = "false" ]; then
                echo "[GEP2] Using HTTP download..."
                
                HTTP_EXIT=0
                gep2_download_with_timeout 14400 enaDataGet.py -f fastq -d . {wildcards.acc} || HTTP_EXIT=$?
                
                echo "[GEP2] HTTP download finished (exit code: $HTTP_EXIT)"
            fi
            
            # TRY SUBMITTED FORMAT (if fastq format produced nothing)
            # -----------------------------------------------------------------
            if ! check_paired_files {wildcards.acc}; then
                echo "[GEP2] No FASTQ files found, trying submitted format..."
                rm -rf {wildcards.acc}/ {wildcards.acc}_*.fastq* 2>/dev/null || true
                
                SUBMITTED_EXIT=0
                gep2_download_with_timeout 14400 enaDataGet.py -f submitted -d . {wildcards.acc} || SUBMITTED_EXIT=$?
                
                echo "[GEP2] Submitted files download finished (exit code: $SUBMITTED_EXIT)"

                # Rename unpredictable submitted files to standard naming
                DIR="."
                [ -d "{wildcards.acc}" ] && DIR="{wildcards.acc}"
                
                FASTQ_FILES=$(find "$DIR" -maxdepth 1 '(' -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" ')' 2>/dev/null | sort)
                NUM_FILES=$(echo "$FASTQ_FILES" | grep -c . || true)
                
                if [ "$NUM_FILES" -ge 2 ]; then
                    R1=$(echo "$FASTQ_FILES" | sed -n '1p')
                    R2=$(echo "$FASTQ_FILES" | sed -n '2p')
                    
                    # Detect extension to avoid renaming .fastq as .fastq.gz
                    EXT1="${{R1##*.fastq}}"
                    if [[ "$R1" == *.gz ]]; then
                        TARGET_R1="$DIR/{wildcards.acc}_1.fastq.gz"
                        TARGET_R2="$DIR/{wildcards.acc}_2.fastq.gz"
                    else
                        TARGET_R1="$DIR/{wildcards.acc}_1.fastq"
                        TARGET_R2="$DIR/{wildcards.acc}_2.fastq"
                    fi
                    
                    echo "[GEP2] Renaming submitted files to standard naming..."
                    echo "[GEP2]   $R1 -> $TARGET_R1"
                    echo "[GEP2]   $R2 -> $TARGET_R2"
                    mv "$R1" "$TARGET_R1"
                    mv "$R2" "$TARGET_R2"
                else
                    echo "[GEP2] Warning: Could not find two submitted fastq files to pair"
                    echo "[GEP2] Files found:"
                    echo "$FASTQ_FILES"
                fi
            fi

            # TRY ENA PORTAL API + HTTPS (last-resort fallback)
            # -----------------------------------------------------------------
             # Handles submissions where enaDataGet.py returns "no files" or
            # crashes, e.g., submissions with FASTQ only in submitted_ftp, not fastq_ftp.
            if ! check_paired_files {wildcards.acc}; then
                echo "[GEP2] Trying ENA portal API + HTTPS as last resort..."
                rm -rf {wildcards.acc}/ {wildcards.acc}_*.fastq* 2>/dev/null || true
                
                API_EXIT=0
                gep2_ena_download_paired {wildcards.acc} || API_EXIT=$?
                
                if [ "$API_EXIT" -eq 0 ] && check_paired_files {wildcards.acc}; then
                    echo "[GEP2] Portal-API download produced files"
                else
                    echo "[GEP2] Portal-API fallback failed (exit $API_EXIT)"
                fi
            fi

            # PROCESS AND VALIDATE FILES
            # -----------------------------------------------------------------

            # Move files from subdirectory if created (fastq format downloads)
            if [ -d "{wildcards.acc}" ]; then
                echo "[GEP2] Moving files from subdirectory..."
                mv {wildcards.acc}/* . 2>/dev/null || true
                rmdir {wildcards.acc} 2>/dev/null || true
            fi
            
            # Compress if needed
            if [ -f "{wildcards.acc}_1.fastq" ]; then
                echo "[GEP2] Compressing R1..."
                pigz -p {threads} "{wildcards.acc}_1.fastq"
            fi
            if [ -f "{wildcards.acc}_2.fastq" ]; then
                echo "[GEP2] Compressing R2..."
                pigz -p {threads} "{wildcards.acc}_2.fastq"
            fi
            
            # Check if we got both files
            if [ -f "{output.r1}" ] && [ -f "{output.r2}" ]; then
                echo "[GEP2] Validating downloaded files..."
                
                R1_SIZE=$(stat -c%s "{output.r1}" 2>/dev/null || echo "0")
                R2_SIZE=$(stat -c%s "{output.r2}" 2>/dev/null || echo "0")
                
                echo "[GEP2] File sizes: R1=$R1_SIZE bytes, R2=$R2_SIZE bytes"
                
                if [ "$R1_SIZE" -lt 1024 ] || [ "$R2_SIZE" -lt 1024 ]; then
                    echo "[GEP2] Downloaded files are suspiciously small"
                else
                    if gzip -t "{output.r1}" 2>/dev/null && gzip -t "{output.r2}" 2>/dev/null; then
                        echo "[GEP2] Downloaded and validated paired reads: {wildcards.acc}"
                        exit 0
                    else
                        echo "[GEP2] Downloaded files failed gzip integrity check"
                    fi
                fi
            else
                echo "[GEP2] Expected files not found after download"
                echo "[GEP2] Looking for: {output.r1}"
                echo "[GEP2]             {output.r2}"
                echo "[GEP2] Directory contents:"
                ls -la {params.outdir}/ 2>/dev/null || echo "(empty)"
            fi
            
            # RETRY LOGIC
            # -----------------------------------------------------------------
            if [ $ATTEMPT -lt $MAX_RETRIES ]; then
                echo ""
                echo "[GEP2] Download attempt $ATTEMPT failed"
                echo "[GEP2] Retrying in $RETRY_DELAY seconds..."
                sleep $RETRY_DELAY
                rm -rf {wildcards.acc}/ {wildcards.acc}_*.fastq* 2>/dev/null || true
            fi
        done
        
        # All retries exhausted
        echo ""
        echo "[GEP2] ------------------------------------------------------------"
        echo "[GEP2] ENA DOWNLOAD FAILED AFTER $MAX_RETRIES ATTEMPTS"
        echo "[GEP2] ------------------------------------------------------------"
        echo "[GEP2]"
        echo "[GEP2] This could be due to:"
        echo "[GEP2]   - ENA/EBI server issues"
        echo "[GEP2]   - Network problems"
        echo "[GEP2]   - Download repeatedly stalling (killed by watchdog)"
        echo "[GEP2]"
        echo "[GEP2] Please verify the accession at:"
        echo "[GEP2]   https://www.ebi.ac.uk/ena/browser/view/{wildcards.acc}"
        echo "[GEP2]"
        echo "[GEP2] Contents of output directory:"
        ls -lh {params.outdir}/ 2>/dev/null || echo "[GEP2] (directory empty or not found)"
        echo ""
        exit 1
        """


rule _00_download_reads_url:
    """Download reads from direct URLs (only matches non-SRA filenames via url_filename constraint)"""
    output:
        reads = "{outdir}/downloaded_data/{species}/reads/{read_type}/{url_filename}"
    params:
        manifest = manifest_path
    threads: cpu_func("download_data")
    resources:
        mem_mb = mem_func("download_data"),
        runtime = time_func("download_data")
    container: CONTAINERS["gep2_base"]
    shell:
        """
        SOURCE=$(python3 -c "
import json, sys
with open('{params.manifest}') as f:
    manifest = json.load(f)
for item in manifest:
    if item.get('type') == 'reads' and item.get('method') == 'curl' and item['destination'] == '{output.reads}':
        print(item['source'])
        sys.exit(0)
print('')
sys.exit(0)
")
        
        if [ -z "$SOURCE" ]; then
            echo "[GEP2] Error: No URL source found in manifest for {output.reads}"
            exit 1
        fi
        
        mkdir -p $(dirname {output.reads})
        
        echo "[GEP2] Downloading reads from URL: $SOURCE"
        curl -L -C - --retry 3 --retry-delay 5 -o {output.reads}.tmp "$SOURCE"
        
        # Validate download
        if [ ! -s {output.reads}.tmp ]; then
            echo "[GEP2] Error: Downloaded file is empty"
            exit 1
        fi
        
        # Check minimum file size (1KB)
        FILE_SIZE=$(stat -c%s "{output.reads}.tmp" 2>/dev/null || echo "0")
        if [ "$FILE_SIZE" -lt 1024 ]; then
            echo "[GEP2] Error: Downloaded file is suspiciously small ($FILE_SIZE bytes)"
            rm -f {output.reads}.tmp
            exit 1
        fi
        
        # Validate gzip integrity for compressed files
        if [[ "{output.reads}" == *.gz ]]; then
            echo "[GEP2] Validating gzip integrity..."
            if ! gzip -t {output.reads}.tmp 2>/dev/null; then
                echo "[GEP2] Error: Downloaded file failed gzip integrity check"
                rm -f {output.reads}.tmp
                exit 1
            fi
        fi
        
        mv {output.reads}.tmp {output.reads}
        echo "[GEP2] Downloaded reads: {output.reads}"
        """