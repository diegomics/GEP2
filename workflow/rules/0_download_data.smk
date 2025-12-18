# ═══════════════════════════════════════════════════════════════════════════════
# GEP2 - Download Data Rules
# ═══════════════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════════════
# LOAD DOWNLOAD MANIFEST
# ═══════════════════════════════════════════════════════════════════════════════

manifest_path = os.path.join(config["OUT_FOLDER"], "GEP2_results", "download_manifest.json")

if os.path.exists(manifest_path):
    with open(manifest_path) as f:
        DOWNLOAD_MANIFEST = json.load(f)
else:
    DOWNLOAD_MANIFEST = []


# ═══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def get_manifest_entry(destination):
    """Get download manifest entry for a given destination path"""
    for entry in DOWNLOAD_MANIFEST:
        if entry["destination"] == destination:
            return entry
    return None


# ═══════════════════════════════════════════════════════════════════════════════
# RULE ORDER
# ═══════════════════════════════════════════════════════════════════════════════
# Prefer paired-end download over single-end when both patterns match
# (e.g., ERR123_1.fastq.gz can match {acc}_1.fastq.gz OR {acc}.fastq.gz)

ruleorder: _00_download_reads_sra > _00_download_reads_sra_single
ruleorder: _00_download_reads_sra_single > _00_download_reads_url
ruleorder: _00_download_reads_sra > _00_download_reads_url


# ═══════════════════════════════════════════════════════════════════════════════
# WILDCARD CONSTRAINTS (global for this module)
# ═══════════════════════════════════════════════════════════════════════════════
# SRA accessions are typically 3 letters + digits (ERR123456, SRR789, DRR001)
# The constraint prevents acc from containing underscores, so ERR123_1 won't match

wildcard_constraints:
    acc = r"[A-Z]{2,3}[0-9]+",
    species = r"[^/]+",
    read_type = r"[^/]+"


# ═══════════════════════════════════════════════════════════════════════════════
# RULES
# ═══════════════════════════════════════════════════════════════════════════════

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
        # Read manifest to get source and method
        MANIFEST_INFO=$(python3 -c "
import json
with open('{params.manifest}') as f:
    manifest = json.load(f)
for item in manifest:
    if item.get('type') == 'assembly' and item['destination'] == '{output.asm}':
        print(item['source'] + '|||' + item['method'])
        break
")
        
        if [ -z \"$MANIFEST_INFO\" ]; then
            echo \"Error: No manifest entry found for {output.asm}\"
            exit 1
        fi
        
        SOURCE=$(echo \"$MANIFEST_INFO\" | cut -d'|' -f1)
        METHOD=$(echo \"$MANIFEST_INFO\" | cut -d'|' -f4)
        
        mkdir -p $(dirname {output.asm})
        
        if [ \"$METHOD\" = \"curl\" ]; then
            echo \"[GEP2] Downloading assembly from URL: $SOURCE\"
            curl -L -C - --retry 3 --retry-delay 5 -o {output.asm}.tmp \"$SOURCE\"
            
            if [ ! -s {output.asm}.tmp ]; then
                echo \"Error: Downloaded file is empty\"
                exit 1
            fi
            
            mv {output.asm}.tmp {output.asm}
            echo \"[GEP2] ✅ Downloaded: {output.asm}\"
            
        elif [ \"$METHOD\" = \"ncbi_assembly\" ]; then
            echo \"[GEP2] Downloading NCBI assembly: $SOURCE\"
            
            # Parse accession using bash regex (e.g., GCA_963854735.1)
            if [[ $SOURCE =~ ^(GC[AF])_([0-9]{{3}})([0-9]{{3}})([0-9]{{3}})\\.([0-9]+)$ ]]; then
                PREFIX=${{BASH_REMATCH[1]}}
                P1=${{BASH_REMATCH[2]}}
                P2=${{BASH_REMATCH[3]}}
                P3=${{BASH_REMATCH[4]}}
                VERSION=${{BASH_REMATCH[5]}}
            else
                echo \"Error: Invalid NCBI accession format: $SOURCE\"
                exit 1
            fi
            
            # Build base FTP directory URL
            BASE_URL=\"https://ftp.ncbi.nlm.nih.gov/genomes/all/${{PREFIX}}/${{P1}}/${{P2}}/${{P3}}\"
            echo \"[GEP2] Looking in: $BASE_URL\"
            
            # Find the assembly directory (contains accession + assembly name)
            ASM_DIR=$(curl -sL \"$BASE_URL/\" | grep -oP \"href=\\\"${{SOURCE}}_[^/\\\"]+\" | head -1 | sed 's/href=\"//')
            
            if [ -z \"$ASM_DIR\" ]; then
                echo \"Error: Could not find assembly directory for $SOURCE\"
                exit 1
            fi
            
            # Construct full URL to genomic.fna.gz
            FULL_URL=\"${{BASE_URL}}/${{ASM_DIR}}/${{ASM_DIR}}_genomic.fna.gz\"
            echo \"[GEP2] Downloading from: $FULL_URL\"
            
            curl -L -C - --retry 5 --retry-delay 10 -o {output.asm}.tmp \"$FULL_URL\"
            
            if [ ! -s {output.asm}.tmp ]; then
                echo \"Error: Downloaded file is empty\"
                exit 1
            fi
            
            # Validate gzip file
            if ! gzip -t {output.asm}.tmp 2>/dev/null; then
                echo \"Error: Downloaded file is not a valid gzip file\"
                exit 1
            fi
            
            mv {output.asm}.tmp {output.asm}
            echo \"[GEP2] ✅ Downloaded NCBI assembly: {output.asm}\"
        else
            echo \"Error: Unknown download method: $METHOD\"
            exit 1
        fi
        """


rule _00_download_reads_sra_single:
    """Download single-end/long reads from SRA/ENA"""
    output:
        reads = "{outdir}/downloaded_data/{species}/reads/{read_type}/{acc}.fastq.gz"
    params:
        outdir = lambda w: os.path.join(w.outdir, "downloaded_data", w.species, "reads", w.read_type)
    threads: cpu_func("download_data")
    resources:
        mem_mb = mem_func("download_data"),
        runtime = time_func("download_data")
    container: CONTAINERS["gep2_base"]
    shell:
        """
        echo "[GEP2] Downloading single-end/long reads: {wildcards.acc}"
        
        mkdir -p {params.outdir}
        cd {params.outdir}
        
        MAX_RETRIES=3
        RETRY_DELAY=60
        
        for ATTEMPT in $(seq 1 $MAX_RETRIES); do
            echo "[GEP2] Download attempt $ATTEMPT of $MAX_RETRIES..."
            
            if enaDataGet.py -f fastq -d . {wildcards.acc}; then
                # Move files from subdirectory if created
                if [ -d "{wildcards.acc}" ]; then
                    mv {wildcards.acc}/* . 2>/dev/null || true
                    rmdir {wildcards.acc} 2>/dev/null || true
                fi
                
                # Compress and rename as needed
                if [ -f "{wildcards.acc}.fastq" ]; then
                    pigz -p {threads} "{wildcards.acc}.fastq"
                elif [ -f "{wildcards.acc}_1.fastq" ]; then
                    pigz -p {threads} "{wildcards.acc}_1.fastq"
                    mv "{wildcards.acc}_1.fastq.gz" "{wildcards.acc}.fastq.gz"
                elif [ -f "{wildcards.acc}_1.fastq.gz" ]; then
                    mv "{wildcards.acc}_1.fastq.gz" "{wildcards.acc}.fastq.gz"
                fi
                
                # Check if we got the file
                if [ -f "{output.reads}" ]; then
                    rm -f "{params.outdir}/{wildcards.acc}_2.fastq.gz"
                    rm -f "{params.outdir}/{wildcards.acc}_2.fastq"
                    echo "[GEP2] ✅ Downloaded: {wildcards.acc}"
                    exit 0
                fi
            fi
            
            # If we get here, download failed
            if [ $ATTEMPT -lt $MAX_RETRIES ]; then
                echo ""
                echo "* *"
                echo "*"
                echo "*  ⚠️  ENA DOWNLOAD FAILED - LIKELY SERVER ISSUE"
                echo "*  Retrying in $RETRY_DELAY seconds... (attempt $ATTEMPT/$MAX_RETRIES)"
                echo "*"
                echo "* *"
                echo ""
                sleep $RETRY_DELAY
                # Clean up partial files
                rm -rf {wildcards.acc}/ {wildcards.acc}.fastq* {wildcards.acc}_*.fastq* 2>/dev/null || true
            fi
        done
        
        # All retries exhausted
        echo ""
        echo "* *"
        echo "*"
        echo "*  ❌ ENA DOWNLOAD FAILED AFTER $MAX_RETRIES ATTEMPTS"
        echo "*"
        echo "*  This is likely due to ENA/EBI server issues."
        echo "*  Please try again later or check:"
        echo "*  https://www.ebi.ac.uk/ena/browser/view/{wildcards.acc}"
        echo "*"
        echo "* *"
        echo ""
        ls -lh {params.outdir}/
        exit 1
        """


rule _00_download_reads_sra:
    """Download paired-end reads from SRA/ENA"""
    output:
        r1 = "{outdir}/downloaded_data/{species}/reads/{read_type}/{acc}_1.fastq.gz",
        r2 = "{outdir}/downloaded_data/{species}/reads/{read_type}/{acc}_2.fastq.gz"
    params:
        outdir = lambda w: os.path.join(w.outdir, "downloaded_data", w.species, "reads", w.read_type)
    threads: cpu_func("download_data")
    resources:
        mem_mb = mem_func("download_data"),
        runtime = time_func("download_data")
    container: CONTAINERS["gep2_base"]
    shell:
        """
        echo "[GEP2] Downloading paired-end reads: {wildcards.acc}"
        
        mkdir -p {params.outdir}
        cd {params.outdir}
        
        MAX_RETRIES=3
        RETRY_DELAY=60
        
        for ATTEMPT in $(seq 1 $MAX_RETRIES); do
            echo "[GEP2] Download attempt $ATTEMPT of $MAX_RETRIES..."
            
            if enaDataGet.py -f fastq -d . {wildcards.acc}; then
                # Move files from subdirectory if created
                if [ -d "{wildcards.acc}" ]; then
                    mv {wildcards.acc}/* . 2>/dev/null || true
                    rmdir {wildcards.acc} 2>/dev/null || true
                fi
                
                # Compress if needed
                if [ -f "{wildcards.acc}_1.fastq" ]; then
                    pigz -p {threads} "{wildcards.acc}_1.fastq"
                fi
                if [ -f "{wildcards.acc}_2.fastq" ]; then
                    pigz -p {threads} "{wildcards.acc}_2.fastq"
                fi
                
                # Check if we got both files
                if [ -f "{output.r1}" ] && [ -f "{output.r2}" ]; then
                    echo "[GEP2] ✅ Downloaded paired reads: {wildcards.acc}"
                    exit 0
                fi
            fi
            
            # If we get here, download failed
            if [ $ATTEMPT -lt $MAX_RETRIES ]; then
                echo ""
                echo "* *"
                echo "*"
                echo "*  ⚠️  ENA DOWNLOAD FAILED - LIKELY SERVER ISSUE"
                echo "*  Retrying in $RETRY_DELAY seconds... (attempt $ATTEMPT/$MAX_RETRIES)"
                echo "*"
                echo "* *"
                echo ""
                sleep $RETRY_DELAY
                # Clean up partial files
                rm -rf {wildcards.acc}/ {wildcards.acc}_*.fastq* 2>/dev/null || true
            fi
        done
        
        # All retries exhausted
        echo ""
        echo "* *"
        echo "*"
        echo "*  ❌ ENA DOWNLOAD FAILED AFTER $MAX_RETRIES ATTEMPTS"
        echo "*"
        echo "*  This is likely due to ENA/EBI server issues."
        echo "*  Please try again later or check:"
        echo "*  https://www.ebi.ac.uk/ena/browser/view/{wildcards.acc}"
        echo "*"
        echo "* *"
        echo ""
        ls -lh {params.outdir}/
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
import json
with open('{params.manifest}') as f:
    manifest = json.load(f)
for item in manifest:
    if item.get('type') == 'reads' and item.get('method') == 'curl' and item['destination'] == '{output.reads}':
        print(item['source'])
        break
")
        
        if [ -z \"$SOURCE\" ]; then
            echo \"Error: No URL source found in manifest for {output.reads}\"
            exit 1
        fi
        
        mkdir -p $(dirname {output.reads})
        
        echo \"[GEP2] Downloading reads from URL: $SOURCE\"
        curl -L -C - --retry 3 --retry-delay 5 -o {output.reads}.tmp \"$SOURCE\"
        
        if [ ! -s {output.reads}.tmp ]; then
            echo \"Error: Downloaded file is empty\"
            exit 1
        fi
        
        mv {output.reads}.tmp {output.reads}
        echo \"[GEP2] ✅ Downloaded reads: {output.reads}\"
        """

