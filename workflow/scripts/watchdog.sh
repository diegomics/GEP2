#!/bin/bash
# GEP2 - Download helpers

gep2_download_with_timeout() {
    local STALL_TIMEOUT=${1:-3600}
    shift

    echo "[GEP2] Running with ${STALL_TIMEOUT}s timeout: $@"
    timeout --signal=TERM --kill-after=30 "$STALL_TIMEOUT" "$@"
    return $?
}


# ---------------------------------------------------------------------
# HTTP fetch (to stdout): curl preferred, wget fallback, python urllib
# last resort. Used for short API responses where we need the body
# inline. Retries built in. Diagnostics to stderr.
# ---------------------------------------------------------------------
gep2_http_fetch_stdout() {
    # gep2_http_fetch_stdout <url>
    local URL=$1

    if command -v curl >/dev/null 2>&1; then
        curl -sS --max-time 60 --retry 3 --retry-delay 10 "$URL"
    elif command -v wget >/dev/null 2>&1; then
        wget -qO- --tries=3 --waitretry=10 --timeout=60 "$URL"
    elif command -v python3 >/dev/null 2>&1 || command -v python >/dev/null 2>&1; then
        local PY
        PY=$(command -v python3 || command -v python)
        "$PY" -c '
import sys, urllib.request
req = urllib.request.Request(sys.argv[1], headers={"User-Agent": "GEP2/urllib"})
sys.stdout.write(urllib.request.urlopen(req, timeout=60).read().decode())
' "$URL"
    else
        echo "[GEP2] No HTTP client available in container" >&2
        return 127
    fi
}

# ---------------------------------------------------------------------
# HTTP fetch (to file): same backend preference as above, designed for
# large downloads. Retries and generous timeouts built in.
# ---------------------------------------------------------------------
gep2_http_fetch() {
    # gep2_http_fetch <url> <output_file>
    # Downloads with retries, returns tool exit status.
    local URL=$1
    local OUT=$2

    if command -v curl >/dev/null 2>&1; then
        curl -sS -L --fail \
            --connect-timeout 30 \
            --max-time 14400 \
            --retry 5 \
            --retry-delay 60 \
            --retry-max-time 10800 \
            -o "$OUT" "$URL"
    elif command -v wget >/dev/null 2>&1; then
        wget -q --tries=5 --waitretry=60 --timeout=14400 -O "$OUT" "$URL"
    elif command -v python3 >/dev/null 2>&1 || command -v python >/dev/null 2>&1; then
        local PY
        PY=$(command -v python3 || command -v python)
        # Stdlib-only: urllib.request with manual retries. No deps.
        "$PY" - "$URL" "$OUT" <<'PYEOF'
import sys, time, urllib.request, urllib.error, shutil
url, out = sys.argv[1], sys.argv[2]
for attempt in range(1, 6):
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "GEP2/urllib"})
        with urllib.request.urlopen(req, timeout=14400) as r, open(out, "wb") as f:
            shutil.copyfileobj(r, f, length=1024 * 1024)
        sys.exit(0)
    except (urllib.error.URLError, OSError) as e:
        sys.stderr.write(f"[GEP2]   python urllib attempt {attempt} failed: {e}\n")
        if attempt < 5:
            time.sleep(60)
sys.exit(1)
PYEOF
    else
        echo "[GEP2] No HTTP client (curl/wget/python) available in container" >&2
        return 127
    fi
}


# ---------------------------------------------------------------------
# Query ENA portal API filereport for a run accession and print the
# HTTPS URLs (one per line) to stdout. Diagnostics go to stderr.
#
# Prefers fastq_ftp; falls back to submitted_ftp when submitted_format
# contains FASTQ.
#
# Exit codes:
#   0 success
#   2 network error reaching ENA portal
#   3 empty response / no data row (accession not found?)
#   4 no FASTQ available (submission may be CRAM/BAM-only)
#   5 URL count doesn't match requested mode
# ---------------------------------------------------------------------
gep2_ena_get_urls() {
    local ACC=$1
    local MODE=${2:-paired}   # paired | single
    local API_URL="https://www.ebi.ac.uk/ena/portal/api/filereport"
    local FIELDS="fastq_ftp,submitted_ftp,submitted_format"
    local META ROW URLS SOURCE COUNT FETCH_EXIT

    echo "[GEP2] Querying ENA portal API for ${ACC} (mode=${MODE})..." >&2

    META=$(gep2_http_fetch_stdout \
        "${API_URL}?accession=${ACC}&result=read_run&fields=${FIELDS}&format=tsv")
    FETCH_EXIT=$?
    if [ "$FETCH_EXIT" -ne 0 ]; then
        echo "[GEP2] ENA portal API query failed (exit ${FETCH_EXIT})" >&2
        return 2
    fi

    ROW=$(echo "$META" | awk 'NR==2')
    if [ -z "$ROW" ]; then
        echo "[GEP2] ENA portal returned no data row for ${ACC}" >&2
        echo "[GEP2] Response was:" >&2
        echo "$META" >&2
        return 3
    fi

    local FASTQ_FTP SUBMITTED_FTP SUBMITTED_FORMAT
    # ENA portal always prepends run_accession as column 1, regardless
    # of what fields= requests. Our requested fields start at column 2.
    FASTQ_FTP=$(echo "$ROW"        | awk -F'\t' '{print $2}')
    SUBMITTED_FTP=$(echo "$ROW"    | awk -F'\t' '{print $3}')
    SUBMITTED_FORMAT=$(echo "$ROW" | awk -F'\t' '{print $4}')

    if [ -n "$FASTQ_FTP" ]; then
        URLS="$FASTQ_FTP"
        SOURCE="fastq_ftp"
    elif [ -n "$SUBMITTED_FTP" ] && echo "$SUBMITTED_FORMAT" | grep -qi "FASTQ"; then
        URLS="$SUBMITTED_FTP"
        SOURCE="submitted_ftp (FASTQ)"
    else
        echo "[GEP2] No FASTQ available via ENA portal for ${ACC}" >&2
        echo "[GEP2]   fastq_ftp='${FASTQ_FTP}'" >&2
        echo "[GEP2]   submitted_ftp='${SUBMITTED_FTP}'" >&2
        echo "[GEP2]   submitted_format='${SUBMITTED_FORMAT}'" >&2
        return 4
    fi

    COUNT=$(echo "$URLS" | tr ';' '\n' | grep -c .)
    if [ "$MODE" = "paired" ] && [ "$COUNT" -ne 2 ]; then
        echo "[GEP2] Expected 2 URLs for paired-end ${ACC}, got ${COUNT}" >&2
        return 5
    fi
    if [ "$MODE" = "single" ] && [ "$COUNT" -lt 1 ]; then
        echo "[GEP2] No URLs for ${ACC}" >&2
        return 5
    fi

    echo "[GEP2] Found ${COUNT} URL(s) for ${ACC} via ${SOURCE}" >&2
    # Emit URLs one per line, prepending https:// (API returns them
    # as bare hostnames, e.g. ftp.sra.ebi.ac.uk/vol1/...)
    echo "$URLS" | tr ';' '\n' | awk 'NF > 0 {print "https://" $0}'
    return 0
}


# ---------------------------------------------------------------------
# Download paired-end reads via portal API + HTTPS.
# CWD must be the destination directory.
# Produces <acc>_1.fastq.gz and <acc>_2.fastq.gz on success.
# ---------------------------------------------------------------------
gep2_ena_download_paired() {
    local ACC=$1
    local URLS URL1 URL2 API_EXIT

    URLS=$(gep2_ena_get_urls "$ACC" paired)
    API_EXIT=$?
    if [ "$API_EXIT" -ne 0 ]; then
        return "$API_EXIT"
    fi

    URL1=$(echo "$URLS" | sed -n '1p')
    URL2=$(echo "$URLS" | sed -n '2p')

    gep2_http_fetch "$URL1" "${ACC}_1.fastq.gz" || return 10
    gep2_http_fetch "$URL2" "${ACC}_2.fastq.gz" || return 11
    return 0
}


# ---------------------------------------------------------------------
# Download single-end/long reads via portal API + HTTPS.
# CWD must be the destination directory.
# Produces <acc>.fastq.gz on success.
# ---------------------------------------------------------------------
gep2_ena_download_single() {
    local ACC=$1
    local URLS URL1 API_EXIT

    URLS=$(gep2_ena_get_urls "$ACC" single)
    API_EXIT=$?
    if [ "$API_EXIT" -ne 0 ]; then
        return "$API_EXIT"
    fi

    URL1=$(echo "$URLS" | sed -n '1p')
    gep2_http_fetch "$URL1" "${ACC}.fastq.gz" || return 10
    return 0
}