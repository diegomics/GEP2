#!/bin/bash
# GEP2 - Download with stall-protected timeout
# Usage: gep2_download_with_timeout <stall_timeout> <command> [args...]
# Uses 'timeout' to enforce a hard time limit on the download.

gep2_download_with_timeout() {
    local STALL_TIMEOUT=${1:-3600}
    shift
    
    echo "[GEP2] Running with ${STALL_TIMEOUT}s timeout: $@"
    timeout --signal=TERM --kill-after=30 "$STALL_TIMEOUT" "$@"
    return $?
}