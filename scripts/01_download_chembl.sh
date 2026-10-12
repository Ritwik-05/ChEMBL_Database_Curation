#!/usr/bin/env bash
# =============================================================================
# Script: 01_download_chembl.sh
# Description: Download ChEMBL 37 SQLite database, schema documentation,
#              UniProt mapping, release notes, and checksums.
#              Supports resume (-c) and verifies checksums after download.
# Usage: bash scripts/01_download_chembl.sh
# =============================================================================

set -euo pipefail

# --- Configuration -----------------------------------------------------------
BASE_URL="https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest"
CHEMBL_VERSION="chembl_37"
RAW_DIR="data/raw"
SCHEMA_DIR="data/schema"
LOG_DIR="logs"
LOGFILE="${LOG_DIR}/download_$(date +%Y%m%d_%H%M%S).log"

# Files to download
FILES=(
    "chembl_37_sqlite.tar.gz"      # Main SQLite DB  (~5.4 GB)
    "chembl_uniprot_mapping.txt"   # ChEMBL -> UniProt mapping (~1.2 MB)
    "chembl_37_release_notes.txt"  # Release notes (~12 KB)
    "schema_documentation.txt"     # Schema documentation (~94 KB)
    "schema_documentation.html"    # Schema documentation HTML (~114 KB)
    "chembl_37_schema.pdf"         # Schema PDF diagram (~800 KB)
    "checksums.txt"                # MD5/SHA checksums
    "README"                       # README
    "LICENSE"                      # License file
)

# --- Setup -------------------------------------------------------------------
mkdir -p "${RAW_DIR}" "${SCHEMA_DIR}" "${LOG_DIR}"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOGFILE}"
}

log "========================================================="
log " ChEMBL 37 Download Script"
log "========================================================="
log "Base URL : ${BASE_URL}"
log "Raw dir  : ${RAW_DIR}"
log "Log file : ${LOGFILE}"
log ""

# --- Check dependencies ------------------------------------------------------
for cmd in wget md5sum tar; do
    if ! command -v "$cmd" &>/dev/null; then
        log "ERROR: Required command '${cmd}' not found. Install it and retry."
        exit 1
    fi
done

# --- Download files ----------------------------------------------------------
download_file() {
    local filename="$1"
    local url="${BASE_URL}/${filename}"
    local outdir="${RAW_DIR}"

    # Schema files go to schema dir
    if [[ "${filename}" == *.pdf || "${filename}" == schema_documentation* ]]; then
        outdir="${SCHEMA_DIR}"
    fi

    local outpath="${outdir}/${filename}"

    log "Downloading: ${filename}"
    wget \
        --continue \
        --tries=5 \
        --waitretry=10 \
        --progress=dot:giga \
        --output-document="${outpath}" \
        --append-output="${LOGFILE}" \
        "${url}" \
        && log "  -> Done: ${outpath}" \
        || { log "  -> FAILED: ${filename}"; return 1; }
}

for file in "${FILES[@]}"; do
    if [[ -f "${RAW_DIR}/${file}" || -f "${SCHEMA_DIR}/${file}" ]]; then
        log "Skipping (already exists): ${file}"
    else
        download_file "${file}"
    fi
done

# --- Verify checksums --------------------------------------------------------
log ""
log "Verifying checksums..."
CHECKSUM_FILE="${RAW_DIR}/checksums.txt"

if [[ -f "${CHECKSUM_FILE}" ]]; then
    cd "${RAW_DIR}"
    # checksums.txt contains md5 hashes for the compressed files
    # Filter only the .tar.gz file we downloaded
    if grep -q "${CHEMBL_VERSION}_sqlite.tar.gz" "${CHECKSUM_FILE}"; then
        grep "${CHEMBL_VERSION}_sqlite.tar.gz" "${CHECKSUM_FILE}" | md5sum -c - 2>&1 | tee -a "../${LOGFILE}"
    fi
    cd - > /dev/null
    log "Checksum verification complete."
else
    log "WARNING: checksums.txt not found. Skipping verification."
fi

# --- Extract SQLite DB -------------------------------------------------------
log ""
TARBALL="${RAW_DIR}/${CHEMBL_VERSION}_sqlite.tar.gz"
SQLITE_DB="${RAW_DIR}/${CHEMBL_VERSION}.db"

if [[ -f "${SQLITE_DB}" ]]; then
    log "SQLite database already extracted: ${SQLITE_DB}"
else
    if [[ -f "${TARBALL}" ]]; then
        log "Extracting SQLite tarball (this may take several minutes)..."
        tar -xzf "${TARBALL}" -C "${RAW_DIR}" 2>&1 | tee -a "${LOGFILE}"
        # ChEMBL extracts as chembl_37/chembl_37_sqlite/chembl_37.db
        # Find and move to flat location
        EXTRACTED_DB=$(find "${RAW_DIR}" -name "*.db" | head -1)
        if [[ -n "${EXTRACTED_DB}" && "${EXTRACTED_DB}" != "${SQLITE_DB}" ]]; then
            mv "${EXTRACTED_DB}" "${SQLITE_DB}"
            # Cleanup empty dirs
            find "${RAW_DIR}" -type d -empty -delete 2>/dev/null || true
        fi
        log "Extracted to: ${SQLITE_DB}"
    else
        log "ERROR: Tarball not found at ${TARBALL}"
        exit 1
    fi
fi

log ""
log "========================================================="
log " Download complete! Summary:"
log "  SQLite DB  : ${SQLITE_DB}"
log "  Schema     : ${SCHEMA_DIR}/"
log "  Mapping    : ${RAW_DIR}/chembl_uniprot_mapping.txt"
log "========================================================="
log ""
log "Next step: python scripts/02_explore_schema.py"
