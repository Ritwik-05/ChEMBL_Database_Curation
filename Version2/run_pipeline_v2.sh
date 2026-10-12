#!/usr/bin/env bash
# =============================================================================
# run_pipeline_v2.sh — ChEMBL 37 Audited Pipeline Runner (Version 2)
# =============================================================================
# Executes the corrected curation pipeline in Version2/ with separated outputs:
#   Step 2: Schema exploration & condition table profiling
#   Step 3: Bioactivity curation with explicit rejection ledger
#   Step 4: Orthogonal compound-target association stratification
#   Step 6: Target directory extraction under strict single-protein predicate
#   Step 7: IC50 condition extraction (time, ATP, format, detection)
#   Step 5: Master reporting & SHA-256 reproducibility manifest
# =============================================================================

set -euo pipefail

LOG_DIR="Version2/logs"
MAIN_LOG="${LOG_DIR}/pipeline_v2_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "${LOG_DIR}"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${MAIN_LOG}"
}

run_step() {
    local step_num="$1"
    local step_name="$2"
    local cmd="$3"

    log ""
    log "=================================================="
    log " STEP ${step_num}: ${step_name}"
    log "=================================================="
    log "Command: ${cmd}"
    log ""

    eval "${cmd}" 2>&1 | tee -a "${MAIN_LOG}"
    log "Step ${step_num} complete."
}

log "=========================================================="
log " ChEMBL 37 Database Curation Pipeline — Version 2"
log " Started: $(date)"
log " Log: ${MAIN_LOG}"
log "=========================================================="

# Check Python environment
if [[ -d ".venv" ]]; then
    source .venv/bin/activate
fi

# Step 2: Schema exploration
run_step 2 "Schema Exploration & Profiling" "python3 Version2/scripts/02_explore_schema.py"

# Step 3: Bioactivity curation & rejection ledger
run_step 3 "Bioactivity Curation & Rejection Accounting" "python3 Version2/scripts/03_curate_bioactivity.py"

# Step 4: Compound-target association stratification
run_step 4 "Compound-Target Association Stratification" "python3 Version2/scripts/04_curate_target_interactions.py"

# Step 6: Well-characterized target extraction (strict predicate)
run_step 6 "Target Partitioning & Endpoint Extraction" "python3 Version2/scripts/06_extract_well_characterized_targets.py"

# Step 7: IC50 assay experimental conditions
run_step 7 "IC50 Experimental Conditions Extraction" "python3 Version2/scripts/07_extract_assay_conditions.py"

# Step 5: Master reporting & manifest (runs last to hash all outputs)
run_step 5 "Master Reporting & Cryptographic Manifest" "python3 Version2/scripts/05_generate_report.py"

log ""
log "=========================================================="
log " Version 2 Pipeline Execution Complete!"
log " Finished: $(date)"
log "=========================================================="
