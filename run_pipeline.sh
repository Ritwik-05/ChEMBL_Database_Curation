#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — ChEMBL 37 Full Curation Pipeline
# =============================================================================
# Runs the complete pipeline:
#   Step 1: Download ChEMBL 37 SQLite + schema files
#   Step 2: Explore schema and generate schema reports
#   Step 3: Curate all compound bioactivity data
#   Step 4: Curate compound-target interactions
#   Step 5: Generate summary report
#   Step 6: Extract well-characterized protein directories and endpoint tables
#   Step 7: Extract IC50 assay experimental conditions and detection methods
#
# Usage:
#   bash run_pipeline.sh                # Run all steps
#   bash run_pipeline.sh --from-step 3  # Resume from step 3
#   bash run_pipeline.sh --skip-download # Skip download (DB already present)
# =============================================================================

set -euo pipefail

# --- Parse arguments ---------------------------------------------------------
FROM_STEP=1
SKIP_DOWNLOAD=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --from-step)   FROM_STEP="$2"; shift 2 ;;
        --skip-download) SKIP_DOWNLOAD=true; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# --- Setup -------------------------------------------------------------------
LOG_DIR="logs"
MAIN_LOG="${LOG_DIR}/pipeline_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "${LOG_DIR}"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${MAIN_LOG}"
}

run_step() {
    local step_num="$1"
    local step_name="$2"
    local cmd="$3"

    if [[ "${step_num}" -lt "${FROM_STEP}" ]]; then
        log "Skipping Step ${step_num}: ${step_name} (--from-step=${FROM_STEP})"
        return 0
    fi

    log ""
    log "=================================================="
    log " STEP ${step_num}: ${step_name}"
    log "=================================================="
    log "Command: ${cmd}"
    log ""

    eval "${cmd}" 2>&1 | tee -a "${MAIN_LOG}"
    log "Step ${step_num} complete."
}

# --- Main pipeline -----------------------------------------------------------
log "=========================================================="
log " ChEMBL 37 Database Curation Pipeline"
log " Started: $(date)"
log "=========================================================="
log "Config:"
log "  FROM_STEP     : ${FROM_STEP}"
log "  SKIP_DOWNLOAD : ${SKIP_DOWNLOAD}"
log "  Log file      : ${MAIN_LOG}"
log ""

# Check Python
if ! command -v python3 &>/dev/null; then
    log "ERROR: python3 not found. Install Python 3.9+."
    exit 1
fi

PYTHON_VERSION=$(python3 --version)
log "Python: ${PYTHON_VERSION}"

# Check if virtual environment exists
if [[ -d ".venv" ]]; then
    log "Activating virtual environment: .venv"
    source .venv/bin/activate
else
    log "WARNING: No .venv found. Running with system Python."
    log "Recommended: python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt"
fi

# --- Run steps ---------------------------------------------------------------

# Step 1: Download
if [[ "${SKIP_DOWNLOAD}" == "false" ]]; then
    run_step 1 "Download ChEMBL 37" "bash scripts/01_download_chembl.sh"
else
    log "Skipping Step 1: Download (--skip-download set)"
fi

# Step 2: Schema exploration
run_step 2 "Explore Database Schema" "python3 scripts/02_explore_schema.py"

# Step 3: Bioactivity curation
run_step 3 "Curate Compound Bioactivity" "python3 scripts/03_curate_bioactivity.py"

# Step 4: Target interaction curation
run_step 4 "Curate Target Interactions" "python3 scripts/04_curate_target_interactions.py"

# Step 5: Report generation
run_step 5 "Generate Summary Report" "python3 scripts/05_generate_report.py"

# Step 6: Well-characterized targets extraction
run_step 6 "Extract Well-Characterized Targets" "python3 scripts/06_extract_well_characterized_targets.py"

# Step 7: Assay conditions extraction
run_step 7 "Extract IC50 Assay Conditions" "python3 scripts/07_extract_assay_conditions.py"

# --- Done --------------------------------------------------------------------
log ""
log "=========================================================="
log " Pipeline Complete!"
log " Finished: $(date)"
log ""
log " Outputs:"
log "   data/curated/chembl_37_bioactivity_broad.tsv.gz"
log "   data/curated/chembl_37_compound_target_interactions_tiered.tsv.gz"
log "   data/curated/chembl_37_drug_dictionary.tsv"
log "   data/curated/chembl_37_target_components.tsv"
log "   data/curated/chembl_37_protein_classification.tsv"
log "   data/curated/assay_conditions/ic50_assay_conditions.tsv.gz"
log "   data/curated/assay_conditions/assay_protocols_conditions_index.tsv"
log "   data/targets/ (902 well-characterized target directories)"
log "   data/targets/targets_index.tsv"
log "   reports/well_characterized_proteins_drug_stats.tsv"
log "   reports/well_characterized_proteins_summary.md"
log "   reports/assay_conditions_stats.tsv"
log "   reports/assay_conditions_summary.md"
log "   reports/chembl_37_curation_summary.md"
log "   reports/run_manifest.json"
log "   METHODOLOGY.md"
log ""
log " Log file: ${MAIN_LOG}"
log "=========================================================="
