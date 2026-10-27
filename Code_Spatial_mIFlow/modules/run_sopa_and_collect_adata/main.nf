process RUN_SOPA_AND_COLLECT_ADATA {

    tag "sopa"

    cache false

    input:
    path sopa_samplesheet
    path resolved_channels
    path sopa_params

    output:
    path "sopa_adata_manifest.csv", emit: adata_manifest

    script:
    def RESUME_OPT = params.sopa_adapter.resume ? "-resume" : ""
    def NXF_VERSION = params.nxf_ver ?: ""
    def SOPA_SKIP = params.sopa_adapter.skip_if_complete ? "true" : "false"
    def PYTHON_BIN = params.python_bin ?: "python3"
    // outdir may be relative (resolved against launchDir); used as a raw bash
    // path below, and passed to the *inner* SOPA pipeline as its own
    // --outdir, so absolutize here rather than leaving it to bash.
    def ABS_OUTDIR = (params.outdir as String).startsWith('/') ? params.outdir : "${launchDir}/${params.outdir}"

    """
    set -euo pipefail

    if [[ -n "${NXF_VERSION}" ]]; then
        export NXF_VER="${NXF_VERSION}"
    fi

    PROJECT_DIR="${projectDir}"
    SOPA_SETUP="\${PROJECT_DIR}/programs_hierarchy/sopa_setup"
    SOPA_PIPELINE="${params.sopa_pipeline}"
    SOPA_CONFIG="${params.sopa_config}"
    SOPA_WORK="${params.sopa_work}"
    SOPA_OUTPUT_DIR="${ABS_OUTDIR}/02_sopa"
    SOPA_LOG_DIR="${params.sopa_logs_dir}"
    SOPA_RUNTIME_INPUTS_DIR="${ABS_OUTDIR}/00_runtime/sopa_inputs"
    SOPA_SKIP_IF_COMPLETE="${SOPA_SKIP}"
    PYTHON_BIN="${PYTHON_BIN}"

    PIPELINE_TMP="${params.tmp_dir}"
    MATPLOTLIB_CACHE="${params.matplotlib_cache_dir}"
    CELLPOSE_MODEL_CACHE="${params.cellpose_model_cache}"

    TASK_WORKDIR="\$(pwd -P)"
    MANIFEST_PATH="\${TASK_WORKDIR}/sopa_adata_manifest.csv"

    SOPA_SAMPLESHEET_FOR_RUN="\${SOPA_RUNTIME_INPUTS_DIR}/sopa_samplesheet.csv"
    SOPA_RESOLVED_CHANNELS_FOR_RUN="\${SOPA_RUNTIME_INPUTS_DIR}/resolved_sopa_channels.json"
    SOPA_PARAMS_FOR_RUN="\${SOPA_RUNTIME_INPUTS_DIR}/sopa_params.generated.yml"

    mkdir -p "\${SOPA_WORK}"
    mkdir -p "\${SOPA_OUTPUT_DIR}"
    mkdir -p "\${SOPA_LOG_DIR}"
    mkdir -p "\${PIPELINE_TMP}"
    mkdir -p "\${MATPLOTLIB_CACHE}"
    mkdir -p "\${CELLPOSE_MODEL_CACHE}"

    export TMPDIR="\${PIPELINE_TMP}"
    export MPLCONFIGDIR="\${MATPLOTLIB_CACHE}"
    export CELLPOSE_LOCAL_MODELS_PATH="\${CELLPOSE_MODEL_CACHE}"

    if ! command -v "\${PYTHON_BIN}" >/dev/null 2>&1; then
        echo "ERROR: Python executable not found: \${PYTHON_BIN}"
        echo "Set it with --python_bin /path/to/python3"
        exit 127
    fi

    for f in "\${SOPA_SAMPLESHEET_FOR_RUN}" "\${SOPA_RESOLVED_CHANNELS_FOR_RUN}" "\${SOPA_PARAMS_FOR_RUN}"; do
        if [[ ! -s "\${f}" ]]; then
            echo "ERROR: expected SOPA runtime file missing or empty:"
            echo "  \${f}"
            echo
            echo "Available runtime files:"
            ls -lah "\${SOPA_RUNTIME_INPUTS_DIR}" || true
            exit 1
        fi
    done

    cat > "\${TASK_WORKDIR}/write_sopa_manifest.py" <<'PY_WRITE_SOPA_MANIFEST'
import csv
import os
import re
import sys
from pathlib import Path

sopa_samplesheet = Path(os.environ["SOPA_SAMPLESHEET_FOR_RUN"])
sopa_output_dir = Path(os.environ["SOPA_OUTPUT_DIR"])
manifest_path = Path(os.environ["MANIFEST_PATH"])
allow_missing = os.environ.get("ALLOW_MISSING", "false").lower() == "true"

def make_safe_name(x):
    x = str(x).strip()
    x = re.sub(r"[^A-Za-z0-9_.-]+", "_", x)
    x = re.sub(r"_+", "_", x).strip("_")
    return x if x else "sample"

with open(sopa_samplesheet, newline="") as f:
    rows = list(csv.DictReader(f))

if not rows:
    raise SystemExit("SOPA samplesheet is empty.")

rows_out = []
missing = []

for i, row in enumerate(rows, start=1):
    sample_id = str(row["sample"]).strip()

    if not sample_id:
        raise SystemExit(f"Empty sample name at row {i}")

    sample_safe = make_safe_name(sample_id)
    sample_dir_name = f"{i:03d}_{sample_safe}"

    expected_adata = sopa_output_dir / f"{sample_id}.explorer" / "adata.h5ad"

    resolved_adata = None

    if expected_adata.exists():
        resolved_adata = expected_adata
    else:
        candidates = sorted(
            sopa_output_dir.glob(f"*{sample_safe}*.explorer/adata.h5ad")
        )

        if len(candidates) == 1:
            resolved_adata = candidates[0]

    if resolved_adata is None:
        missing.append({
            "sample_id": sample_id,
            "expected": str(expected_adata),
        })
        continue

    explorer_dir = resolved_adata.parent

    rows_out.append({
        "sample_index": i,
        "sample_id": sample_id,
        "sample_dir_name": sample_dir_name,
        "adata_path": str(resolved_adata.resolve()),
        "sopa_explorer_dir": str(explorer_dir.resolve()),
    })

if missing:
    print("Missing SOPA AnnData outputs:")
    for item in missing:
        print(f"  sample_id: {item['sample_id']}")
        print(f"  expected:  {item['expected']}")

    if allow_missing:
        sys.exit(42)

    raise SystemExit("Missing required SOPA AnnData output(s).")

with open(manifest_path, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "sample_index",
            "sample_id",
            "sample_dir_name",
            "adata_path",
            "sopa_explorer_dir",
        ],
    )
    writer.writeheader()
    writer.writerows(rows_out)

print("SOPA AnnData manifest written:")
print(manifest_path)
print(manifest_path.read_text())
PY_WRITE_SOPA_MANIFEST

    echo "============================================================"
    echo "GLOBAL SOPA ADAPTER"
    echo "============================================================"
    echo "SOPA setup:           \${SOPA_SETUP}"
    echo "SOPA pipeline:        \${SOPA_PIPELINE}"
    echo "SOPA config:          \${SOPA_CONFIG}"
    echo "SOPA params:          \${SOPA_PARAMS_FOR_RUN}"
    echo "SOPA samplesheet:     \${SOPA_SAMPLESHEET_FOR_RUN}"
    echo "SOPA output dir:      \${SOPA_OUTPUT_DIR}"
    echo "SOPA work dir:        \${SOPA_WORK}"
    echo "SOPA log dir:         \${SOPA_LOG_DIR}"
    echo "SOPA skip complete:   \${SOPA_SKIP_IF_COMPLETE}"
    echo "TMPDIR:               \${TMPDIR}"
    echo "MPLCONFIGDIR:         \${MPLCONFIGDIR}"
    echo "CELLPOSE cache:       \${CELLPOSE_LOCAL_MODELS_PATH}"
    echo "Manifest path:        \${MANIFEST_PATH}"
    echo "Python executable:    \${PYTHON_BIN}"
    "\${PYTHON_BIN}" --version
    echo "============================================================"

    echo "SOPA runtime samplesheet:"
    cat "\${SOPA_SAMPLESHEET_FOR_RUN}"

    echo
    echo "SOPA runtime generated params:"
    cat "\${SOPA_PARAMS_FOR_RUN}"

    if [[ "\${SOPA_SKIP_IF_COMPLETE}" == "true" ]]; then
        echo
        echo "Checking for existing complete SOPA outputs..."

        set +e
        ALLOW_MISSING=true \
        SOPA_SAMPLESHEET_FOR_RUN="\${SOPA_SAMPLESHEET_FOR_RUN}" \
        SOPA_OUTPUT_DIR="\${SOPA_OUTPUT_DIR}" \
        MANIFEST_PATH="\${MANIFEST_PATH}" \
        "\${PYTHON_BIN}" "\${TASK_WORKDIR}/write_sopa_manifest.py"
        CHECK_STATUS=\$?
        set -e

        if [[ "\${CHECK_STATUS}" -eq 0 ]]; then
            echo "Existing SOPA outputs detected and validated. Skipping inner SOPA run."
            exit 0
        elif [[ "\${CHECK_STATUS}" -eq 42 ]]; then
            echo "Some SOPA outputs are missing. Inner SOPA will run."
        else
            echo "SOPA output check failed with status \${CHECK_STATUS}."
            exit "\${CHECK_STATUS}"
        fi
    fi

    cd "\${SOPA_SETUP}"

    RUN_NAME="global_sopa_\${SLURM_JOB_ID:-manual}_\$(date +%Y%m%d_%H%M%S)"

    echo "============================================================"
    echo "RUNNING INNER NF-CORE/SOPA"
    echo "============================================================"
    echo "Run name: \${RUN_NAME}"
    echo "Resume:   ${RESUME_OPT}"
    echo "============================================================"

    echo "Nextflow version:"
    nextflow -version

    NF_RESUME_OPT="${RESUME_OPT}"

    NEXTFLOW_CMD=(
        nextflow
        -c "\${SOPA_CONFIG}"
        run "\${SOPA_PIPELINE}"
        -name "\${RUN_NAME}"
        -profile "${params.sopa_adapter.profile}"
        -params-file "\${SOPA_PARAMS_FOR_RUN}"
        --input "\${SOPA_SAMPLESHEET_FOR_RUN}"
        --outdir "\${SOPA_OUTPUT_DIR}"
        --cpu_queue "${params.cpu_queue}"
        --gpu_queue "${params.gpu_queue}"
        -work-dir "\${SOPA_WORK}"
    )

    if [[ -n "\${NF_RESUME_OPT}" ]]; then
        NEXTFLOW_CMD+=("\${NF_RESUME_OPT}")
    fi

    NEXTFLOW_CMD+=(
        -with-report "\${SOPA_LOG_DIR}/\${RUN_NAME}_report.html"
        -with-trace "\${SOPA_LOG_DIR}/\${RUN_NAME}_trace.txt"
        -with-timeline "\${SOPA_LOG_DIR}/\${RUN_NAME}_timeline.html"
    )

    echo "Executing inner SOPA command:"
    printf ' %q' "\${NEXTFLOW_CMD[@]}"
    echo

    "\${NEXTFLOW_CMD[@]}"

    echo "Inner SOPA completed. Validating and collecting AnnData outputs..."

    ALLOW_MISSING=false \
    SOPA_SAMPLESHEET_FOR_RUN="\${SOPA_SAMPLESHEET_FOR_RUN}" \
    SOPA_OUTPUT_DIR="\${SOPA_OUTPUT_DIR}" \
    MANIFEST_PATH="\${MANIFEST_PATH}" \
    "\${PYTHON_BIN}" "\${TASK_WORKDIR}/write_sopa_manifest.py"

    echo "SOPA adapter completed successfully."
    """
}
