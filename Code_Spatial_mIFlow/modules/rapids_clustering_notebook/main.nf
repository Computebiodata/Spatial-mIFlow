process RAPIDS_CLUSTERING_NOTEBOOK {

    tag "${sample_dir_name}"

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/03_RAPIDS_clustering",
        mode: 'copy',
        overwrite: true,
        failOnError: true

    input:
    tuple val(sample_index),
          val(sample_id),
          val(sample_dir_name),
          path(input_adata)

    output:
    path "${sample_dir_name}", emit: rapids_outputs

    tuple val(sample_index),
          val(sample_id),
          val(sample_dir_name),
          path("${sample_dir_name}/clustering/rapids_clustered.h5ad"),
          emit: rapids_adata

    script:
    def RENDER_HTML = params.notebooks.rapids.render_quarto_html ? "true" : "false"
    def RENDER_PDF = params.notebooks.rapids.render_quarto_pdf ? "true" : "false"
    def NOTEBOOK = "${moduleDir}/rapids_clustering_template.ipynb"

    """
    set -euo pipefail

    mkdir -p "${sample_dir_name}/report"

    EXECUTED_NOTEBOOK="${sample_dir_name}/report/${sample_dir_name}_rapids_clustering_executed.ipynb"

    echo "============================================================"
    echo "RAPIDS CLUSTERING NOTEBOOK"
    echo "============================================================"
    echo "Sample index: ${sample_index}"
    echo "Sample ID:    ${sample_id}"
    echo "Output dir:   ${sample_dir_name}"
    echo "Input AnnData:"
    readlink -f "${input_adata}"
    echo "Notebook:"
    echo "${NOTEBOOK}"
    echo "Executed notebook:"
    echo "\${EXECUTED_NOTEBOOK}"
    echo "============================================================"

    echo "Python executable:"
    command -v python
    python --version

    echo "Checking notebook tools:"
    python - <<'PY_CHECK'
import nbformat
import papermill
import ipykernel

print("nbformat:", nbformat.__version__)
print("papermill:", papermill.__version__)
print("ipykernel:", ipykernel.__version__)
PY_CHECK

    echo "Ensuring python3 kernel exists in the conda environment..."
    # Installing into --sys-prefix means every parallel task sharing this
    # same conda env (one per sample) targets the identical destination
    # directory. ipykernel's installer rmtree's-then-recreates it, so two
    # samples' tasks landing here at once race: one's rmtree can hit an
    # entry the other already removed (verified: FileNotFoundError from
    # shutil.rmtree). flock serializes the install; the existence check
    # then makes every task after the first a no-op instead of a reinstall.
    KERNEL_SPEC_DIR="\$(python -c 'import sys; print(sys.prefix)')/share/jupyter/kernels/python3"
    mkdir -p "\$(dirname "\${KERNEL_SPEC_DIR}")"
    (
        flock -x 200
        if [[ ! -d "\${KERNEL_SPEC_DIR}" ]]; then
            python -m ipykernel install --sys-prefix --name python3 --display-name "Python 3"
        fi
    ) 200>"\$(dirname "\${KERNEL_SPEC_DIR}")/.install.lock"

    echo "Running notebook with Papermill..."

    python -m papermill "${NOTEBOOK}" "\${EXECUTED_NOTEBOOK}" \\
        -k python3 \\
        --log-output \\
        -p SAMPLE_INDEX "${sample_index}" \\
        -p SAMPLE_ID "${sample_id}" \\
        -p SAMPLE_DIR_NAME "${sample_dir_name}" \\
        -p INPUT_PATH "\$(readlink -f "${input_adata}")" \\
        -p OUTPUT_DIR "\$(readlink -f "${sample_dir_name}")" \\
        -p DAPI_QC_PERCENTILE ${params.notebooks.rapids.dapi_qc_percentile}

    echo "Validating expected notebook outputs..."

    python - <<'PY_VALIDATE'
from pathlib import Path
import sys

sample_dir = Path("${sample_dir_name}")

expected_files = [
    (sample_dir / "dapi_qc" / "dapi_qc_summary.png", "DAPI QC summary PNG"),
    (sample_dir / "dapi_qc" / "dapi_histogram.png", "DAPI histogram PNG"),
    (sample_dir / "clustering" / "rapids_clustered.h5ad", "final clustered AnnData"),
    (sample_dir / "report" / "session_summary.png", "session summary PNG"),
    (sample_dir / "report" / "${sample_dir_name}_rapids_clustering_executed.ipynb", "executed notebook"),
]

missing = []

for path, label in expected_files:
    if (not path.exists()) or path.stat().st_size == 0:
        missing.append((path, label))

spatial_candidates = [
    sample_dir / "dapi_qc" / "dapi_spatial_qc.png",
    sample_dir / "dapi_qc" / "dapi_spatial_kept_cells.png",
]

spatial_found = [
    path for path in spatial_candidates
    if path.exists() and path.stat().st_size > 0
]

if len(spatial_found) == 0:
    missing.append((
        sample_dir / "dapi_qc" / "dapi_spatial_qc.png",
        "DAPI spatial QC PNG",
    ))

if missing:
    print("ERROR: missing expected RAPIDS notebook output(s):")
    for path, label in missing:
        print(f"  - {label}: {path}")

    print()
    print(f"Existing files under {sample_dir}:")
    for path in sorted(sample_dir.rglob("*")):
        if path.is_file():
            print(f"  {path}")

    sys.exit(1)

print("Validated RAPIDS notebook outputs:")
for path, label in expected_files:
    print(f"  - {label}: {path}")

print(f"  - DAPI spatial QC PNG: {spatial_found[0]}")
PY_VALIDATE

    RENDER_HTML="${RENDER_HTML}"
    RENDER_PDF="${RENDER_PDF}"

    if command -v quarto >/dev/null 2>&1; then
        echo "Quarto executable:"
        command -v quarto
        quarto --version

        if [[ "\${RENDER_HTML}" == "true" ]]; then
            quarto render "\${EXECUTED_NOTEBOOK}" \\
                --to html \\
                --output "${sample_dir_name}_rapids_clustering_report.html" \\
                --output-dir "\$(readlink -f "${sample_dir_name}/report")" || \\
                echo "WARNING: Quarto HTML rendering failed, but executed notebook was saved."
        fi

        if [[ "\${RENDER_PDF}" == "true" ]]; then
            quarto render "\${EXECUTED_NOTEBOOK}" \\
                --to pdf \\
                --output "${sample_dir_name}_rapids_clustering_report.pdf" \\
                --output-dir "\$(readlink -f "${sample_dir_name}/report")" || \\
                echo "WARNING: Quarto PDF rendering failed, but executed notebook was saved."
        fi
    else
        echo "WARNING: quarto command not found. Skipping optional report rendering."
    fi

    echo "RAPIDS notebook step completed."
    """
}
