process CONTINGENCY_CLUSTER_POPULATIONS_NOTEBOOK {

    tag "${sample_dir_name}"

    // Reuses the RAPIDS clustering environment; contingency_setup has no envs/ of its own.
    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/06_CONTINGENCY_clusters_populations",
        mode: 'copy',
        overwrite: true,
        failOnError: true

    input:
    tuple val(sample_index),
          val(sample_id),
          val(sample_dir_name),
          path(input_adata)

    output:
    path "${sample_dir_name}", emit: contingency_outputs

    script:
    def NOTEBOOK = "${moduleDir}/contingency_clusters_populations_template.ipynb"

    """
    set -euo pipefail

    mkdir -p "${sample_dir_name}/report"

    # Temporary local copy.
    # The notebook can patch this copy if needed.
    # The upstream SCIMAP AnnData is not modified.
    cp -L "${input_adata}" "input_for_contingency.h5ad"

    EXECUTED_NOTEBOOK="${sample_dir_name}/report/${sample_dir_name}_contingency_clusters_populations_executed.ipynb"

    echo "============================================================"
    echo "CLUSTER × GATED-POPULATION CONTINGENCY NOTEBOOK"
    echo "============================================================"
    echo "Sample index: ${sample_index}"
    echo "Sample ID:    ${sample_id}"
    echo "Output dir:   ${sample_dir_name}"
    echo "Input SCIMAP AnnData:"
    readlink -f "${input_adata}"
    echo "Temporary work copy:"
    readlink -f "input_for_contingency.h5ad"
    echo "Notebook:"
    echo "${NOTEBOOK}"
    echo "Executed notebook:"
    echo "\${EXECUTED_NOTEBOOK}"
    echo "Top populations: ${params.notebooks.contingency.top_n_populations}"
    echo "Output: one count PDF per Leiden clustering"
    echo "============================================================"

    echo "Python executable:"
    command -v python
    python --version

    echo "Checking contingency notebook environment:"
    python - <<'PY_CHECK'
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib
import seaborn as sns
import h5py
import papermill
import ipykernel
import nbformat

print("anndata:", ad.__version__)
print("pandas:", pd.__version__)
print("numpy:", np.__version__)
print("matplotlib:", matplotlib.__version__)
print("seaborn:", sns.__version__)
print("h5py:", h5py.__version__)
print("papermill:", papermill.__version__)
print("ipykernel:", ipykernel.__version__)
print("nbformat:", nbformat.__version__)
PY_CHECK

    echo "Ensuring python3 kernel exists in the active environment..."
    # See rapids_clustering_notebook/main.nf for why this is flock-guarded:
    # parallel per-sample tasks share this same conda env's --sys-prefix
    # kernel destination, and ipykernel's rmtree-then-recreate install races
    # across concurrent tasks without serialization.
    KERNEL_SPEC_DIR="\$(python -c 'import sys; print(sys.prefix)')/share/jupyter/kernels/python3"
    mkdir -p "\$(dirname "\${KERNEL_SPEC_DIR}")"
    (
        flock -x 200
        if [[ ! -d "\${KERNEL_SPEC_DIR}" ]]; then
            python -m ipykernel install \\
                --sys-prefix \\
                --name python3 \\
                --display-name "Python 3"
        fi
    ) 200>"\$(dirname "\${KERNEL_SPEC_DIR}")/.install.lock"

    echo "Running contingency notebook with Papermill..."

    python -m papermill \\
        "${NOTEBOOK}" \\
        "\${EXECUTED_NOTEBOOK}" \\
        -k python3 \\
        --log-output \\
        -p SAMPLE_INDEX "${sample_index}" \\
        -p SAMPLE_ID "${sample_id}" \\
        -p SAMPLE_DIR_NAME "${sample_dir_name}" \\
        -p INPUT_PATH "\$(readlink -f "input_for_contingency.h5ad")" \\
        -p OUTPUT_DIR "\$(readlink -f "${sample_dir_name}")" \\
        -p TOP_N_POPULATIONS ${params.notebooks.contingency.top_n_populations}

    echo "Checking expected contingency PDF outputs..."

    python - <<'PY_VALIDATE'
from pathlib import Path
import sys

sample_dir = Path("${sample_dir_name}")
cont_dir = sample_dir / "plots" / "contingency"

executed_nb = (
    sample_dir
    / "report"
    / "${sample_dir_name}_contingency_clusters_populations_executed.ipynb"
)

missing = []

if not executed_nb.exists() or executed_nb.stat().st_size == 0:
    missing.append(executed_nb)

pdfs = []

if cont_dir.exists():
    pdfs = sorted(cont_dir.glob("*_counts.pdf"))

if len(pdfs) == 0:
    missing.append(cont_dir / "*_counts.pdf")

if missing:
    print("ERROR: missing expected contingency output(s):")
    for path in missing:
        print(f"  - {path}")

    print()
    print(f"Existing files under {sample_dir}:")
    for path in sorted(sample_dir.rglob("*")):
        if path.is_file():
            print(f"  {path}")

    sys.exit(1)

print("Found contingency outputs:")
print(f"  - executed notebook: {executed_nb}")
print("  - count PDFs:")
for path in pdfs:
    print(f"    {path}")
PY_VALIDATE

    echo "Contingency notebook completed."
    """
}
