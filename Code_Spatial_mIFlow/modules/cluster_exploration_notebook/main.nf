process CLUSTER_EXPLORATION_NOTEBOOK {

    tag "${sample_dir_name}"

    // Reuses the RAPIDS clustering environment; cluster_exploration_setup has no envs/ of its own.
    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/05_CLUSTER_exploration",
        mode: 'copy',
        overwrite: true,
        failOnError: true

    input:
    tuple val(sample_index),
          val(sample_id),
          val(sample_dir_name),
          path(input_adata)

    output:
    path "${sample_dir_name}", emit: cluster_exploration_outputs

    script:
    def NOTEBOOK = "${moduleDir}/cluster_exploration_template.ipynb"

    """
    set -euo pipefail

    mkdir -p "${sample_dir_name}/report"

    EXECUTED_NOTEBOOK="${sample_dir_name}/report/${sample_dir_name}_cluster_exploration_executed.ipynb"

    echo "============================================================"
    echo "CLUSTER EXPLORATION NOTEBOOK"
    echo "============================================================"
    echo "Sample index: ${sample_index}"
    echo "Sample ID:    ${sample_id}"
    echo "Output dir:   ${sample_dir_name}"
    echo "Input RAPIDS AnnData:"
    readlink -f "${input_adata}"
    echo "Notebook:"
    echo "${NOTEBOOK}"
    echo "Executed notebook:"
    echo "\${EXECUTED_NOTEBOOK}"
    echo "Selected normalizations: ${params.notebooks.cluster_exploration.selected_normalizations}"
    echo "Selected resolutions:    ${params.notebooks.cluster_exploration.selected_resolutions}"
    echo "Selected UMAP min_dist:  ${params.notebooks.cluster_exploration.selected_umap_min_dists}"
    echo "============================================================"

    echo "Python executable:"
    command -v python
    python --version

    echo "Checking cluster exploration environment:"
    python - <<'PY_CHECK'
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
import papermill
import ipykernel
import nbformat

print("anndata:", ad.__version__)
print("scanpy:", sc.__version__)
print("pandas:", pd.__version__)
print("numpy:", np.__version__)
print("matplotlib:", matplotlib.__version__)
print("papermill:", papermill.__version__)
print("ipykernel:", ipykernel.__version__)
print("nbformat:", nbformat.__version__)
PY_CHECK

    echo "Ensuring python3 kernel exists in the active RAPIDS environment..."
    # See rapids_clustering_notebook/main.nf for why this is flock-guarded:
    # parallel per-sample tasks share this same conda env's --sys-prefix
    # kernel destination, and ipykernel's rmtree-then-recreate install races
    # across concurrent tasks without serialization.
    KERNEL_SPEC_DIR="\$(python -c 'import sys; print(sys.prefix)')/share/jupyter/kernels/python3"
    mkdir -p "\$(dirname "\${KERNEL_SPEC_DIR}")"
    (
        flock -x 200
        if [[ ! -d "\${KERNEL_SPEC_DIR}" ]]; then
            python -m ipykernel install --sys-prefix --name python3 --display-name "Python 3"
        fi
    ) 200>"\$(dirname "\${KERNEL_SPEC_DIR}")/.install.lock"

    echo "Running cluster exploration notebook with Papermill..."

    python -m papermill "${NOTEBOOK}" "\${EXECUTED_NOTEBOOK}" \
        -k python3 \
        --log-output \
        -p SAMPLE_INDEX "${sample_index}" \
        -p SAMPLE_ID "${sample_id}" \
        -p SAMPLE_DIR_NAME "${sample_dir_name}" \
        -p INPUT_PATH "\$(readlink -f "${input_adata}")" \
        -p OUTPUT_DIR "\$(readlink -f "${sample_dir_name}")" \
        -p SELECTED_NORMALIZATIONS "${params.notebooks.cluster_exploration.selected_normalizations}" \
        -p SELECTED_RESOLUTIONS "${params.notebooks.cluster_exploration.selected_resolutions}" \
        -p SELECTED_UMAP_MIN_DISTS "${params.notebooks.cluster_exploration.selected_umap_min_dists}" \
        -p TOP_N_MARKERS ${params.notebooks.cluster_exploration.top_n_markers} \
        -p MARKER_COLORMAP "${params.notebooks.cluster_exploration.marker_colormap}"

    echo "Validating expected cluster exploration outputs..."

    python - <<'PY_VALIDATE'
from pathlib import Path
import sys

sample_dir = Path("${sample_dir_name}")

expected_files = [
    sample_dir / "tables" / "available_cluster_choices.pdf",
    sample_dir / "tables" / "selected_cluster_choices.pdf",
    sample_dir / "report" / "${sample_dir_name}_cluster_exploration_executed.ipynb",
    sample_dir / "report" / "session_summary.pdf",
]

missing = [
    f for f in expected_files
    if (not f.exists()) or f.stat().st_size == 0
]

if missing:
    print("ERROR: missing expected cluster exploration output(s):")
    for f in missing:
        print(f"  - {f}")

    print()
    print(f"Existing files under {sample_dir}:")
    for path in sorted(sample_dir.rglob("*")):
        if path.is_file():
            print(f"  {path}")

    sys.exit(1)

pdfs = sorted(sample_dir.rglob("*.pdf"))

if not pdfs:
    print("ERROR: no PDF outputs generated.")
    sys.exit(1)

print("Cluster exploration PDFs generated:")
for f in pdfs:
    print(f"  - {f}")
PY_VALIDATE

    echo "Cluster exploration notebook completed."
    """
}
