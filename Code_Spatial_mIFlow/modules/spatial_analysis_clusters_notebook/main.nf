process SPATIAL_ANALYSIS_CLUSTERS_NOTEBOOK {

    tag "${sample_dir_name}"

    cache false

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/07_SPATIAL_cluster_analysis",
        mode: 'copy',
        overwrite: true,
        failOnError: true

    input:
    tuple val(sample_index),
          val(sample_id),
          val(sample_dir_name),
          path(input_adata)

    output:
    path "${sample_dir_name}", emit: spatial_analysis_outputs

    script:
    def NOTEBOOK = "${moduleDir}/spatial_analysis_clusters_template.ipynb"

    """
    set -euo pipefail

    mkdir -p "${sample_dir_name}/report"

    # Temporary local copy.
    # The notebook can patch this copy if needed.
    # The upstream SCIMAP AnnData is not modified.
    cp -L "${input_adata}" "input_for_spatial_analysis.h5ad"

    EXECUTED_NOTEBOOK="${sample_dir_name}/report/${sample_dir_name}_spatial_analysis_clusters_executed.ipynb"

    echo "============================================================"
    echo "SPATIAL ANALYSIS OF LEIDEN CLUSTERS"
    echo "============================================================"
    echo "Sample index: ${sample_index}"
    echo "Sample ID:    ${sample_id}"
    echo "Output dir:   ${sample_dir_name}"
    echo "Input SCIMAP AnnData:"
    readlink -f "${input_adata}"
    echo "Temporary work copy:"
    readlink -f "input_for_spatial_analysis.h5ad"
    echo "Notebook:"
    echo "${NOTEBOOK}"
    echo "Executed notebook:"
    echo "\${EXECUTED_NOTEBOOK}"
    echo "Cluster normalization: ${params.notebooks.spatial.cluster_normalization}"
    echo "Cluster resolution:    ${params.notebooks.spatial.cluster_resolution}"
    echo "============================================================"

    echo "Python executable:"
    command -v python
    python --version

    echo "Checking spatial-analysis notebook environment:"
    python - <<'PY_CHECK'
import squidpy as sq
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib
import scipy
import h5py
import papermill
import ipykernel
import nbformat

print("squidpy:", sq.__version__)
print("anndata:", ad.__version__)
print("pandas:", pd.__version__)
print("numpy:", np.__version__)
print("matplotlib:", matplotlib.__version__)
print("scipy:", scipy.__version__)
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

    echo "Running spatial-analysis notebook with Papermill..."

    python -m papermill \\
        "${NOTEBOOK}" \\
        "\${EXECUTED_NOTEBOOK}" \\
        -k python3 \\
        --log-output \\
        -p SAMPLE_INDEX "${sample_index}" \\
        -p SAMPLE_ID "${sample_id}" \\
        -p SAMPLE_DIR_NAME "${sample_dir_name}" \\
        -p INPUT_PATH "\$(readlink -f "input_for_spatial_analysis.h5ad")" \\
        -p OUTPUT_DIR "\$(readlink -f "${sample_dir_name}")" \\
        -p CLUSTER_NORMALIZATION "${params.notebooks.spatial.cluster_normalization}" \\
        -p CLUSTER_RESOLUTION "${params.notebooks.spatial.cluster_resolution}"

    echo "Checking expected spatial-analysis PDF outputs..."

    python - <<'PY_VALIDATE'
from pathlib import Path
import sys

sample_dir = Path("${sample_dir_name}")

executed_nb = (
    sample_dir
    / "report"
    / "${sample_dir_name}_spatial_analysis_clusters_executed.ipynb"
)

plot_dir = sample_dir / "plots"

missing = []

if not executed_nb.exists() or executed_nb.stat().st_size == 0:
    missing.append(executed_nb)

pdfs = sorted(plot_dir.rglob("*.pdf")) if plot_dir.exists() else []

if len(pdfs) == 0:
    missing.append(plot_dir / "*.pdf")

if missing:
    print("ERROR: missing expected spatial-analysis output(s):")
    for path in missing:
        print(f"  - {path}")

    print()
    print(f"Existing files under {sample_dir}:")
    for path in sorted(sample_dir.rglob("*")):
        if path.is_file():
            print(f"  {path}")

    sys.exit(1)

print("Found spatial-analysis outputs:")
print(f"  - executed notebook: {executed_nb}")
print("  - PDF plots:")
for path in pdfs:
    print(f"    {path}")
PY_VALIDATE

    echo "Spatial cluster analysis completed."
    """
}
