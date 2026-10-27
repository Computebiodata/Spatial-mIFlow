process SCIMAP_AUTOGATING_NOTEBOOK {

    tag "${sample_dir_name}"

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/04_SCIMAP_autogating",
        mode: 'copy',
        overwrite: true,
        failOnError: true

    input:
    tuple val(sample_index),
          val(sample_id),
          val(sample_dir_name),
          path(input_adata)

    output:
    path "${sample_dir_name}", emit: scimap_outputs

    tuple val(sample_index),
          val(sample_id),
          val(sample_dir_name),
          path("${sample_dir_name}/clustering/scimap_gated.h5ad"),
          emit: scimap_adata

    script:
    def NOTEBOOK = "${moduleDir}/scimap_autogating_template.ipynb"

    """
    set -euo pipefail

    mkdir -p "${sample_dir_name}/report"

    EXECUTED_NOTEBOOK="${sample_dir_name}/report/${sample_dir_name}_scimap_autogating_executed.ipynb"

    echo "============================================================"
    echo "SCIMAP AUTOMATIC GATING NOTEBOOK"
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
    echo "============================================================"

    echo "Python executable:"
    command -v python
    python --version

    echo "Checking SCIMAP notebook environment:"
    python - <<'PY_CHECK'
import scimap as sm
import anndata as ad
import papermill
import ipykernel
import nbformat

print("scimap:", getattr(sm, "__version__", "version_not_reported"))
print("anndata:", ad.__version__)
print("papermill:", papermill.__version__)
print("ipykernel:", ipykernel.__version__)
print("nbformat:", nbformat.__version__)
PY_CHECK

    echo "Ensuring python3 kernel exists in the active SCIMAP environment..."
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

    echo "Running SCIMAP notebook with Papermill..."

    python -m papermill "${NOTEBOOK}" "\${EXECUTED_NOTEBOOK}" \
        -k python3 \
        --log-output \
        -p SAMPLE_INDEX "${sample_index}" \
        -p SAMPLE_ID "${sample_id}" \
        -p SAMPLE_DIR_NAME "${sample_dir_name}" \
        -p INPUT_PATH "\$(readlink -f "${input_adata}")" \
        -p OUTPUT_DIR "\$(readlink -f "${sample_dir_name}")" \
        -p SCIMAP_GMM_COMPONENTS ${params.notebooks.scimap.gmm_components} \
        -p SCIMAP_LOG ${params.notebooks.scimap.log} \
        -p SCIMAP_METHOD "${params.notebooks.scimap.method}" \
        -p SCIMAP_RANDOM_STATE ${params.notebooks.scimap.random_state} \
        -p SCIMAP_POSITIVITY_THRESHOLD ${params.notebooks.scimap.positivity_threshold} \
        -p TOP_N_POPULATIONS ${params.notebooks.scimap.top_n_populations} \
        -p TOP_N_SPATIAL_POPULATIONS ${params.notebooks.scimap.top_n_spatial_populations}

    echo "Validating expected SCIMAP outputs..."

    python - <<'PY_VALIDATE'
from pathlib import Path
import sys

sample_dir = Path("${sample_dir_name}")

expected_files = [
    (sample_dir / "clustering" / "scimap_gated.h5ad", "final SCIMAP-gated AnnData"),
    (sample_dir / "tables" / "marker_table.pdf", "marker table PDF"),
    (sample_dir / "tables" / "scimap_gate_summary.pdf", "SCIMAP gate summary PDF"),
    (sample_dir / "tables" / "scimap_marker_positive_summary.pdf", "SCIMAP marker positivity summary PDF"),
    (sample_dir / "tables" / "scimap_population_counts_top30.pdf", "top population table PDF"),
    (sample_dir / "plots" / "scimap_population_counts_top30.pdf", "top population barplot PDF"),
    (sample_dir / "plots" / "scimap_populations_by_n_positive_markers.pdf", "population by marker-count PDF"),
    (sample_dir / "plots" / "spatial_by_n_positive_lineage_markers.pdf", "spatial marker-count PDF"),
    (sample_dir / "plots" / "spatial_top_scimap_populations.pdf", "spatial top-population PDF"),
    (sample_dir / "report" / "session_summary.pdf", "session summary PDF"),
    (sample_dir / "report" / "${sample_dir_name}_scimap_autogating_executed.ipynb", "executed SCIMAP notebook"),
]

missing = []

for path, label in expected_files:
    if (not path.exists()) or path.stat().st_size == 0:
        missing.append((path, label))

if missing:
    print("ERROR: missing expected SCIMAP output(s):")
    for path, label in missing:
        print(f"  - {label}: {path}")

    print()
    print(f"Existing files under {sample_dir}:")
    for path in sorted(sample_dir.rglob("*")):
        if path.is_file():
            print(f"  {path}")

    sys.exit(1)

print("Validated SCIMAP outputs:")
for path, label in expected_files:
    print(f"  - {label}: {path}")
PY_VALIDATE

    echo "SCIMAP automatic gating step completed."
    """
}
