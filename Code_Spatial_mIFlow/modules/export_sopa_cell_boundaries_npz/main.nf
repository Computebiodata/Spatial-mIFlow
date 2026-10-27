process EXPORT_SOPA_CELL_BOUNDARIES_NPZ {

    tag "${sample_dir_name}"

    conda "${moduleDir}/environment.yml"

    publishDir "${params.outdir}/03_SEGMENTATION_export",
        mode: 'copy',
        overwrite: true,
        failOnError: true

    input:
    tuple val(sample_index),
          val(sample_id),
          val(sample_dir_name),
          path(input_adata)

    output:
    tuple val(sample_index),
          val(sample_id),
          val(sample_dir_name),
          path("${sample_dir_name}/segmentation/cell_boundaries_for_napari.npz"),
          emit: cell_boundaries_npz

    script:
    // outdir may be relative (resolved against launchDir); used as a raw bash
    // path below, so it's absolutized here rather than left to bash, which
    // would otherwise resolve it against this task's own scratch work dir.
    def ABS_OUTDIR = (params.outdir as String).startsWith('/') ? params.outdir : "${launchDir}/${params.outdir}"

    """
    set -euo pipefail

    mkdir -p "${sample_dir_name}/segmentation"

    SOPA_OUTDIR="${ABS_OUTDIR}/02_sopa"
    SAMPLE_ID="${sample_id}"
    INPUT_ADATA_ABS="\$(readlink -f "${input_adata}")"

    export SOPA_OUTDIR
    export SAMPLE_ID
    export INPUT_ADATA_ABS

    echo "============================================================"
    echo "EXPORT SOPA CELL BOUNDARIES NPZ"
    echo "============================================================"
    echo "Sample index: ${sample_index}"
    echo "Sample ID:    ${sample_id}"
    echo "Output dir:   ${sample_dir_name}"
    echo "SOPA outdir:  \${SOPA_OUTDIR}"
    echo "Input AnnData:"
    readlink -f "${input_adata}"
    echo "============================================================"

    echo "Finding SOPA SpatialData .zarr..."

    python - <<'PY_FIND'
from pathlib import Path
import os
import re
import sys

sopa_outdir = Path(os.environ["SOPA_OUTDIR"])
sample_id = os.environ["SAMPLE_ID"]

if not sopa_outdir.exists():
    raise SystemExit(f"SOPA output directory does not exist: {sopa_outdir}")

candidates = []

exact = sopa_outdir / f"{sample_id}.zarr"
if exact.exists():
    candidates.append(exact)

safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", sample_id).strip("_")

for pattern in [
    f"*{sample_id}*.zarr",
    f"*{safe}*.zarr",
    "*.zarr",
]:
    for p in sorted(sopa_outdir.glob(pattern)):
        if p not in candidates:
            candidates.append(p)

if len(candidates) == 0:
    print("Could not find any .zarr directory in:", sopa_outdir)
    raise SystemExit(1)

print("Zarr candidates:")
for p in candidates:
    print("  -", p)

chosen = candidates[0]
Path("sopa_sdata_path.txt").write_text(str(chosen.resolve()) + "\\n")
print("Chosen:", chosen.resolve())
PY_FIND

    SDATA_PATH="\$(cat sopa_sdata_path.txt)"
    OUT_NPZ="${sample_dir_name}/segmentation/cell_boundaries_for_napari.npz"

    echo "Chosen SpatialData:"
    echo "\${SDATA_PATH}"

    echo "Output NPZ:"
    echo "\${OUT_NPZ}"

    python "${moduleDir}/bin/export_sopa_cell_boundaries_npz.py" \
        --sdata-path "\${SDATA_PATH}" \
        --output-npz "\${OUT_NPZ}" \
        --shape-key "${params.segmentation_export.shape_key ?: ''}"

    echo "Validating NPZ..."

    python - <<'PY_VALIDATE'
from pathlib import Path
import numpy as np
import sys

p = Path("${sample_dir_name}/segmentation/cell_boundaries_for_napari.npz")

if not p.exists() or p.stat().st_size == 0:
    raise SystemExit(f"Missing or empty NPZ: {p}")

data = np.load(p, allow_pickle=True)

required = ["polygons", "cell_ids", "shape_key", "bounds"]

missing = [k for k in required if k not in data.files]
if missing:
    raise SystemExit(f"NPZ missing keys: {missing}")

n_polygons = len(data["polygons"])
n_cell_ids = len(data["cell_ids"])

print("NPZ keys:", data.files)
print("n_polygons:", n_polygons)
print("n_cell_ids:", n_cell_ids)
print("shape_key:", data["shape_key"])
print("bounds:", data["bounds"])

if n_polygons == 0:
    raise SystemExit("NPZ contains zero polygons.")

if n_polygons != n_cell_ids:
    raise SystemExit(
        f"Mismatch polygons/cell_ids: {n_polygons} vs {n_cell_ids}"
    )

print("NPZ validation passed.")
PY_VALIDATE

    echo "SOPA cell-boundary NPZ export completed."
    """
}
