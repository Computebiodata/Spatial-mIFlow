process PUBLISH_PREPARED_IMAGE {

    tag "${sample}"

    cache false

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample),
          path(prepared_ome),
          path(channel_map),
          path(summary_json)

    output:
    path "${sample}_prepared_manifest.csv", emit: manifest

    script:
    // outdir may be relative (resolved against launchDir); used as a raw bash
    // path below, so it's absolutized here rather than left to bash, which
    // would otherwise resolve it against this task's own scratch work dir.
    def ABS_OUTDIR = (params.outdir as String).startsWith('/') ? params.outdir : "${launchDir}/${params.outdir}"

    """
    set -euo pipefail

    OUTDIR="${ABS_OUTDIR}/01_prepared_images/${sample}"

    OUT_OME="\${OUTDIR}/${sample}_plus_composite_MIP.ome.tif"
    OUT_MAP="\${OUTDIR}/${sample}_channel_name_map.csv"
    OUT_SUMMARY="\${OUTDIR}/${sample}_image_preparation_summary.json"

    mkdir -p "\${OUTDIR}"

    SRC_OME="\$(readlink -f "${prepared_ome}")"
    SRC_MAP="\$(readlink -f "${channel_map}")"
    SRC_SUMMARY="\$(readlink -f "${summary_json}")"

    echo "============================================================"
    echo "PUBLISH PREPARED IMAGE"
    echo "============================================================"
    echo "Sample:      ${sample}"
    echo "Source OME:  \${SRC_OME}"
    echo "Output OME:  \${OUT_OME}"
    echo "============================================================"

    link_or_copy() {
        src="\$1"
        dst="\$2"

        rm -f "\${dst}"

        if ln "\${src}" "\${dst}" 2>/dev/null; then
            echo "Hard-linked: \${dst}"
        else
            echo "Hard-link failed. Falling back to rsync copy."
            tmp="\${dst}.tmp.\$\$"
            rm -f "\${tmp}"
            rsync -ah --partial --progress "\${src}" "\${tmp}"
            mv -f "\${tmp}" "\${dst}"
            echo "Copied: \${dst}"
        fi
    }

    link_or_copy "\${SRC_OME}" "\${OUT_OME}"

    cp -f "\${SRC_MAP}" "\${OUT_MAP}"
    cp -f "\${SRC_SUMMARY}" "\${OUT_SUMMARY}"

    echo "Validating published OME-TIFF..."

    python - "\${OUT_OME}" <<'PY'
import sys
from pathlib import Path
from collections import Counter
import tifffile
import xml.etree.ElementTree as ET

p = Path(sys.argv[1])

with tifffile.TiffFile(p) as tif:
    if not tif.is_ome:
        raise SystemExit(f"ERROR: not OME-TIFF: {p}")

    if tif.ome_metadata is None:
        raise SystemExit(f"ERROR: missing OME metadata: {p}")

    if len(tif.series) != 1:
        raise SystemExit(f"ERROR: expected one series, found {len(tif.series)}")

    if tif.series[0].axes != "CYX":
        raise SystemExit(
            f"ERROR: expected CYX, got axes={tif.series[0].axes}, "
            f"shape={tif.series[0].shape}"
        )

    root = ET.fromstring(tif.ome_metadata)
    image = root.find(".//{*}Image")

    if image is None:
        raise SystemExit("ERROR: OME metadata contains no Image element.")

    channels = image.findall(".//{*}Channel")
    names = [ch.get("Name") for ch in channels]

    if len(names) != tif.series[0].shape[0]:
        raise SystemExit(
            f"ERROR: OME channel count {len(names)} != image C {tif.series[0].shape[0]}"
        )

dups = {k: v for k, v in Counter(names).items() if v > 1}

if dups:
    raise SystemExit(f"ERROR: duplicate channel names: {dups}")

print("Published OME-TIFF validation passed.")
print(f"Path: {p}")
PY

    cat > "${sample}_prepared_manifest.csv" <<EOF2
sample,prepared_image,channel_map,summary
${sample},\${OUT_OME},\${OUT_MAP},\${OUT_SUMMARY}
EOF2

    echo "Prepared manifest:"
    cat "${sample}_prepared_manifest.csv"
    """
}
