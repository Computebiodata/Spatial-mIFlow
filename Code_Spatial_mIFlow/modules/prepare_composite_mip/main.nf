process PREPARE_COMPOSITE_MIP {

    tag "${sample}"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sample), path(input_ome)

    output:
    tuple val(sample),
          path("${sample}/${sample}_plus_composite_MIP.ome.tif"),
          path("${sample}/${sample}_channel_name_map.csv"),
          path("${sample}/${sample}_image_preparation_summary.json"),
          emit: prepared

    script:
    def PREPARE_SKIP = params.image_preparation.prepare_skip_if_complete ? "true" : "false"
    // outdir may be relative (resolved against launchDir); used as a raw bash
    // path below, so it's absolutized here rather than left to bash, which
    // would otherwise resolve it against this task's own scratch work dir.
    def ABS_OUTDIR = (params.outdir as String).startsWith('/') ? params.outdir : "${launchDir}/${params.outdir}"

    """
    set -euo pipefail

    mkdir -p "${sample}"

    PREPARE_SKIP_IF_COMPLETE="${PREPARE_SKIP}"

    EXISTING_DIR="${ABS_OUTDIR}/01_prepared_images/${sample}"
    EXISTING_OME="\${EXISTING_DIR}/${sample}_plus_composite_MIP.ome.tif"
    EXISTING_MAP="\${EXISTING_DIR}/${sample}_channel_name_map.csv"
    EXISTING_SUMMARY="\${EXISTING_DIR}/${sample}_image_preparation_summary.json"

    WORK_OME="${sample}/${sample}_plus_composite_MIP.ome.tif"
    WORK_MAP="${sample}/${sample}_channel_name_map.csv"
    WORK_SUMMARY="${sample}/${sample}_image_preparation_summary.json"

    link_or_copy() {
        src="\$1"
        dst="\$2"

        rm -f "\${dst}"

        if ln "\${src}" "\${dst}" 2>/dev/null; then
            echo "Hard-linked existing file:"
            echo "  \${src}"
            echo "  -> \${dst}"
        else
            echo "Hard-link failed. Copying existing file:"
            echo "  \${src}"
            echo "  -> \${dst}"
            rsync -ah --partial "\${src}" "\${dst}"
        fi
    }

    validate_prepared_ome() {
        ome_path="\$1"

        python - "\${ome_path}" <<'PY'
import sys
from pathlib import Path
from collections import Counter
import tifffile
import xml.etree.ElementTree as ET

p = Path(sys.argv[1])

if not p.exists():
    raise SystemExit(f"ERROR: file does not exist: {p}")

if p.stat().st_size == 0:
    raise SystemExit(f"ERROR: file is empty: {p}")

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

    c = tif.series[0].shape[0]

    root = ET.fromstring(tif.ome_metadata)
    image = root.find(".//{*}Image")

    if image is None:
        raise SystemExit("ERROR: OME metadata contains no Image element.")

    channels = image.findall(".//{*}Channel")
    names = [ch.get("Name") for ch in channels]

    if len(names) != c:
        raise SystemExit(
            f"ERROR: OME channel count {len(names)} != image C {c}"
        )

dups = {k: v for k, v in Counter(names).items() if v > 1}

if dups:
    raise SystemExit(f"ERROR: duplicate channel names: {dups}")

print("Prepared OME validation passed.")
print(f"Path: {p}")
PY
    }

    INPUT_OME_ABS="\$(readlink -f "${input_ome}")"

    echo "============================================================"
    echo "PREPARE COMPOSITE MIP"
    echo "============================================================"
    echo "Sample:                   ${sample}"
    echo "Input OME:                ${input_ome}"
    echo "Input OME absolute:       \${INPUT_OME_ABS}"
    echo "Reuse existing output:    \${PREPARE_SKIP_IF_COMPLETE}"
    echo "Existing prepared OME:    \${EXISTING_OME}"
    echo "============================================================"

    if [[ "\${PREPARE_SKIP_IF_COMPLETE}" == "true" ]] && \
       [[ -s "\${EXISTING_OME}" ]] && \
       [[ -s "\${EXISTING_MAP}" ]] && \
       [[ -s "\${EXISTING_SUMMARY}" ]]; then

        echo "Existing prepared image detected. Validating before reuse..."
        validate_prepared_ome "\${EXISTING_OME}"

        echo "Reusing existing prepared image."
        link_or_copy "\${EXISTING_OME}" "\${WORK_OME}"
        cp -f "\${EXISTING_MAP}" "\${WORK_MAP}"
        cp -f "\${EXISTING_SUMMARY}" "\${WORK_SUMMARY}"

        echo "PREPARE_COMPOSITE_MIP reused existing valid output."
        exit 0
    fi

    echo "No valid existing prepared image found. Running image-preparation script."

    python "${moduleDir}/bin/prepare_composite_mip.py" \
        --input "\${INPUT_OME_ABS}" \
        --output "\${WORK_OME}" \
        --sample-id "${sample}" \
        --mip-start ${params.image_preparation.mip_start ?: 1} \
        --mip-end ${params.image_preparation.mip_end} \
        --mip-name "${params.image_preparation.mip_name ?: 'composite_MIP'}" \
        --drop-channels-enabled "${params.image_preparation.drop_channels_enabled ?: 'false'}" \
        --drop-channel-names "${params.image_preparation.drop_channel_names}" \
        --write-clean-channel-names-to-ome "${params.image_preparation.write_clean_channel_names_to_ome}" \
        --channel-map-csv "\${WORK_MAP}" \
        --summary-json "\${WORK_SUMMARY}"

    echo "PREPARE_COMPOSITE_MIP generated new prepared output."
    """
}
