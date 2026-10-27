#!/usr/bin/env python

"""
build_sopa_samplesheet.py

Build a SOPA samplesheet from validated prepared-image manifests.

Input manifest CSV:
    sample,prepared_image,channel_map,summary

Output SOPA samplesheet:
    sample,data_path

This script validates:
    - prepared image exists
    - image is OME-TIFF
    - axes are CYX
    - channel names are unique
    - required SOPA channels can be found
"""

import argparse
import csv
import json
import re
from pathlib import Path
from collections import Counter

import tifffile
import xml.etree.ElementTree as ET


def normalize_name(x: str) -> str:
    s = str(x).strip().lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def parse_csv_list(x):
    if x is None:
        return []

    return [v.strip() for v in str(x).split(",") if v.strip()]


def parse_ome_xml(ome_xml):
    if ome_xml is None:
        raise ValueError("OME metadata is missing.")

    if isinstance(ome_xml, bytes):
        ome_xml = ome_xml.decode("utf-8")

    return ET.fromstring(ome_xml)


def get_ome_channel_names(path: Path):
    with tifffile.TiffFile(path) as tif:
        if not tif.is_ome:
            raise ValueError(f"Not an OME-TIFF: {path}")

        if tif.ome_metadata is None:
            raise ValueError(f"OME metadata is missing: {path}")

        if len(tif.series) != 1:
            raise ValueError(f"Expected one series. Found {len(tif.series)} in {path}")

        series = tif.series[0]
        axes = series.axes
        shape = series.shape

        if axes != "CYX":
            raise ValueError(
                f"Expected CYX axes for SOPA input. "
                f"Got axes={axes}, shape={shape}, file={path}"
            )

        root = parse_ome_xml(tif.ome_metadata)
        image = root.find(".//{*}Image")

        if image is None:
            raise ValueError(f"No Image block found in OME metadata: {path}")

        channels = image.findall(".//{*}Channel")
        names = [ch.get("Name") for ch in channels]

        if len(names) != shape[0]:
            raise ValueError(
                f"OME channel count does not match image C dimension. "
                f"OME names={len(names)}, shape={shape}, file={path}"
            )

    return names, axes, shape


def resolve_required_channel(requested, available_names):
    """
    Resolve requested channel against available OME names.

    Priority:
        1. exact string match
        2. normalized exact match
        3. unique normalized loose match
    """

    if requested in available_names:
        return requested, "exact"

    requested_norm = normalize_name(requested)

    norm_to_names = {}
    for name in available_names:
        norm_to_names.setdefault(normalize_name(name), []).append(name)

    if requested_norm in norm_to_names:
        matches = norm_to_names[requested_norm]

        if len(matches) == 1:
            return matches[0], "normalized_exact"

        raise ValueError(
            f"Requested channel '{requested}' has ambiguous normalized matches: {matches}"
        )

    loose_matches = [
        name for name in available_names
        if requested_norm in normalize_name(name)
        or normalize_name(name) in requested_norm
    ]

    loose_matches = sorted(set(loose_matches))

    if len(loose_matches) == 1:
        return loose_matches[0], "normalized_loose"

    if len(loose_matches) > 1:
        raise ValueError(
            f"Requested channel '{requested}' is ambiguous after loose matching: "
            f"{loose_matches}"
        )

    raise ValueError(
        f"Requested channel '{requested}' not found. "
        f"Available channels: {available_names}"
    )


def validate_prepared_image(path: Path, required_channels):
    names, axes, shape = get_ome_channel_names(path)

    duplicates = {k: v for k, v in Counter(names).items() if v > 1}

    if duplicates:
        raise ValueError(
            f"Duplicate OME channel names detected in {path}: {duplicates}. "
            "SOPA requires unique channel names."
        )

    resolved = []

    for requested in required_channels:
        exact_name, method = resolve_required_channel(requested, names)
        resolved.append({
            "requested": requested,
            "resolved": exact_name,
            "method": method,
        })

    return {
        "path": str(path),
        "axes": axes,
        "shape": list(shape),
        "channel_names": names,
        "resolved_required_channels": resolved,
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build and validate SOPA samplesheet from prepared image manifests."
    )

    parser.add_argument(
        "--prepared-manifest",
        nargs="+",
        required=True,
        help="One or more manifest CSV files produced by PUBLISH_PREPARED_IMAGE."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output SOPA samplesheet CSV."
    )

    parser.add_argument(
        "--resolved-channels-json",
        required=True,
        help="Output JSON with resolved channel names."
    )

    parser.add_argument(
        "--required-channels",
        default="composite_MIP,DAPI",
        help="Comma-separated required SOPA channels."
    )

    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail on invalid prepared images."
    )

    return parser.parse_args()


def read_manifest_rows(manifest_paths):
    rows = []

    for manifest in manifest_paths:
        manifest = Path(manifest).resolve()

        if not manifest.exists():
            raise FileNotFoundError(f"Prepared manifest not found: {manifest}")

        with open(manifest, newline="") as f:
            reader = csv.DictReader(f)

            required_cols = {"sample", "prepared_image"}
            missing = required_cols - set(reader.fieldnames or [])

            if missing:
                raise ValueError(f"Manifest {manifest} missing columns: {missing}")

            for row in reader:
                rows.append(row)

    return rows


def main():
    args = parse_args()

    output_path = Path(args.output).resolve()
    resolved_json = Path(args.resolved_channels_json).resolve()
    required_channels = parse_csv_list(args.required_channels)

    manifest_rows = read_manifest_rows(args.prepared_manifest)

    if len(manifest_rows) == 0:
        raise ValueError("No rows found in prepared manifests.")

    sopa_rows = []
    validation_records = []

    for row in manifest_rows:
        sample = str(row["sample"]).strip()
        prepared_image = Path(row["prepared_image"]).resolve()

        if sample == "":
            raise ValueError("Empty sample name in prepared manifest.")

        if not prepared_image.exists():
            msg = f"Prepared image does not exist for sample '{sample}': {prepared_image}"
            if args.strict:
                raise FileNotFoundError(msg)
            print(f"WARNING: {msg}")
            continue

        print("Validating prepared image for SOPA:")
        print(f"  sample: {sample}")
        print(f"  image:  {prepared_image}")

        record = validate_prepared_image(
            prepared_image,
            required_channels=required_channels,
        )

        record["sample"] = sample

        print(f"  axes:   {record['axes']}")
        print(f"  shape:  {record['shape']}")
        print("  resolved channels:")
        for x in record["resolved_required_channels"]:
            print(f"    {x['requested']} -> {x['resolved']} ({x['method']})")

        sopa_rows.append({
            "sample": sample,
            "data_path": str(prepared_image),
        })

        validation_records.append(record)

    if len(sopa_rows) == 0:
        raise ValueError("No valid prepared images found. SOPA samplesheet not created.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    resolved_json.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample", "data_path"])
        writer.writeheader()
        writer.writerows(sopa_rows)

    with open(resolved_json, "w") as f:
        json.dump(
            {
                "required_channels_requested": required_channels,
                "samples": validation_records,
            },
            f,
            indent=2,
        )

    print()
    print("SOPA samplesheet written:")
    print(output_path)
    print()
    print("Resolved channels JSON written:")
    print(resolved_json)
    print()
    print("SOPA samplesheet content:")
    with open(output_path) as f:
        print(f.read())


if __name__ == "__main__":
    main()
