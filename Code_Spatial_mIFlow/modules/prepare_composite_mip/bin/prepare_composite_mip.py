#!/usr/bin/env python

"""
prepare_composite_mip.py

Streaming image-preparation script for large multiplex OME-TIFF images.

Main goals
----------
1. Avoid loading the full image into RAM.
2. Optionally drop technical channels by substring tokens such as TRITC,Cy5.
3. Compute composite_MIP using one-pass channel streaming.
4. Write a SOPA-compatible OME-TIFF.
5. Validate the written OME-TIFF before finalizing.
6. Save a rich channel map CSV and summary JSON.

Memory strategy
---------------
The script does NOT do:

    image_stack = series.asarray()
    kept_stack = image_stack[keep_indices]
    final_stack = np.concatenate(...)

Instead, it uses:

    one TiffWriter.write() call
    data = generator yielding one output tile at a time
    declared shape/dtype up front

The generator:
    - reads each kept channel once
    - updates the composite_MIP accumulator if the channel is in the MIP range
    - splits each kept channel plane into row-major tiles and yields them
    - yields the finished composite_MIP last, also tiled

The output is written tiled and compressed (zstd + horizontal predictor) so
that file size and read performance are comparable to the compressed, tiled
source images instead of ballooning to the fully decompressed pixel volume.

Important
---------
This streaming version requires direct CYX OME-TIFF input.

For images with singleton extra dimensions such as TCYX or ZCYX, do not silently
guess the page-to-channel mapping. Convert them first or add a verified fallback.
"""

import argparse
import json
import os
import re
from pathlib import Path
import xml.etree.ElementTree as ET
from collections import Counter

import numpy as np
import pandas as pd
import tifffile


# ==============================================================================
# BASIC HELPERS
# ==============================================================================

def parse_bool(x) -> bool:
    if isinstance(x, bool):
        return x

    s = str(x).strip().lower()

    if s in {"true", "1", "yes", "y"}:
        return True

    if s in {"false", "0", "no", "n"}:
        return False

    raise ValueError(f"Cannot parse boolean value: {x}")


def parse_csv_list(x):
    if x is None:
        return []

    return [v.strip() for v in str(x).split(",") if v.strip()]


def compact_for_detection(x: str) -> str:
    """
    Compact a channel name for technical-token detection only.

    This is NOT renaming.

    Examples:
        "TRITC Autofluorescence" -> "tritcautofluorescence"
        "CD3 (TRITC = 700)"      -> "cd3tritc700"
        "Cy 5"                   -> "cy5"
    """

    s = str(x).strip().lower()
    s = re.sub(r"[^a-z0-9]+", "", s)
    return s


def simplify_marker_name(name: str) -> str:
    """
    Create a clean marker name for the channel map and optional OME renaming.

    Examples:
        "DAPI (color=white)"      -> "DAPI"
        "PD-1 - Cy5 (color=700)"  -> "PD_1"
        "HLA-DR"                  -> "HLA_DR"
        "Collagen 3"              -> "Collagen_3"
    """

    s = str(name).strip()

    s = re.sub(r"\s*\(.*?\)\s*", "", s)
    s = re.sub(r"\s*\[.*?\]\s*", "", s)

    s = re.split(r"\s-\s", s, maxsplit=1)[0].strip()

    s = re.sub(r"[^a-zA-Z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")

    return s if s else "Unnamed"


def make_unique(names):
    """
    Ensure names are unique.

    Example:
        TRITC, TRITC, TRITC -> TRITC, TRITC_1, TRITC_2
    """

    seen = {}
    out = []

    for name in names:
        if name not in seen:
            seen[name] = 0
            out.append(name)
        else:
            seen[name] += 1
            out.append(f"{name}_{seen[name]}")

    return out


def is_likely_technical_channel(input_channel_index, clean_name: str, dropped: bool) -> bool:
    """
    Mark channels excluded from downstream phenotyping.

    Physical dropping is controlled separately.

    Rules:
      - input channel 0 is considered DAPI/nuclear/technical
      - dropped channels are technical
      - pure fluorophore-only clean names are technical
      - composite_MIP is added separately and marked technical in the MIP row
    """

    if input_channel_index == 0:
        return True

    if dropped:
        return True

    s = str(clean_name).strip().lower()

    if s == "":
        return True

    technical_terms = {
        "composite",
        "composite_mip",
        "mip",
    }

    if any(term in s for term in technical_terms):
        return True

    pure_fluorophore_names = {
        "tritc",
        "cy5",
        "cy3",
        "cy7",
        "fitc",
    }

    base = re.sub(r"_\d+$", "", s)

    if base in pure_fluorophore_names:
        return True

    return False


# ==============================================================================
# OME / TIFF HELPERS
# ==============================================================================

def parse_ome_xml(ome_xml):
    if ome_xml is None:
        raise ValueError("OME metadata is missing.")

    if isinstance(ome_xml, bytes):
        ome_xml = ome_xml.decode("utf-8")

    return ET.fromstring(ome_xml)


def extract_channel_names_from_ome(tif: tifffile.TiffFile):
    root = parse_ome_xml(tif.ome_metadata)

    first_image = root.find(".//{*}Image")

    if first_image is None:
        raise ValueError("Could not find an <Image> block in OME-XML metadata.")

    channels = first_image.findall(".//{*}Channel")

    if len(channels) == 0:
        raise ValueError("Could not find any <Channel> entries in OME-XML metadata.")

    channel_ids = [
        ch.get("ID", f"Channel_{i}")
        for i, ch in enumerate(channels)
    ]

    channel_names = [
        ch.get("Name", f"Unnamed_{i}")
        for i, ch in enumerate(channels)
    ]

    return channel_ids, channel_names


def inspect_tiff_layout(input_path: Path):
    """
    Inspect image metadata without reading the full image.

    This streaming path requires:
        - OME-TIFF
        - one series
        - direct CYX axes
        - one page per channel
    """

    with tifffile.TiffFile(input_path) as tif:
        if not tif.is_ome:
            raise ValueError(f"Input is not detected as OME-TIFF: {input_path}")

        if tif.ome_metadata is None:
            raise ValueError(f"Input OME metadata is missing: {input_path}")

        if len(tif.series) != 1:
            raise ValueError(
                f"Expected one OME image series. Found {len(tif.series)} series. "
                "Streaming path assumes one CYX image series."
            )

        series = tif.series[0]
        axes = series.axes
        shape = tuple(series.shape)
        dtype = np.dtype(series.dtype)

        if axes != "CYX":
            raise ValueError(
                f"Streaming mode currently requires direct CYX axes. "
                f"Got axes={axes}, shape={shape}, file={input_path}. "
                "For singleton dimensions such as TCYX or ZCYX, convert first "
                "or add a verified fallback."
            )

        if len(shape) != 3:
            raise ValueError(
                f"Expected 3D CYX image. Got axes={axes}, shape={shape}"
            )

        n_channels, height, width = shape

        channel_ids, original_names = extract_channel_names_from_ome(tif)

        if len(original_names) != n_channels:
            raise ValueError(
                "Number of OME channel names does not match C dimension.\n"
                f"C dimension: {n_channels}\n"
                f"OME channel names: {len(original_names)}"
            )

        if len(series.pages) < n_channels:
            raise ValueError(
                "Streaming mode requires at least one TIFF page per channel.\n"
                f"Image channels: {n_channels}\n"
                f"Series pages: {len(series.pages)}"
            )

        first_page_shape = tuple(series.pages[0].shape)

        if first_page_shape != (height, width):
            raise ValueError(
                f"First TIFF page shape {first_page_shape} does not match expected YX {(height, width)}"
            )

    return {
        "axes": axes,
        "shape": shape,
        "dtype": dtype,
        "n_channels": int(n_channels),
        "height": int(height),
        "width": int(width),
        "channel_ids": channel_ids,
        "original_names": original_names,
    }


def validate_written_ome(path: Path, expected_channel_names):
    """
    Reopen written OME-TIFF and validate SOPA compatibility.
    """

    path = Path(path)

    with tifffile.TiffFile(path) as tif:
        if not tif.is_ome:
            raise ValueError(f"Written file is not detected as OME-TIFF: {path}")

        if tif.ome_metadata is None:
            raise ValueError(f"Written file has no OME metadata: {path}")

        if len(tif.series) != 1:
            raise ValueError(
                f"Written file should contain one series. Found {len(tif.series)} series."
            )

        series = tif.series[0]

        if series.axes != "CYX":
            raise ValueError(
                f"Written OME-TIFF must be CYX. "
                f"Got axes={series.axes}, shape={series.shape}, path={path}"
            )

        if series.shape[0] != len(expected_channel_names):
            raise ValueError(
                f"Written channel count mismatch. "
                f"Image C={series.shape[0]}, "
                f"expected names={len(expected_channel_names)}, "
                f"path={path}"
            )

        root = parse_ome_xml(tif.ome_metadata)
        image = root.find(".//{*}Image")

        if image is None:
            raise ValueError(f"No Image block found in written OME metadata: {path}")

        channels = image.findall(".//{*}Channel")
        observed_names = [ch.get("Name") for ch in channels]

    if observed_names != list(expected_channel_names):
        raise ValueError(
            "Written OME channel names do not match expected names.\n"
            f"Expected: {expected_channel_names}\n"
            f"Observed: {observed_names}"
        )

    duplicates = {k: v for k, v in Counter(observed_names).items() if v > 1}

    if duplicates:
        raise ValueError(
            f"Written OME-TIFF has duplicate channel names: {duplicates}. "
            "SOPA requires unique channel names."
        )

    return True


# ==============================================================================
# CHANNEL DECISION LOGIC
# ==============================================================================

def should_drop_channel_by_original_name(original_name: str, drop_tokens, drop_enabled: bool):
    """
    Physical dropping logic.

    If enabled, a channel is dropped when any compacted token appears anywhere
    inside the compacted original OME channel name.

    Example with token TRITC:
        TRITC                         -> drop
        TRITC Autofluorescence         -> drop
        CD3 (TRITC = 700)              -> drop if enabled

    Therefore only enable this for images where these tokens identify unwanted
    technical/autofluorescence channels.
    """

    if not drop_enabled:
        return False, ""

    compact_name = compact_for_detection(original_name)

    for token in drop_tokens:
        compact_token = compact_for_detection(token)

        if compact_token and compact_token in compact_name:
            return True, f"token '{compact_token}' found in original channel name"

    return False, ""


def cast_nullable_int_columns(df, columns):
    for col in columns:
        df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    return df


def iter_plane_tiles(plane, tile_shape):
    """
    Yield row-major tiles from a single 2D plane.

    Edge tiles are yielded at their true, smaller size. tifffile zero-pads
    them to the declared tile shape internally, so no manual padding here.
    """

    tile_height, tile_width = tile_shape
    height, width = plane.shape

    for y in range(0, height, tile_height):
        y_end = min(y + tile_height, height)
        for x in range(0, width, tile_width):
            x_end = min(x + tile_width, width)
            yield plane[y:y_end, x:x_end]


def stream_output_planes(input_path, keep_indices, mip_source_indices, original_names, tile_shape, maxworkers):
    """
    Single-pass streaming generator.

    Reads each kept channel once.
    Updates the MIP accumulator when the channel is in mip_source_indices.
    Yields each kept channel plane to the writer as row-major tiles.
    Yields the finished composite_MIP plane last, also tiled.
    """

    mip_source_set = set(mip_source_indices)

    with tifffile.TiffFile(input_path) as tif:
        if len(tif.series) != 1:
            raise ValueError(
                f"Expected one series, found {len(tif.series)}."
            )

        series = tif.series[0]

        if series.axes != "CYX":
            raise ValueError(
                f"Streaming path requires axes == CYX. Got axes={series.axes}."
            )

        accumulator = None

        for output_i, input_i in enumerate(keep_indices):
            print(
                f"  Writing output C{output_i + 1}: "
                f"input {input_i:02d} / C{input_i + 1}: {original_names[input_i]}"
            )

            plane = series.pages[input_i].asarray(maxworkers=maxworkers)

            if input_i in mip_source_set:
                if accumulator is None:
                    accumulator = plane.copy()
                else:
                    np.maximum(accumulator, plane, out=accumulator)

            yield from iter_plane_tiles(plane, tile_shape)

            del plane

        if accumulator is None:
            raise ValueError("No channels left for composite MIP after drop filtering.")

        print(f"  Writing output C{len(keep_indices) + 1}: composite_MIP")
        yield from iter_plane_tiles(accumulator, tile_shape)


# ==============================================================================
# CLI
# ==============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Prepare one OME-TIFF for SOPA using streaming channel reading, "
            "optional channel dropping, and composite MIP generation."
        )
    )

    parser.add_argument("--input", required=True, help="Input OME-TIFF image.")
    parser.add_argument("--output", required=True, help="Output prepared OME-TIFF.")
    parser.add_argument("--sample-id", required=True, help="Sample ID.")
    parser.add_argument("--channel-map-csv", required=True, help="Output channel map CSV.")
    parser.add_argument("--summary-json", required=True, help="Output summary JSON.")

    parser.add_argument(
        "--mip-name",
        default="composite_MIP",
        help="Name of the appended composite MIP channel."
    )

    parser.add_argument(
        "--mip-start",
        type=int,
        default=1,
        help="First input channel index for MIP, 0-based inclusive."
    )

    parser.add_argument(
        "--mip-end",
        type=int,
        default=None,
        help="End input channel index for MIP, 0-based exclusive."
    )

    parser.add_argument(
        "--drop-channels-enabled",
        default="false",
        help="true/false. If true, drop channels matching --drop-channel-names tokens."
    )

    parser.add_argument(
        "--drop-channel-names",
        default="",
        help="Comma-separated tokens for substring-based dropping. Example: TRITC,Cy5"
    )

    parser.add_argument(
        "--write-clean-channel-names-to-ome",
        default="true",
        help=(
            "true/false. If true, write unique cleaned names in output OME. "
            "If false, keep original names for remaining channels; validation fails if duplicated."
        )
    )

    return parser.parse_args()


# ==============================================================================
# MAIN
# ==============================================================================

# Multiple of 16, as required by the TIFF tile spec. Large enough to keep the
# tile count (and IFD overhead) manageable on ~40k x 40k images, small enough
# to still allow cheap random-access reads of sub-regions downstream.
OUTPUT_TILE_SHAPE = (1024, 1024)


def available_worker_count():
    """Cores actually usable by this task, not the whole node.

    os.cpu_count()/multiprocessing.cpu_count() report the physical node's
    total core count regardless of SLURM/cgroup restriction (verified on
    this cluster: a job with --cpus-per-task=2 still reports 384). tifffile
    derives its own internal thread-pool size from that same uncapped
    count, so leaving maxworkers unset causes it to spin up far more
    threads than the job actually has cores for. sched_getaffinity reflects
    the real cgroup-restricted CPU set and is what we want instead.
    """
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


WORKER_COUNT = available_worker_count()


def main():
    args = parse_args()

    input_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()
    channel_map_csv = Path(args.channel_map_csv).resolve()
    summary_json = Path(args.summary_json).resolve()

    tmp_output_path = Path(str(output_path) + ".tmp")

    drop_channels_enabled = parse_bool(args.drop_channels_enabled)
    write_clean_channel_names_to_ome = parse_bool(args.write_clean_channel_names_to_ome)
    drop_tokens = parse_csv_list(args.drop_channel_names)

    if not input_path.exists():
        raise FileNotFoundError(f"Input image does not exist: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    channel_map_csv.parent.mkdir(parents=True, exist_ok=True)
    summary_json.parent.mkdir(parents=True, exist_ok=True)

    if tmp_output_path.exists():
        tmp_output_path.unlink()

    print("============================================================")
    print("PREPARE OME-TIFF STREAMING MODE")
    print("============================================================")
    print(f"Sample ID:                         {args.sample_id}")
    print(f"Input image:                       {input_path}")
    print(f"Temporary output image:            {tmp_output_path}")
    print(f"Final output image:                {output_path}")
    print(f"Channel map CSV:                   {channel_map_csv}")
    print(f"Summary JSON:                      {summary_json}")
    print(f"MIP name:                          {args.mip_name}")
    print(f"MIP input range:                   {args.mip_start}:{args.mip_end}")
    print(f"Drop channels enabled:             {drop_channels_enabled}")
    print(f"Drop channel tokens:               {drop_tokens}")
    print(f"Write clean channel names to OME:  {write_clean_channel_names_to_ome}")
    print("============================================================")

    layout = inspect_tiff_layout(input_path)

    axes = layout["axes"]
    input_shape = layout["shape"]
    input_dtype = layout["dtype"]
    n_channels = layout["n_channels"]
    height = layout["height"]
    width = layout["width"]
    channel_ids = layout["channel_ids"]
    original_names = layout["original_names"]

    print(f"\nInput axes:  {axes}")
    print(f"Input shape: {input_shape}")
    print(f"Input dtype: {input_dtype}")
    print(f"Channels:    {n_channels}")

    print("\nOriginal channel names:")
    for i, name in enumerate(original_names):
        print(f"  {i:02d} / C{i + 1}: {name}")

    mip_start = args.mip_start
    mip_end = args.mip_end if args.mip_end is not None else n_channels

    if mip_start < 0 or mip_start >= n_channels:
        raise ValueError(
            f"Invalid --mip-start {mip_start}. Image has {n_channels} channels."
        )

    if mip_end <= mip_start or mip_end > n_channels:
        raise ValueError(
            f"Invalid --mip-end {mip_end}. Must be > mip_start and <= {n_channels}."
        )

    clean_names_raw = [simplify_marker_name(name) for name in original_names]
    clean_names_unique = make_unique(clean_names_raw)

    print("\nClean unique names:")
    for i, name in enumerate(clean_names_unique):
        print(f"  {i:02d} / C{i + 1}: {original_names[i]} -> {name}")

    drop_flags = []
    drop_reasons = []

    for original in original_names:
        flag, reason = should_drop_channel_by_original_name(
            original_name=original,
            drop_tokens=drop_tokens,
            drop_enabled=drop_channels_enabled,
        )
        drop_flags.append(flag)
        drop_reasons.append(reason)

    keep_indices = [i for i in range(n_channels) if not drop_flags[i]]
    dropped_indices = [i for i in range(n_channels) if drop_flags[i]]

    if len(keep_indices) == 0:
        raise ValueError("All channels were dropped. Refusing to write empty image.")

    print("\nChannels dropped from output:")
    if dropped_indices:
        for i in dropped_indices:
            print(
                f"  {i:02d} / C{i + 1}: {original_names[i]} -> "
                f"{clean_names_unique[i]} ({drop_reasons[i]})"
            )
    else:
        print("  None")

    mip_source_indices = [
        i for i in range(mip_start, mip_end)
        if not drop_flags[i]
    ]

    if len(mip_source_indices) == 0:
        raise ValueError(
            "No channels left for composite MIP after applying drop-channel rules."
        )

    print("\nMIP source channels after drop filtering:")
    for i in mip_source_indices:
        print(f"  {i:02d} / C{i + 1}: {original_names[i]} -> {clean_names_unique[i]}")

    if write_clean_channel_names_to_ome:
        output_channel_names = [clean_names_unique[i] for i in keep_indices] + [args.mip_name]
    else:
        output_channel_names = [original_names[i] for i in keep_indices] + [args.mip_name]

    duplicates = {k: v for k, v in Counter(output_channel_names).items() if v > 1}

    if duplicates:
        raise ValueError(
            "Output OME channel names are not unique.\n"
            f"Duplicates: {duplicates}\n"
            "SOPA requires unique channel names.\n"
            "Fix options:\n"
            "  1. enable channel dropping for duplicated technical channels\n"
            "  2. set --write-clean-channel-names-to-ome true\n"
            "  3. manually correct channel names before running the pipeline"
        )

    final_shape = (len(output_channel_names), height, width)

    print("\nOutput channel names:")
    for i, name in enumerate(output_channel_names):
        print(f"  {i:02d} / C{i + 1}: {name}")

    print("\nWriting prepared OME-TIFF in streaming mode...")
    print(f"Final shape: {final_shape}")
    print(f"Final dtype: {input_dtype}")
    print(f"Tile shape: {OUTPUT_TILE_SHAPE}")
    print("Compression: zstd (level 10) with horizontal predictor")
    print(f"Worker threads (decode + encode): {WORKER_COUNT}")

    with tifffile.TiffWriter(tmp_output_path, bigtiff=True, ome=True) as tw:
        tw.write(
            stream_output_planes(
                input_path=input_path,
                keep_indices=keep_indices,
                mip_source_indices=mip_source_indices,
                original_names=original_names,
                tile_shape=OUTPUT_TILE_SHAPE,
                maxworkers=WORKER_COUNT,
            ),
            shape=final_shape,
            dtype=input_dtype,
            metadata={
                "axes": "CYX",
                "Channel": {"Name": output_channel_names},
            },
            tile=OUTPUT_TILE_SHAPE,
            compression="zstd",
            compressionargs={"level": 10},
            predictor=True,
            maxworkers=WORKER_COUNT,
        )

    print("\nValidating temporary written OME-TIFF...")
    validate_written_ome(tmp_output_path, output_channel_names)
    print("Temporary OME-TIFF is valid.")

    print("\nFinalizing output file with atomic rename...")
    os.replace(tmp_output_path, output_path)

    print("\nValidating final OME-TIFF...")
    validate_written_ome(output_path, output_channel_names)
    print("Final OME-TIFF is valid and SOPA-compatible.")

    # --------------------------------------------------------------------------
    # CHANNEL MAP
    # --------------------------------------------------------------------------

    input_to_output_index = {
        input_i: output_i
        for output_i, input_i in enumerate(keep_indices)
    }

    channel_map = pd.DataFrame({
        "sample_id": args.sample_id,
        "input_channel_index": list(range(n_channels)),
        "input_channel_number_1based": list(range(1, n_channels + 1)),
        "output_channel_index": [
            input_to_output_index.get(i, None)
            for i in range(n_channels)
        ],
        "channel_id": channel_ids,
        "channel_original": original_names,
        "channel_clean_raw": clean_names_raw,
        "channel_clean": clean_names_unique,
        "dropped_from_output": drop_flags,
        "drop_reason": drop_reasons,
        "exclude_from_phenotyping": [
            is_likely_technical_channel(
                input_channel_index=i,
                clean_name=clean_names_unique[i],
                dropped=drop_flags[i],
            )
            for i in range(n_channels)
        ],
        "used_for_composite_mip": [
            i in mip_source_indices
            for i in range(n_channels)
        ],
    })

    mip_channel_index = len(keep_indices)

    mip_row = pd.DataFrame({
        "sample_id": [args.sample_id],
        "input_channel_index": [None],
        "input_channel_number_1based": [None],
        "output_channel_index": [mip_channel_index],
        "channel_id": [f"Channel:0:{mip_channel_index}"],
        "channel_original": [args.mip_name],
        "channel_clean_raw": [args.mip_name],
        "channel_clean": [args.mip_name],
        "dropped_from_output": [False],
        "drop_reason": [""],
        "exclude_from_phenotyping": [True],
        "used_for_composite_mip": [False],
    })

    channel_map = pd.concat([channel_map, mip_row], ignore_index=True)

    channel_map = cast_nullable_int_columns(
        channel_map,
        [
            "input_channel_index",
            "input_channel_number_1based",
            "output_channel_index",
        ],
    )

    print("\nSaving channel map CSV...")
    channel_map.to_csv(channel_map_csv, index=False, na_rep="")

    excluded_clean_names = (
        channel_map.loc[
            channel_map["exclude_from_phenotyping"].astype(bool),
            "channel_clean"
        ]
        .astype(str)
        .tolist()
    )

    summary = {
        "sample_id": args.sample_id,
        "input_path": str(input_path),
        "output_path": str(output_path),
        "channel_map_csv": str(channel_map_csv),

        "streaming_mode": True,
        "full_image_loaded_into_ram": False,

        "input_axes": axes,
        "input_shape": list(input_shape),
        "input_dtype": str(input_dtype),

        "n_input_channels": int(n_channels),
        "n_kept_input_channels": int(len(keep_indices)),
        "n_dropped_input_channels": int(len(dropped_indices)),
        "n_output_channels": int(final_shape[0]),

        "output_shape": list(final_shape),
        "output_dtype": str(input_dtype),
        "output_tile_shape": list(OUTPUT_TILE_SHAPE),
        "output_compression": "zstd",
        "output_compression_level": 10,
        "output_predictor": True,

        "mip_name": args.mip_name,
        "mip_channel_index": int(mip_channel_index),
        "mip_start": int(mip_start),
        "mip_end": int(mip_end),

        "mip_source_input_channel_indices_after_drop": [
            int(i) for i in mip_source_indices
        ],
        "mip_source_channel_original_names_after_drop": [
            original_names[i] for i in mip_source_indices
        ],
        "mip_source_channel_clean_names_after_drop": [
            clean_names_unique[i] for i in mip_source_indices
        ],

        "drop_channels_enabled": bool(drop_channels_enabled),
        "drop_channel_tokens": drop_tokens,
        "drop_matching_mode": "substring_on_compacted_original_channel_name",
        "dropped_input_channel_indices": [
            int(i) for i in dropped_indices
        ],
        "dropped_channel_original_names": [
            original_names[i] for i in dropped_indices
        ],
        "dropped_channel_clean_names": [
            clean_names_unique[i] for i in dropped_indices
        ],
        "drop_reasons": [
            drop_reasons[i] for i in dropped_indices
        ],

        "kept_input_channel_indices": [
            int(i) for i in keep_indices
        ],
        "kept_channel_original_names": [
            original_names[i] for i in keep_indices
        ],
        "kept_channel_clean_names": [
            clean_names_unique[i] for i in keep_indices
        ],

        "write_clean_channel_names_to_ome": bool(write_clean_channel_names_to_ome),
        "channel_names_cleaned_inside_ome": bool(write_clean_channel_names_to_ome),
        "output_channel_names_written_to_ome": output_channel_names,
        "output_channel_names_unique": len(output_channel_names) == len(set(output_channel_names)),

        "excluded_from_phenotyping_clean_names": excluded_clean_names,

        "sopa_compatible_validation_passed": True,
    }

    print("\nSaving summary JSON...")
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)

    print("============================================================")
    print("DONE")
    print("============================================================")


if __name__ == "__main__":
    main()
