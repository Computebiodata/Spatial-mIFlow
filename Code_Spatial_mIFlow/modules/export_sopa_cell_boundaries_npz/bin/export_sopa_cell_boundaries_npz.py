#!/usr/bin/env python3

from pathlib import Path
import argparse
import numpy as np

from spatialdata import read_zarr
from shapely.geometry import Polygon, MultiPolygon, GeometryCollection


def pick_shape_key(sdata, preferred=None):
    shape_keys = list(sdata.shapes.keys())
    print("Available shape layers:", shape_keys)

    if preferred is not None and preferred != "":
        if preferred not in shape_keys:
            raise ValueError(f"{preferred!r} not found. Available: {shape_keys}")
        return preferred

    priority = [
        "cell_boundaries",
        "cell_boundary",
        "cells",
        "cell",
        "segmentation",
        "boundaries",
        "boundary",
    ]

    for word in priority:
        for key in shape_keys:
            if word.lower() in key.lower():
                return key

    if len(shape_keys) == 1:
        return shape_keys[0]

    raise ValueError(f"Could not choose shape layer automatically. Available: {shape_keys}")


def geometry_to_polygons_xy(geom):
    """
    Return list of polygons as x,y coordinates from a Shapely geometry.
    We save x,y here and decide later in Napari/notebooks whether to swap to y,x.
    """
    if geom is None or geom.is_empty:
        return []

    if isinstance(geom, Polygon):
        geoms = [geom]

    elif isinstance(geom, MultiPolygon):
        geoms = list(geom.geoms)

    elif isinstance(geom, GeometryCollection):
        out = []
        for g in geom.geoms:
            out.extend(geometry_to_polygons_xy(g))
        return out

    else:
        return []

    polygons = []

    for poly in geoms:
        coords = np.asarray(poly.exterior.coords, dtype=np.float32)[:, :2]
        if coords.shape[0] >= 3:
            polygons.append(coords)

    return polygons


def main():
    parser = argparse.ArgumentParser(
        description="Export SOPA SpatialData cell boundary polygons to compressed NPZ."
    )

    parser.add_argument(
        "--sdata-path",
        required=True,
        help="Path to SOPA SpatialData .zarr directory.",
    )

    parser.add_argument(
        "--output-npz",
        required=True,
        help="Output .npz file.",
    )

    parser.add_argument(
        "--shape-key",
        default="",
        help="Optional shape layer name. If empty, choose automatically.",
    )

    args = parser.parse_args()

    sdata_path = Path(args.sdata_path)
    out_path = Path(args.output_npz)

    if not sdata_path.exists():
        raise FileNotFoundError(f"SpatialData path does not exist: {sdata_path}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("EXPORT SOPA CELL BOUNDARIES TO NPZ")
    print("=" * 80)
    print("SpatialData path:", sdata_path)
    print("Output NPZ:       ", out_path)
    print("Requested shape:  ", args.shape_key if args.shape_key else "auto")
    print("=" * 80)

    print("Reading SpatialData...")
    sdata = read_zarr(str(sdata_path))

    shape_key = pick_shape_key(sdata, preferred=args.shape_key)
    print("Using shape layer:", shape_key)

    gdf = sdata.shapes[shape_key]

    if hasattr(gdf, "compute"):
        gdf = gdf.compute()

    print("GeoDataFrame shape:", gdf.shape)
    print("Geometry types:")
    print(gdf.geometry.geom_type.value_counts())
    print("Bounds [minx, miny, maxx, maxy]:", gdf.total_bounds)

    polygons = []
    cell_ids = []

    for idx, row in gdf.iterrows():
        polys = geometry_to_polygons_xy(row.geometry)

        for p in polys:
            polygons.append(p)
            cell_ids.append(str(idx))

    print("Exported polygons:", len(polygons))

    if len(polygons) == 0:
        raise RuntimeError("No polygons were exported from the selected shape layer.")

    np.savez_compressed(
        out_path,
        polygons=np.asarray(polygons, dtype=object),
        cell_ids=np.asarray(cell_ids, dtype=str),
        shape_key=np.asarray(shape_key),
        bounds=np.asarray(gdf.total_bounds, dtype=float),
        source_sdata_path=np.asarray(str(sdata_path)),
    )

    print("Saved:", out_path)
    print("=" * 80)


if __name__ == "__main__":
    main()
