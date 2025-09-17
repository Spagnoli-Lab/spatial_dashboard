#!/usr/bin/env python3
"""Apply hand-drawn masks to DAPI images within spatial sample folders.

The script walks over sample directories, looks for matching DAPI and mask TIFF
files, applies the mask (values > 0 are kept), crops the result to the masked
region to remove black borders, and stores the masked output next to the
original DAPI file. Output files are named ``DAPI_<stage>_<sample>_masked.tif``
when the mask encodes a sample identifier and ``DAPI_<stage>_masked.tif``
otherwise.

Example
-------
Run over the default spatial test data::

    python utils/apply_masked_dapi.py

Run over a custom directory and overwrite previously generated outputs::

    python utils/apply_masked_dapi.py --root /path/to/spatial --overwrite
"""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np

try:
    import tifffile
except ImportError as exc:  # pragma: no cover - safeguard for missing dependency
    raise SystemExit(
        "The 'tifffile' package is required. Install it with `pip install tifffile`."
    ) from exc


MASK_HINT = re.compile(r"(e\d+(?:\.\d+)?)_(\d+)", re.IGNORECASE)
STAGE_HINT = re.compile(r"(e\d+(?:\.\d+)?)", re.IGNORECASE)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("test_data") / "spatial",
        help="Root directory that contains spatial sample folders (default: test_data/spatial)",
    )
    parser.add_argument(
        "paths",
        nargs="*",
        type=Path,
        help="Optional explicit sample directories to process; defaults to walking under --root",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing *_masked.tif outputs",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Inspect matches without writing any files",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )
    return parser.parse_args(argv)


def configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(message)s")


def discover_directories(root: Path, overrides: Sequence[Path]) -> List[Path]:
    candidates: List[Path] = []
    if overrides:
        for path in overrides:
            if not path.exists():
                logging.warning("Skipping missing path %s", path)
                continue
            if path.is_file():
                logging.warning("Skipping file %s (expected a directory)", path)
                continue
            candidates.append(path)
    else:
        if not root.exists():
            raise SystemExit(f"Root directory {root} does not exist")
        candidates.extend(sorted(p for p in root.rglob("*") if p.is_dir()))

    # Ensure deterministic order and deduplicate while keeping first encounter
    seen = set()
    unique_candidates: List[Path] = []
    for path in sorted(candidates):
        if path in seen:
            continue
        seen.add(path)
        unique_candidates.append(path)
    return unique_candidates


def extract_sample_token(path: Path) -> Optional[str]:
    match = MASK_HINT.search(path.stem)
    if match:
        return match.group(0).lower()
    match = STAGE_HINT.search(path.stem)
    if match:
        return match.group(1).lower()
    return None


def find_candidate_pairs(directory: Path) -> List[Tuple[Path, Path]]:
    tif_files = list(directory.glob("*.tif"))
    if not tif_files:
        return []

    dapi_files = [p for p in tif_files if "dapi" in p.name.lower()]
    mask_files = [p for p in tif_files if "8bit" in p.stem.lower() and p not in dapi_files]

    if not dapi_files or not mask_files:
        return []

    pairs: List[Tuple[Path, Path]] = []
    used_dapi: set[Path] = set()

    for mask_path in sorted(mask_files):
        mask_token = extract_sample_token(mask_path)
        matched_dapi: Optional[Path] = None

        if mask_token:
            for dapi_path in sorted(dapi_files):
                if dapi_path in used_dapi:
                    continue
                if mask_token in dapi_path.stem.lower():
                    matched_dapi = dapi_path
                    break

        if matched_dapi is None:
            remaining = [p for p in sorted(dapi_files) if p not in used_dapi]
            if len(remaining) == 1:
                matched_dapi = remaining[0]

        if matched_dapi is None:
            logging.warning(
                "Could not match mask %s to a DAPI image in %s", mask_path.name, directory
            )
            continue

        used_dapi.add(matched_dapi)
        pairs.append((matched_dapi, mask_path))

    if not pairs and len(dapi_files) == 1 and len(mask_files) == 1:
        pairs.append((dapi_files[0], mask_files[0]))

    return pairs


def load_image(path: Path) -> np.ndarray:
    logging.debug("Loading %s", path)
    return tifffile.imread(path)


def normalize_mask(mask: np.ndarray, target_shape: Tuple[int, ...]) -> np.ndarray:
    mask_bool = np.squeeze(mask > 0)
    if mask_bool.shape == target_shape:
        return mask_bool

    if mask_bool.shape == target_shape[:-1]:
        mask_bool = mask_bool[..., np.newaxis]

    try:
        mask_bool = np.broadcast_to(mask_bool, target_shape)
    except ValueError as exc:
        raise ValueError(
            f"Mask shape {mask_bool.shape} is incompatible with DAPI shape {target_shape}"
        ) from exc
    return mask_bool.astype(bool, copy=False)


def apply_mask(dapi: np.ndarray, mask: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    mask_bool = normalize_mask(mask, dapi.shape)
    masked = np.zeros_like(dapi)
    masked[mask_bool] = dapi[mask_bool]
    return masked, mask_bool


def crop_to_mask_region(image: np.ndarray, mask_bool: np.ndarray) -> np.ndarray:
    """Crop the image to the smallest bounding box covering the mask."""
    if not np.any(mask_bool):
        return image

    # Collapse any trailing dimensions so we crop only spatial axes
    if mask_bool.ndim > 2:
        spatial_axes = tuple(range(2, mask_bool.ndim))
        spatial_mask = np.any(mask_bool, axis=spatial_axes)
    else:
        spatial_mask = mask_bool

    coords = np.argwhere(spatial_mask)
    if coords.size == 0:
        return image

    mins = coords.min(axis=0)
    maxs = coords.max(axis=0) + 1

    slices = [slice(start, stop) for start, stop in zip(mins, maxs)]
    if image.ndim > spatial_mask.ndim:
        slices.extend([slice(None)] * (image.ndim - spatial_mask.ndim))

    cropped = image[tuple(slices)]
    return np.ascontiguousarray(cropped)


def derive_stage(source: Path) -> Optional[str]:
    match = STAGE_HINT.search(source.stem)
    if match:
        return match.group(1).upper()
    return None


def ensure_output_path(dapi_path: Path, mask_path: Path) -> Path:
    mask_match = MASK_HINT.search(mask_path.stem)
    if mask_match:
        stage = mask_match.group(1).upper()
        sample = mask_match.group(2)
        suffix = f"_{sample}"
    else:
        stage = derive_stage(mask_path) or derive_stage(dapi_path) or dapi_path.stem
        suffix = ""

    return dapi_path.with_name(f"DAPI_{stage}{suffix}_masked.tif")


def process_directory(directory: Path, overwrite: bool, dry_run: bool) -> None:
    pairs = find_candidate_pairs(directory)
    if not pairs:
        logging.debug("No matching DAPI/mask pair found in %s", directory)
        return

    for dapi_path, mask_path in pairs:
        output_path = ensure_output_path(dapi_path, mask_path)
        if output_path.exists() and not overwrite:
            logging.info("Skipping %s (exists); use --overwrite to regenerate", output_path)
            continue

        if dry_run:
            logging.info(
                "[DRY-RUN] Would create %s using %s and %s",
                output_path,
                dapi_path,
                mask_path,
            )
            continue

        dapi_image = load_image(dapi_path)
        mask_image = load_image(mask_path)

        try:
            masked, mask_bool = apply_mask(dapi_image, mask_image)
        except ValueError as exc:
            logging.warning("Skipping %s due to shape mismatch: %s", directory, exc)
            continue

        cropped = crop_to_mask_region(masked, mask_bool)

        logging.info("Writing %s", output_path)
        tifffile.imwrite(output_path, cropped)


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    configure_logging(args.verbose)

    directories = discover_directories(args.root, args.paths)
    if not directories:
        logging.info("No directories to process")
        return

    for directory in directories:
        process_directory(directory, overwrite=args.overwrite, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
