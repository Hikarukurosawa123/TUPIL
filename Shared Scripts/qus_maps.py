#!/usr/bin/env python3
"""Generate Gaussian ESD/EAC and MBF/SS/SI maps from Clarius RF data.

Python port of ``generate_mbf_ss_si_power_spectra.m``.  The numerical core
reproduces the active harmonic-1 estimator_gui path:

    RF load -> inverse TGC -> regional STFT attenuation correction
    -> phantom_NPS -> reference-phantom BSC -> spectral_analysis
    -> Gaussian form-factor fit -> overlapping-window parameter maps

The command-line defaults and QUS_* environment variables mirror the MATLAB
script.  A coordinator invocation creates a timestamped run directory and
launches one fresh Python process per 50 cases.  Worker output is HDF5 rather
than MATLAB v7.3: ``all_parametric_maps`` has logical dimension order
axial x lateral x parameter x case, and the parameter names are recorded in
``map_order``.  Case-level settings, ROI positions, vectors, axes, status,
and source paths are stored below each file's ``cases/<zero-padded index>``
group.  CSV manifests have the same principal columns as the MATLAB output.

Examples
--------
Run the biopsy cohort with defaults::

    python3 qus_maps.py

Run a small mapped cohort or one RF file::

    QUS_MAPPING_PATH=subset.csv QUS_COHORT_MODE=all_raw \
      python3 qus_maps.py --case-rf-file P001_..._rf.raw

Exercise the estimator core without project data::

    python3 qus_maps.py --self-test
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import struct
import subprocess
import sys
import tempfile
import time
import traceback
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.io import loadmat


SCRIPT_PATH = Path(__file__).resolve()
ANALYSIS_DIR = SCRIPT_PATH.parent
REPO_ROOT = ANALYSIS_DIR.parent

MAP_ORDER = (
    "Gaussian_ESD_um",
    "Gaussian_EAC_dB_per_cm",
    "MBF_dB",
    "SS_dB_per_MHz",
    "SI_dB",
)
MAP_DIMENSION_ORDER = ("axial", "lateral", "parameter", "case")
MANIFEST_COLUMNS = (
    "PatientID", "VisitType", "VisitNumber", "ImageNumber", "CICt",
    "CaseName", "RFFile", "Phantom", "TGCRemoval",
    "AttenuationCorrection", "SampleAttenuationMode",
    "SampleAttenuation_dB_per_cm_MHz",
    "NonKidneyAttenuation_dB_per_cm_MHz",
    "KidneyAttenuation_dB_per_cm_MHz",
    "PhantomAttenuation_dB_per_cm_MHz", "RFFrame", "ROIFile",
    "NumROIs", "F0_MHz", "F1_MHz", "MBF_dB", "SS_dB_per_MHz",
    "SI_dB", "MeanR2", "GaussianESD_um",
    "GaussianEAC_dB_per_cm", "GaussianMeanR2", "GaussianMeanChi2",
    "SampleID", "MapRows", "MapColumns", "MapGroup", "GroupCaseIndex",
    "MapHDF5File", "Status", "ErrorMessage",
)


def env_path(name: str, default: Path) -> Path:
    value = os.environ.get(name, "").strip()
    return Path(value).expanduser().resolve() if value else default


@dataclass(frozen=True)
class Config:
    biopsy_path: Path
    mapping_path: Path
    raw_root: Path
    phantom_root: Path
    reference_bsc_path: Path
    cohort_mode: str = "biopsy"
    output_folder_name: str = (
        "QUS_TGCremoved_STFTAtten_kidney_0p50_non_kidney_0p80_"
        "6_by_6_ROI_minus18dB_python"
    )
    apply_tgc_removal: bool = True
    apply_attenuation_correction: bool = True
    sample_attenuation_mode: str = "regional"
    non_kidney_alpha: float = 0.80
    kidney_alpha: float = 0.50
    sample_alpha: float = 0.80
    phantom_alpha: float = 0.439
    stft_taper: bool = True
    nfft: int = 4096
    cutoff_db: float = -18.0
    sound_speed: float = 1540.0
    roi_size_axial_mm: float = 6.0
    roi_size_lateral_mm: float = 6.0
    roi_overlap_axial: float = 0.9
    roi_overlap_lateral: float = 0.9
    water_temperature_c: float = 37.0
    roi_dimension_suffix: str = ""
    chunk_size: int = 50
    max_folder_bytes: int = 3_000_000_000


@dataclass(frozen=True)
class RawHeader:
    file_id: int
    frames: int
    lines: int
    samples: int
    sample_size: int


def find_reference_bsc() -> Path:
    candidates = (
        ANALYSIS_DIR / "SMH_LF_BSC.txt",
        REPO_ROOT / "ignored_files" / "#BSC GUI" / "bsc_gui_TB"
        / "BSC_GUI_v109 - ORIGINAL" / "SMH_LF_BSC.txt",
        REPO_ROOT / "ignored_files" / "#BSC GUI" / "BSC_GUI_Clarius"
        / "BSC_GUI_v109 - ORIGINAL" / "SMH_LF_BSC.txt",
        REPO_ROOT / "ignored_files" / "#BSC GUI"
        / "BSC_GUI_v109 - ORIGINAL" / "SMH_LF_BSC.txt",
    )
    gui_override = os.environ.get("QUS_GUI_DIR", "").strip()
    if gui_override:
        candidates = (Path(gui_override).expanduser() / "SMH_LF_BSC.txt",)
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    return candidates[0].resolve()


def default_config() -> Config:
    cohort_mode = os.environ.get("QUS_COHORT_MODE", "biopsy").strip().lower()
    if cohort_mode not in {"biopsy", "all_raw"}:
        raise ValueError("QUS_COHORT_MODE must be biopsy or all_raw")
    if cohort_mode == "all_raw":
        default_output = (
            "QUS_TGCremoved_STFTAtten_kidney_0p50_non_kidney_0p80_"
            "6_by_6_ROI_all_raw_cases_minus18dB_python"
        )
    else:
        default_output = (
            "QUS_TGCremoved_STFTAtten_kidney_0p50_non_kidney_0p80_"
            "6_by_6_ROI_minus18dB_python"
        )
    return Config(
        biopsy_path=REPO_ROOT / "csv"
        / "20260506_hscanV3_image_biopsy_2p5MHz_upToP655_attenuation.xlsx",
        mapping_path=env_path(
            "QUS_MAPPING_PATH",
            REPO_ROOT / "Data" / "Raw_data_cleaned"
            / "patient_rf_with_phantoms.csv",
        ),
        raw_root=env_path("QUS_RAW_ROOT", REPO_ROOT / "Data" / "Raw_data_cleaned"),
        phantom_root=env_path(
            "QUS_PHANTOM_ROOT", REPO_ROOT / "ignored_files" / "#phantom_new"
        ),
        reference_bsc_path=find_reference_bsc(),
        cohort_mode=cohort_mode,
        output_folder_name=os.environ.get(
            "QUS_OUTPUT_FOLDER_NAME", default_output
        ).strip(),
    )


def require_inputs(config: Config) -> None:
    if config.cohort_mode == "biopsy" and not config.biopsy_path.is_file():
        raise FileNotFoundError(f"Biopsy workbook not found: {config.biopsy_path}")
    for label, path, is_dir in (
        ("mapping", config.mapping_path, False),
        ("raw-data root", config.raw_root, True),
        ("phantom root", config.phantom_root, True),
        ("reference BSC", config.reference_bsc_path, False),
    ):
        valid = path.is_dir() if is_dir else path.is_file()
        if not valid:
            raise FileNotFoundError(f"{label} not found: {path}")
    if not config.output_folder_name:
        raise ValueError("Output folder name must not be empty")


def read_raw_header(path: Path) -> RawHeader:
    with path.open("rb") as stream:
        raw = stream.read(20)
    if len(raw) != 20:
        raise ValueError(f"Invalid Clarius header: {path}")
    return RawHeader(*struct.unpack("<5i", raw))


def read_rf_frame(path: Path, frame_number: int) -> tuple[np.ndarray, int, RawHeader]:
    """Read one 1-based Clarius RF frame, matching rdataread's layout."""
    header = read_raw_header(path)
    if header.file_id != 2:
        raise ValueError(f"Not a Clarius RF file (id={header.file_id}): {path}")
    if not 1 <= frame_number <= header.frames:
        raise ValueError(
            f"Frame {frame_number} is outside 1-{header.frames}: {path}"
        )
    if header.sample_size != 2:
        raise ValueError(f"Expected 2-byte RF samples, got {header.sample_size}: {path}")
    frame_bytes = 8 + header.samples * header.lines * header.sample_size
    with path.open("rb") as stream:
        stream.seek(20 + (frame_number - 1) * frame_bytes)
        timestamp_raw = stream.read(8)
        if len(timestamp_raw) != 8:
            raise EOFError(f"Missing timestamp for frame {frame_number}: {path}")
        timestamp = struct.unpack("<q", timestamp_raw)[0]
        values = np.fromfile(
            stream, dtype="<i2", count=header.samples * header.lines
        )
    if values.size != header.samples * header.lines:
        raise EOFError(f"Truncated frame {frame_number}: {path}")
    return (
        values.reshape((header.samples, header.lines), order="F").astype(np.float64),
        timestamp,
        header,
    )


def yaml_number(text: str, key: str) -> float:
    match = re.search(
        rf"^\s*{re.escape(key)}:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)",
        text,
        flags=re.MULTILINE,
    )
    if not match:
        raise ValueError(f"YAML field {key!r} was not found")
    return float(match.group(1))


def read_clarius_metadata(raw_path: Path, expected_lines: int) -> dict[str, float]:
    yml_path = raw_path.with_suffix(".yml")
    text = yml_path.read_text(encoding="utf-8", errors="replace")
    rx = [float(value) for value in re.findall(r"rx element:\s*([-+]?\d*\.?\d+)", text)]
    if len(rx) < expected_lines:
        raise ValueError(
            f"Expected {expected_lines} rx-element entries, found {len(rx)}: {yml_path}"
        )
    return {
        "probe_elements": yaml_number(text, "elements"),
        "probe_pitch_um": yaml_number(text, "pitch"),
        "frames": yaml_number(text, "frames"),
        "transmit_frequency_mhz": yaml_number(text, "transmit frequency"),
        "imaging_depth_mm": yaml_number(text, "imaging depth"),
        "focal_depth_mm": yaml_number(text, "focal depth"),
        "sampling_rate_mhz": yaml_number(text, "sampling rate"),
        "delay_samples": yaml_number(text, "delay samples"),
        "first_rx_element": rx[0],
        "last_rx_element": rx[expected_lines - 1],
    }


def load_rf_data(path: Path, frame_number: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, RawHeader, dict[str, float], int]:
    """Port of the active Clarius branch of load_RF_data.m."""
    rf, timestamp, header = read_rf_frame(path, frame_number)
    metadata = read_clarius_metadata(path, header.lines)
    rate_hz = metadata["sampling_rate_mhz"] * 1e6
    delay_mm = (1_538_500.0 / 2.0) * metadata["delay_samples"] / rate_hz
    depth_mm = (1_538_500.0 / 2.0) * header.samples / rate_hz
    axial = np.linspace(delay_mm, depth_mm + delay_mm, header.samples)
    # The MATLAB script replaces load_RF_data's geometric width with these
    # absolute receive-element coordinates immediately after loading.
    lateral = np.linspace(
        metadata["first_rx_element"], metadata["last_rx_element"], header.lines
    )
    return rf, axial, lateral, 15.0, header, metadata, timestamp


def tgc_path_for(raw_path: Path) -> Path:
    return raw_path.with_name(re.sub(r"_rf\.raw$", "_env.tgc.yml", raw_path.name, flags=re.I))


def apply_clarius_tgc(
    rf: np.ndarray, axial_mm: np.ndarray, tgc_path: Path, target_timestamp: int | None
) -> np.ndarray:
    if not tgc_path.is_file():
        print(f"Warning: TGC file not found; RF unchanged: {tgc_path}", file=sys.stderr)
        return rf
    text = tgc_path.read_text(encoding="utf-8", errors="replace")
    blocks = re.findall(
        r"timestamp:\s*(\d+)\s*(.*?)(?=\r?\ntimestamp:|\Z)", text, flags=re.S
    )
    active = text
    if blocks:
        chosen = 0
        if target_timestamp is not None:
            eligible = [i for i, (stamp, _) in enumerate(blocks) if int(stamp) <= target_timestamp]
            if eligible:
                chosen = eligible[-1]
        active = blocks[chosen][1]
    points = re.findall(
        r"\{\s*([-+]?\d*\.?\d+)mm,\s*([-+]?\d*\.?\d+)dB\s*\}", active
    )
    if not points:
        print(f"Warning: no TGC points; RF unchanged: {tgc_path}", file=sys.stderr)
        return rf
    depth = np.asarray([float(p[0]) for p in points])
    gain = np.asarray([float(p[1]) for p in points])
    _, first = np.unique(depth, return_index=True)
    order = first[np.argsort(depth[first])]
    depth, gain = depth[order], gain[order]
    query = np.clip(np.asarray(axial_mm, dtype=float), depth[0], depth[-1])
    gain_db = np.full(query.shape, gain[0]) if len(depth) == 1 else np.interp(query, depth, gain)
    return np.asarray(rf, dtype=float) / (10.0 ** (gain_db[:, None] / 20.0))


def stft_starts(signal_length: int, window_length: int = 256, hop: int = 64) -> tuple[np.ndarray, int]:
    pad = window_length // 2
    padded_length = signal_length + 2 * pad
    last = padded_length - window_length
    starts = np.arange(0, last + 1, hop, dtype=int)
    if starts[-1] != last:
        starts = np.append(starts, last)
    return starts, pad


def apply_stft_attenuation(
    rf: np.ndarray,
    sampling_frequency_mhz: float,
    axial_mm: np.ndarray,
    kidney_mask: np.ndarray,
    non_kidney_alpha: float,
    kidney_alpha: float,
    taper: bool = True,
) -> np.ndarray:
    """Vectorized port of apply_estimator_stft_attenuation_by_region.m (taper=False drops the 2.5 MHz Gaussian gain weight)."""
    rf = np.where(np.isfinite(rf), rf, 0.0).astype(np.float64, copy=False)
    kidney_mask = np.asarray(kidney_mask, dtype=bool)
    if kidney_mask.shape != rf.shape:
        raise ValueError("Kidney mask and RF frame sizes must match")
    if len(axial_mm) != rf.shape[0]:
        raise ValueError("Axial coordinates must match RF rows")
    nperseg, hop, nfft = 256, 64, 1024
    window = 0.5 - 0.5 * np.cos(2.0 * np.pi * np.arange(nperseg) / nperseg)
    starts, pad = stft_starts(rf.shape[0], nperseg, hop)
    padded = np.pad(rf, ((pad, pad), (0, 0)))
    indices = starts[:, None] + np.arange(nperseg)[None, :]
    spectra = np.fft.rfft(
        padded[indices, :] * window[None, :, None], n=nfft, axis=1
    )
    fs_hz = sampling_frequency_mhz * 1e6 if sampling_frequency_mhz < 1e5 else sampling_frequency_mhz
    frequency_hz = np.arange(nfft // 2 + 1) * fs_hz / nfft
    frequency_mhz = frequency_hz / 1e6
    sigma_hz = ((1.0 + 0.6) * 2.5e6 - (1.0 - 0.6) * 2.5e6) / 3.0
    frequency_weight = np.exp(-0.5 * ((frequency_hz - 2.5e6) / sigma_hz) ** 2)

    positions_1based = np.clip(
        starts - pad + (nperseg + 1) / 2.0, 1.0, float(rf.shape[0])
    )
    path_cm = np.r_[0.0, np.cumsum(np.abs(np.diff(axial_mm)) / 10.0)]
    alpha = np.where(kidney_mask, kidney_alpha, non_kidney_alpha)
    cumulative = cumulative_trapezoid(alpha, path_cm, axis=0, initial=0.0)
    lo = np.floor(positions_1based).astype(int)
    hi = np.minimum(lo + 1, rf.shape[0])
    fraction = positions_1based - lo
    path_attenuation = (
        (1.0 - fraction)[:, None] * cumulative[lo - 1, :]
        + fraction[:, None] * cumulative[hi - 1, :]
    )
    attenuation_db = 2.0 * frequency_mhz[None, :, None] * path_attenuation[:, None, :]
    gain = np.minimum(10.0 ** (attenuation_db / 20.0), 100.0)
    if taper:
        gain *= frequency_weight[None, :, None]

    frames = np.fft.irfft(spectra * gain, n=nfft, axis=1)[:, :nperseg, :]
    frames *= window[None, :, None]
    reconstructed = np.zeros_like(padded)
    normalization = np.zeros(padded.shape[0])
    for frame_index, start in enumerate(starts):
        reconstructed[start : start + nperseg, :] += frames[frame_index]
        normalization[start : start + nperseg] += window**2
    valid = normalization > np.finfo(float).eps * normalization.max()
    reconstructed[valid, :] /= normalization[valid, None]
    return reconstructed[pad : pad + rf.shape[0], :]


ROI_FIELDS = ("top", "bottom", "left", "right")


def matlab_struct_to_roi(value: Any) -> dict[str, np.ndarray]:
    result: dict[str, np.ndarray] = {}
    for field in ROI_FIELDS:
        if hasattr(value, field):
            raw = getattr(value, field)
        elif isinstance(value, np.ndarray) and value.dtype.names and field in value.dtype.names:
            raw = value[field].item()
        else:
            raise ValueError(f"ROI_positions is missing field {field}")
        result[field] = np.atleast_1d(np.asarray(raw, dtype=float).squeeze())
    counts = {len(v) for v in result.values()}
    if len(counts) != 1 or not counts or next(iter(counts)) == 0:
        raise ValueError("ROI boundary fields are empty or have unequal lengths")
    if not all(np.isfinite(v).all() for v in result.values()):
        raise ValueError("ROI boundary fields contain nonfinite values")
    if np.any(result["top"] > result["bottom"]) or np.any(result["left"] > result["right"]):
        raise ValueError("ROI boundaries are reversed")
    return result


def load_roi(path: Path) -> dict[str, np.ndarray]:
    data = loadmat(path, variable_names=["ROI_positions"], squeeze_me=True, struct_as_record=False)
    if "ROI_positions" not in data:
        raise ValueError(f"ROI_positions was not found: {path}")
    return matlab_struct_to_roi(data["ROI_positions"])


def roi_pattern(sample_base: str, suffix: str = "") -> re.Pattern[str]:
    tail = rf"_{re.escape(suffix)}" if suffix else ""
    return re.compile(rf"^{re.escape(sample_base)}_ROI_(\d+){tail}\.mat$")


def choose_standard_roi(roi_dir: Path, sample_base: str, raw_path: Path, suffix: str = "") -> tuple[Path, int]:
    header = read_raw_header(raw_path)
    pattern = roi_pattern(sample_base, suffix)
    candidates: list[tuple[int, Path]] = []
    if roi_dir.is_dir():
        for path in roi_dir.iterdir():
            match = pattern.match(path.name)
            if match and 1 <= int(match.group(1)) <= header.frames:
                candidates.append((int(match.group(1)), path))
    if not candidates:
        raise FileNotFoundError(f"No standard ROI with a valid RF frame in {roi_dir}")
    frame, path = max(candidates, key=lambda item: item[0])
    return path, frame


def has_valid_standard_roi(roi_dir: Path, sample_base: str, raw_path: Path, suffix: str = "") -> bool:
    try:
        path, _ = choose_standard_roi(roi_dir, sample_base, raw_path, suffix)
        load_roi(path)
        return True
    except Exception:
        return False


def has_standard_roi(roi_dir: Path, sample_base: str, suffix: str = "") -> bool:
    """MATLAB select_all_cases eligibility check (name/shape not loaded yet)."""
    pattern = roi_pattern(sample_base, suffix)
    return roi_dir.is_dir() and any(
        path.is_file() and pattern.match(path.name) for path in roi_dir.iterdir()
    )


def kidney_mask_path(roi_dir: Path, sample_base: str) -> Path | None:
    expected = roi_dir / f"{sample_base}_raw_kidney.mat"
    if expected.is_file():
        return expected
    candidates = sorted(roi_dir.glob("*rf_raw_kidney.mat")) if roi_dir.is_dir() else []
    return candidates[0] if len(candidates) == 1 else None


def load_kidney_mask(roi_dir: Path, sample_base: str, shape: tuple[int, int]) -> np.ndarray:
    path = kidney_mask_path(roi_dir, sample_base)
    if path is None:
        raise FileNotFoundError(f"Expected one kidney-mask file in {roi_dir}")
    data = loadmat(path, variable_names=["kidney_mask"])
    if "kidney_mask" not in data:
        raise ValueError(f"kidney_mask was not found: {path}")
    mask = np.asarray(data["kidney_mask"], dtype=bool)
    if mask.shape == shape[::-1]:
        mask = mask.T
    if mask.shape != shape:
        # Nearest-neighbour resampling on normalized coordinates.
        row = np.rint(np.linspace(0, mask.shape[0] - 1, shape[0])).astype(int)
        col = np.rint(np.linspace(0, mask.shape[1] - 1, shape[1])).astype(int)
        mask = mask[np.ix_(row, col)]
    if not mask.any():
        raise ValueError(f"Kidney mask is empty: {path}")
    return mask


def transform_roi_lateral(roi: dict[str, np.ndarray], old_axis: np.ndarray, new_axis: np.ndarray) -> dict[str, np.ndarray]:
    # In this Python port load_rf_data already returns the new receive-element
    # axis.  This helper remains explicit for saved ROI files whose coordinates
    # were made on the GUI geometric axis.
    result = {key: value.copy() for key, value in roi.items()}
    if np.array_equal(old_axis, new_axis):
        return result
    transform = interp1d(old_axis, new_axis, kind="linear", fill_value="extrapolate")
    result["left"] = transform(result["left"])
    result["right"] = transform(result["right"])
    swapped = result["left"] > result["right"]
    result["left"][swapped], result["right"][swapped] = (
        result["right"][swapped].copy(), result["left"][swapped].copy()
    )
    return result


def gui_geometric_lateral_axis(metadata: dict[str, float], nlines: int, nsamples: int) -> np.ndarray:
    rate_hz = metadata["sampling_rate_mhz"] * 1e6
    delay_mm = (1_538_500.0 / 2.0) * metadata["delay_samples"] / rate_hz
    depth_mm = (1_538_500.0 / 2.0) * nsamples / rate_hz
    arc_length = metadata["probe_pitch_um"] / 1e3 * nlines
    fov_degrees = arc_length * 360.0 / (2.0 * np.pi * 45.0)
    width = 2.0 * (depth_mm + delay_mm) * np.tan(np.deg2rad(fov_degrees / 2.0))
    return np.linspace(0.0, width, nlines)


def keep_rois_inside_grid(roi: dict[str, np.ndarray], axial: np.ndarray, lateral: np.ndarray) -> dict[str, np.ndarray]:
    valid = (
        (roi["top"] >= np.min(axial))
        & (roi["bottom"] <= np.max(axial))
        & (roi["left"] >= np.min(lateral))
        & (roi["right"] <= np.max(lateral))
    )
    return {field: roi[field][valid] for field in ROI_FIELDS}


def nearest_indices(axis: np.ndarray, values: np.ndarray) -> np.ndarray:
    return np.abs(np.asarray(axis)[:, None] - np.asarray(values)[None, :]).argmin(axis=0)


def roi_indices(roi: dict[str, np.ndarray], axial: np.ndarray, lateral: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "top": nearest_indices(axial, roi["top"]),
        "bottom": nearest_indices(axial, roi["bottom"]),
        "left": nearest_indices(lateral, roi["left"]),
        "right": nearest_indices(lateral, roi["right"]),
    }


def match_sample_to_reference_grid(
    sample: np.ndarray,
    sample_axial: np.ndarray,
    sample_lateral: np.ndarray,
    reference_shape: tuple[int, int],
    reference_axial: np.ndarray,
    reference_lateral: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if (
        sample.shape == reference_shape
        and np.array_equal(sample_axial, reference_axial)
        and np.array_equal(sample_lateral, reference_lateral)
    ):
        return sample, reference_axial, reference_lateral
    interpolator = RegularGridInterpolator(
        (sample_axial, sample_lateral), sample, method="linear",
        bounds_error=False, fill_value=0.0,
    )
    query_axial, query_lateral = np.meshgrid(
        reference_axial, reference_lateral, indexing="ij"
    )
    points = np.column_stack((query_axial.ravel(), query_lateral.ravel()))
    matched = interpolator(points).reshape(reference_shape)
    return matched, reference_axial, reference_lateral


def phantom_nps(
    rf: np.ndarray,
    axial: np.ndarray,
    lateral: np.ndarray,
    fs_mhz: float,
    roi: dict[str, np.ndarray],
    nfft: int,
    gain: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Port of phantom_NPS.m; columns are ROIs and rows are FFT bins."""
    indices = roi_indices(roi, axial, lateral)
    spectra = np.empty((nfft, len(indices["left"])), dtype=float)
    for n, (top, bottom, left, right) in enumerate(
        zip(indices["top"], indices["bottom"], indices["left"], indices["right"])
    ):
        block = rf[top : bottom + 1, left : right + 1]
        window = np.hanning(block.shape[0])[:, None]
        power = np.mean(np.abs(np.fft.fft(block * window, n=nfft, axis=0)) ** 2, axis=1)
        if gain:
            power += np.sign(gain) * 10.0 ** (7.4425 + abs(gain) * 0.0025)
        spectra[:, n] = power
    frequency = np.linspace(0.0, fs_mhz, nfft)
    return frequency, spectra


def estimator_bandwidth(
    phantom_rf: np.ndarray,
    sample_axial: np.ndarray,
    sample_lateral: np.ndarray,
    frequency_mhz: np.ndarray,
    focal_depth_mm: float,
    wavelength_mm: float,
    nfft: int,
    cutoff_db: float,
) -> tuple[int, int]:
    up = int(np.argmin(np.abs(sample_axial - (focal_depth_mm - wavelength_mm * 10.0))))
    down = int(np.argmin(np.abs(sample_axial - (focal_depth_mm + wavelength_mm * 10.0))))
    first_row, last_row = sorted((up, down))
    window = np.hanning(last_row - first_row + 1)[:, None]
    # MATLAB round is half-away-from-zero. These quantities are positive.
    first_column_1 = int(math.floor(len(sample_lateral) * 0.1 + 0.5))
    last_column_1 = int(math.floor(len(sample_lateral) * 0.9 + 0.5))
    focus = phantom_rf[
        first_row : last_row + 1, first_column_1 - 1 : last_column_1
    ]
    spectrum = np.mean(np.abs(np.fft.fft(focus * window, n=nfft, axis=0)) ** 2, axis=1)
    search = np.arange(100, nfft // 2)  # MATLAB 101:N_fft/2
    peak = search[np.argmax(spectrum[search])]
    with np.errstate(divide="ignore"):
        keep = 10.0 * np.log10(spectrum[search]) > 10.0 * np.log10(spectrum[peak]) + cutoff_db
    bandwidth = search[keep]
    if not bandwidth.size:
        raise ValueError("Estimator bandwidth search returned no frequency bins")
    f0, f1 = frequency_mhz[bandwidth[0]], frequency_mhz[bandwidth[-1]]
    return int(np.searchsorted(frequency_mhz, f0, side="left")), int(np.searchsorted(frequency_mhz, f1, side="right") - 1)


def reference_phantom_bsc(
    sample_nps: np.ndarray,
    phantom_nps_values: np.ndarray,
    reference_bsc: np.ndarray,
    band: slice,
) -> tuple[np.ndarray, np.ndarray]:
    if sample_nps.shape != phantom_nps_values.shape:
        raise ValueError("Sample and phantom NPS arrays do not match")
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = sample_nps / phantom_nps_values
    return ratio, ratio[band, :] * reference_bsc[:, None]


def spectral_analysis(
    normalized_ratio: np.ndarray, frequency_mhz: np.ndarray, f0: int, f1: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    band_frequency = frequency_mhz[f0 : f1 + 1]
    with np.errstate(divide="ignore", invalid="ignore"):
        db = 10.0 * np.log10(normalized_ratio[f0 : f1 + 1, :])
    count = len(band_frequency)
    matlab_mid_1 = int(math.floor(count / 2.0 + 0.5))
    midpoint_frequency = band_frequency[matlab_mid_1 - 1]
    coefficients = np.full((db.shape[1], 2), np.nan)
    mbf = np.full(db.shape[1], np.nan)
    r2 = np.full(db.shape[1], np.nan)
    for roi_index in range(db.shape[1]):
        values = db[:, roi_index]
        finite = np.isfinite(values) & np.isfinite(band_frequency)
        if finite.sum() < 2:
            continue
        slope, intercept = np.polyfit(band_frequency[finite], values[finite], 1)
        fit = slope * band_frequency[finite] + intercept
        serr = np.sum((fit - values[finite]) ** 2)
        stot = np.sum((values[finite] - np.mean(values[finite])) ** 2)
        coefficients[roi_index] = (slope, intercept)
        mbf[roi_index] = slope * midpoint_frequency + intercept
        r2[roi_index] = 1.0 - serr / stot if stot != 0 else np.nan
    return mbf, coefficients, r2, db


def gaussian_form_factor_fit(
    bsc: np.ndarray,
    frequency_mhz: np.ndarray,
    sound_speed_m_s: float,
    diameters_um: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Port of estimate_FF_simple + fit_ESD_EAC for Gaussian FF."""
    diameters = np.arange(1.0, 201.0) if diameters_um is None else np.asarray(diameters_um, dtype=float)
    ka = np.pi * frequency_mhz[:, None] * diameters[None, :] / sound_speed_m_s
    theory = frequency_mhz[:, None] ** 4 * np.exp(-0.827 * ka**2)
    n_rois = bsc.shape[1]
    esd = np.full(n_rois, np.nan)
    eac = np.full(n_rois, np.nan)
    r2 = np.full(n_rois, np.nan)
    chi2 = np.full(n_rois, np.nan)
    c_constant = np.pi**4 / (36.0 * (sound_speed_m_s * 100.0 / 1e6) ** 4)
    for roi_index in range(n_rois):
        measured = bsc[:, roi_index]
        valid = np.isfinite(measured) & (measured > 0)
        if valid.sum() < 2:
            continue
        with np.errstate(divide="ignore", invalid="ignore"):
            residual_db = 10.0 * np.log10(measured[valid, None] / theory[valid, :])
        means = np.mean(residual_db, axis=0)
        mse = np.sum((residual_db - means[None, :]) ** 2, axis=0)
        best = int(np.nanargmin(mse))
        diameter = diameters[best]
        scale_denominator = theory[valid, best] * c_constant * (diameter * 1e-4) ** 6
        concentration = 10.0 * np.log10(np.mean(measured[valid] / scale_denominator))
        fitted = theory[valid, best] * c_constant * 10.0 ** (concentration / 10.0) * (diameter * 1e-4) ** 6
        ss_error = np.sum((fitted - measured[valid]) ** 2)
        ss_mean = np.sum((measured[valid] - np.mean(measured[valid])) ** 2)
        esd[roi_index], eac[roi_index] = diameter, concentration
        r2[roi_index] = 1.0 - ss_error / ss_mean if ss_mean else np.nan
        chi2[roi_index] = ss_error
    return esd, eac, r2, chi2


def create_estimator_map(
    shape: tuple[int, int], axial: np.ndarray, lateral: np.ndarray,
    roi: dict[str, np.ndarray], estimates: np.ndarray,
) -> np.ndarray:
    indices = roi_indices(roi, axial, lateral)
    total = np.zeros(shape, dtype=float)
    count = np.zeros(shape, dtype=np.uint16)
    for top, bottom, left, right, estimate in zip(
        indices["top"], indices["bottom"], indices["left"], indices["right"], estimates
    ):
        if np.isnan(estimate):
            continue
        area = np.s_[top : bottom + 1, left : right + 1]
        total[area] += estimate
        count[area] += 1
    with np.errstate(divide="ignore", invalid="ignore"):
        output = total / count
    output[(~np.isfinite(output)) & (count > 0)] = 0.0
    return output


def normalize_ci_ct(values: pd.Series) -> pd.Series:
    labels = values.astype("string").str.strip().str.replace(" ", "", regex=False)
    return labels.where(labels.isin(["0", "1", "1+", "2", "2+", "3"]), "").fillna("")


def parse_image_number(value: str) -> float:
    match = re.search(r"Image_?(\d+)", str(value), flags=re.I)
    return float(match.group(1)) if match else np.nan


def parse_case_metadata(case_name: str) -> tuple[float, str, float]:
    match = re.match(r"^P0*(\d+)_([^_]+)_0*(\d+)_", str(case_name))
    if not match:
        return np.nan, "", np.nan
    return float(match.group(1)), match.group(2).strip().upper(), float(match.group(3))


def clean_mapping(path: Path) -> pd.DataFrame:
    mapping = pd.read_csv(path, dtype="string", keep_default_na=True)
    required = ["Case_Name", "Image_Folder", "RF_File", "AssignedPhantom"]
    missing = [name for name in required if name not in mapping.columns]
    if missing:
        raise ValueError(f"Mapping CSV is missing columns: {', '.join(missing)}")
    for source, target in zip(required, ["CaseName", "ImageFolder", "RFFile", "Phantom"]):
        mapping[target] = mapping[source].astype("string").str.strip()
    mapping["ImageNumeric"] = mapping["ImageFolder"].map(parse_image_number)
    return mapping


def mapped_row_eligible(row: pd.Series, config: Config, require_mask: bool, validate_roi: bool) -> bool:
    if any(pd.isna(row[name]) or not str(row[name]).strip() for name in ("CaseName", "ImageFolder", "RFFile", "Phantom")):
        return False
    image_dir = config.raw_root / str(row.CaseName) / str(row.ImageFolder)
    raw = image_dir / str(row.RFFile)
    phantom = config.phantom_root / f"{row.Phantom}_rf.raw"
    base = re.sub(r"\.raw$", "", str(row.RFFile), flags=re.I)
    roi_dir = image_dir / "ROIs"
    if not raw.is_file() or not phantom.is_file():
        return False
    valid = (
        has_valid_standard_roi(
            roi_dir, base, raw, config.roi_dimension_suffix
        )
        if validate_roi
        else has_standard_roi(roi_dir, base, config.roi_dimension_suffix)
    )
    if not valid:
        return False
    return not require_mask or kidney_mask_path(roi_dir, base) is not None


def select_cases(config: Config) -> pd.DataFrame:
    mapping = clean_mapping(config.mapping_path)
    columns = ["PatientID", "VisitType", "VisitNumber", "ImageNumber", "CICtLabel", "CaseName", "ImageFolder", "RFFile", "Phantom"]
    selected: list[dict[str, Any]] = []
    require_mask = config.apply_attenuation_correction and config.sample_attenuation_mode == "regional"
    if config.cohort_mode == "all_raw":
        mapping = mapping.sort_values(["CaseName", "ImageNumeric", "ImageFolder", "RFFile"], na_position="last")
        for _, row in mapping.iterrows():
            patient, visit_type, visit = parse_case_metadata(str(row.CaseName))
            if not (np.isfinite(patient) and visit_type and np.isfinite(visit) and np.isfinite(row.ImageNumeric)):
                continue
            if not mapped_row_eligible(row, config, require_mask, validate_roi=True):
                continue
            selected.append(dict(zip(columns, [patient, visit_type, visit, row.ImageNumeric, "not_filtered", row.CaseName, row.ImageFolder, row.RFFile, row.Phantom])))
        result = pd.DataFrame(selected, columns=columns).drop_duplicates(["CaseName", "ImageFolder", "RFFile"], keep="first")
    else:
        biopsy = pd.read_excel(config.biopsy_path)
        required = ["Patient ID", "Visit Type", "Visit Number", "Image Number", "ci/ct"]
        missing = [name for name in required if name not in biopsy.columns]
        if missing:
            raise ValueError(f"Biopsy workbook is missing columns: {', '.join(missing)}")
        biopsy["Patient ID"] = biopsy["Patient ID"].ffill()
        biopsy["CICtLabel"] = normalize_ci_ct(biopsy["ci/ct"])
        biopsy["PatientNumeric"] = pd.to_numeric(biopsy["Patient ID"], errors="coerce")
        biopsy["VisitNumeric"] = pd.to_numeric(biopsy["Visit Number"], errors="coerce")
        biopsy["ImageNumeric"] = pd.to_numeric(biopsy["Image Number"], errors="coerce")
        biopsy["VisitTypeClean"] = biopsy["Visit Type"].astype("string").str.strip().str.upper()
        valid = biopsy[["PatientNumeric", "VisitNumeric", "ImageNumeric"]].notna().all(axis=1) & biopsy.CICtLabel.ne("")
        biopsy = biopsy[valid].sort_values(["PatientNumeric", "VisitTypeClean", "VisitNumeric", "ImageNumeric"])
        for _, biopsy_row in biopsy.iterrows():
            pattern = re.compile(rf"^P0*{int(biopsy_row.PatientNumeric)}_{re.escape(str(biopsy_row.VisitTypeClean))}_0*{int(biopsy_row.VisitNumeric)}_", re.I)
            matches = mapping[mapping.CaseName.fillna("").map(lambda value: bool(pattern.match(str(value)))) & mapping.ImageNumeric.eq(biopsy_row.ImageNumeric)]
            for _, row in matches.iterrows():
                if mapped_row_eligible(row, config, False, validate_roi=False):
                    selected.append(dict(zip(columns, [biopsy_row.PatientNumeric, biopsy_row.VisitTypeClean, biopsy_row.VisitNumeric, biopsy_row.ImageNumeric, biopsy_row.CICtLabel, row.CaseName, row.ImageFolder, row.RFFile, row.Phantom])))
                    break
        result = pd.DataFrame(selected, columns=columns)
        if len(result) != len(biopsy):
            raise RuntimeError(
                f"Only {len(result)} of {len(biopsy)} valid-ci/ct workbook images had an eligible RF, ROI, and phantom"
            )
        if result.duplicated(["PatientID", "VisitType", "VisitNumber", "ImageNumber"]).any():
            raise RuntimeError("Selected cohort contains duplicate patient/visit/image keys")
    if result.empty:
        raise RuntimeError("No eligible cases were found")
    return result.reset_index(drop=True)


def initialize_dimension_groups(selected: pd.DataFrame, config: Config, output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected = selected.copy()
    dimensions: list[tuple[int, int]] = []
    for phantom in selected.Phantom:
        header = read_raw_header(config.phantom_root / f"{phantom}_rf.raw")
        dimensions.append((header.samples, header.lines))
    selected["MapRows"] = [item[0] for item in dimensions]
    selected["MapColumns"] = [item[1] for item in dimensions]
    selected["MapGroup"] = 0
    selected["GroupCaseIndex"] = 0
    selected["GroupCaseCount"] = 0
    selected["GroupFolder"] = ""
    selected["GroupHDF5File"] = ""
    groups: list[dict[str, Any]] = []
    map_group = 0
    dimension_values = sorted(set(dimensions))
    for dimension_group, (rows, columns) in enumerate(dimension_values, start=1):
        indexes = selected.index[(selected.MapRows == rows) & (selected.MapColumns == columns)].tolist()
        bytes_per_case = rows * columns * len(MAP_ORDER) * 4
        payload_limit = int(config.max_folder_bytes * 0.95)
        per_folder = max(1, payload_limit // bytes_per_case)
        parts = math.ceil(len(indexes) / per_folder)
        for part in range(parts):
            map_group += 1
            group_indexes = indexes[part * per_folder : (part + 1) * per_folder]
            folder_name = f"maps_{rows}x{columns}" if parts == 1 else f"maps_{rows}x{columns}_part{part + 1:03d}"
            folder = output_dir / folder_name
            folder.mkdir(parents=True, exist_ok=True)
            h5_path = folder / "all_parametric_maps.h5"
            selected.loc[group_indexes, "MapGroup"] = map_group
            selected.loc[group_indexes, "GroupCaseIndex"] = np.arange(1, len(group_indexes) + 1)
            selected.loc[group_indexes, "GroupCaseCount"] = len(group_indexes)
            selected.loc[group_indexes, "GroupFolder"] = str(folder)
            selected.loc[group_indexes, "GroupHDF5File"] = str(h5_path)
            group_selected = selected.loc[group_indexes]
            with h5py.File(h5_path, "w") as handle:
                handle.create_dataset(
                    "all_parametric_maps",
                    shape=(rows, columns, len(MAP_ORDER), len(group_indexes)),
                    dtype="f4", chunks=(min(rows, 256), min(columns, 32), len(MAP_ORDER), 1),
                    fillvalue=np.nan,
                )
                strings = h5py.string_dtype("utf-8")
                handle.create_dataset("map_order", data=np.asarray(MAP_ORDER, dtype=strings))
                handle.create_dataset("map_dimension_order", data=np.asarray(MAP_DIMENSION_ORDER, dtype=strings))
                handle.attrs["aggregate_format_version"] = 1
                handle.attrs["logical_dimension_order"] = "axial,lateral,parameter,case"
                handle.attrs["source_matlab_script"] = "generate_mbf_ss_si_power_spectra.m"
                handle.create_group("cases")
            group_selected.to_csv(folder / "selected_cases.csv", index=False)
            groups.append({
                "MapGroup": map_group, "DimensionGroup": dimension_group,
                "PartIndex": part + 1, "NumParts": parts, "MapRows": rows,
                "MapColumns": columns, "NumCases": len(group_indexes),
                "BytesPerCase": bytes_per_case,
                "EstimatedMapBytes": bytes_per_case * len(group_indexes),
                "MaxFolderBytes": config.max_folder_bytes,
                "MapPayloadLimitBytes": payload_limit,
                "Phantoms": ", ".join(pd.unique(group_selected.Phantom.astype(str))),
                "Folder": str(folder), "HDF5File": str(h5_path),
            })
    return selected, pd.DataFrame(groups)


def json_ready(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Cannot JSON-encode {type(value).__name__}")


def write_case_hdf5(
    h5_path: Path, slot_1based: int, map_cube: np.ndarray,
    parameter_vectors: np.ndarray | None, axial: np.ndarray | None,
    lateral: np.ndarray | None, roi: dict[str, np.ndarray] | None,
    settings: dict[str, Any] | None, metadata: dict[str, Any],
) -> None:
    with h5py.File(h5_path, "r+") as handle:
        handle["all_parametric_maps"][:, :, :, slot_1based - 1] = map_cube.astype(np.float32)
        cases = handle["cases"]
        name = f"{slot_1based:06d}"
        if name in cases:
            del cases[name]
        case = cases.create_group(name)
        if parameter_vectors is not None:
            case.create_dataset("parameter_vectors", data=parameter_vectors.astype(np.float32), compression="gzip")
        if axial is not None:
            case.create_dataset("axial_axis", data=np.asarray(axial, dtype=np.float32))
        if lateral is not None:
            case.create_dataset("lateral_axis", data=np.asarray(lateral, dtype=np.float32))
        if roi is not None:
            roi_group = case.create_group("roi_positions")
            for field in ROI_FIELDS:
                roi_group.create_dataset(field, data=np.asarray(roi[field], dtype=float))
        if settings is not None:
            case.attrs["settings_json"] = json.dumps(settings, default=json_ready, sort_keys=True)
        for key, value in metadata.items():
            case.attrs[key] = "" if value is None else value


def nanmean(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.nanmean(values)) if np.isfinite(values).any() else np.nan


def interpolate_reference_bsc(table: np.ndarray, frequency: np.ndarray) -> np.ndarray:
    function = interp1d(table[:, 0], table[:, 1], kind="linear", fill_value="extrapolate")
    return np.asarray(function(frequency), dtype=float)


def process_case(row: pd.Series, config: Config, reference_table: np.ndarray, phantom_cache: dict[str, Any]) -> tuple[dict[str, Any], np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray], dict[str, Any]]:
    image_dir = config.raw_root / str(row.CaseName) / str(row.ImageFolder)
    raw_path = image_dir / str(row.RFFile)
    sample_base = re.sub(r"\.raw$", "", str(row.RFFile), flags=re.I)
    roi_path, frame = choose_standard_roi(image_dir / "ROIs", sample_base, raw_path, config.roi_dimension_suffix)
    roi = load_roi(roi_path)
    sample_rf, sample_axial, sample_lateral, sample_fs, sample_header, metadata, timestamp = load_rf_data(raw_path, frame)
    old_lateral = gui_geometric_lateral_axis(metadata, sample_header.lines, sample_header.samples)
    roi = transform_roi_lateral(roi, old_lateral, sample_lateral)
    mask = None
    if config.sample_attenuation_mode == "regional":
        mask = load_kidney_mask(image_dir / "ROIs", sample_base, sample_rf.shape)
    if config.apply_tgc_removal:
        sample_rf = apply_clarius_tgc(sample_rf, sample_axial, tgc_path_for(raw_path), timestamp)
    if config.apply_attenuation_correction:
        if config.sample_attenuation_mode == "regional":
            sample_rf = apply_stft_attenuation(sample_rf, sample_fs, sample_axial, mask, config.non_kidney_alpha, config.kidney_alpha, config.stft_taper)
        else:
            sample_rf = apply_stft_attenuation(sample_rf, sample_fs, sample_axial, np.zeros_like(sample_rf, dtype=bool), config.sample_alpha, config.sample_alpha, config.stft_taper)
    if len(roi["left"]) == 0:
        raise ValueError(f"Saved GUI ROI contains no windows: {roi_path}")

    phantom_name = str(row.Phantom)
    if phantom_name not in phantom_cache:
        phantom_path = config.phantom_root / f"{phantom_name}_rf.raw"
        phantom_rf, phantom_axial, phantom_lateral, phantom_fs, _, _, phantom_timestamp = load_rf_data(phantom_path, 1)
        if config.apply_tgc_removal:
            phantom_rf = apply_clarius_tgc(phantom_rf, phantom_axial, tgc_path_for(phantom_path), phantom_timestamp)
        if config.apply_attenuation_correction:
            phantom_rf = apply_stft_attenuation(phantom_rf, phantom_fs, phantom_axial, np.zeros_like(phantom_rf, dtype=bool), config.phantom_alpha, config.phantom_alpha, config.stft_taper)
        phantom_cache[phantom_name] = (phantom_rf, phantom_axial, phantom_lateral, phantom_fs)
    phantom_rf, phantom_axial, phantom_lateral, phantom_fs = phantom_cache[phantom_name]

    sample_rf, sample_axial, sample_lateral = match_sample_to_reference_grid(
        sample_rf, sample_axial, sample_lateral, phantom_rf.shape,
        phantom_axial, phantom_lateral,
    )
    roi = keep_rois_inside_grid(roi, sample_axial, sample_lateral)
    if len(roi["left"]) == 0:
        raise ValueError("No ROI windows remain after matching to the phantom grid")
    frequency_ref, phantom_spectra = phantom_nps(
        phantom_rf, phantom_axial, phantom_lateral, phantom_fs, roi, config.nfft
    )
    frequency, sample_spectra = phantom_nps(
        sample_rf, sample_axial, sample_lateral, sample_fs, roi, config.nfft
    )
    if not np.allclose(frequency, frequency_ref, rtol=0, atol=1e-9):
        raise ValueError("Sample and phantom frequency axes do not match")
    wavelength = config.sound_speed / (metadata["transmit_frequency_mhz"] * 1e6) * 1000.0
    f0, f1 = estimator_bandwidth(
        phantom_rf, sample_axial, sample_lateral, frequency_ref,
        metadata["focal_depth_mm"], wavelength, config.nfft, config.cutoff_db,
    )
    band = slice(f0, f1 + 1)
    reference_bsc = interpolate_reference_bsc(reference_table, frequency[band])
    normalized_ratio, bsc = reference_phantom_bsc(sample_spectra, phantom_spectra, reference_bsc, band)
    mbf, coefficients, r2, _ = spectral_analysis(normalized_ratio, frequency, f0, f1)
    ss, si = coefficients[:, 0], coefficients[:, 1]
    esd, eac, gaussian_r2, gaussian_chi2 = gaussian_form_factor_fit(bsc, frequency[band], config.sound_speed)
    vectors = np.column_stack((esd, eac, mbf, ss, si))
    maps = [create_estimator_map(sample_rf.shape, sample_axial, sample_lateral, roi, values) for values in (esd, eac, mbf, ss, si)]
    cube = np.stack(maps, axis=2).astype(np.float32)
    settings = {
        "N_fft": config.nfft, "window": "hann", "cutoff": config.cutoff_db,
        "c": config.sound_speed, "harm": 1, "f0": frequency[f0], "f1": frequency[f1],
        "freq": metadata["transmit_frequency_mhz"] * 1e6,
        "lambda": wavelength, "fd": metadata["focal_depth_mm"],
        "ROI_size_axial": config.roi_size_axial_mm,
        "ROI_size_lateral": config.roi_size_lateral_mm,
        "window_overlap_axial": config.roi_overlap_axial,
        "window_overlap_lateral": config.roi_overlap_lateral,
        "use_tgc": int(config.apply_tgc_removal), "atten_correct": 0,
        "use_atten_const_freq": int(config.apply_attenuation_correction),
        "sample_attenuation_mode": config.sample_attenuation_mode,
        "sample_atten": config.sample_alpha,
        "sample_non_kidney_atten": config.non_kidney_alpha,
        "sample_kidney_atten": config.kidney_alpha,
        "ref_atten": config.phantom_alpha,
        "attenuation_correction_method": "STFT_RF",
    }
    sample_id = f"{sample_base}_{frame}"
    result = {
        "PatientID": row.PatientID, "VisitType": row.VisitType,
        "VisitNumber": row.VisitNumber, "ImageNumber": row.ImageNumber,
        "CICt": row.CICtLabel, "CaseName": row.CaseName,
        "RFFile": row.RFFile, "Phantom": row.Phantom,
        "TGCRemoval": config.apply_tgc_removal,
        "AttenuationCorrection": config.apply_attenuation_correction,
        "SampleAttenuationMode": config.sample_attenuation_mode,
        "SampleAttenuation_dB_per_cm_MHz": config.sample_alpha,
        "NonKidneyAttenuation_dB_per_cm_MHz": config.non_kidney_alpha,
        "KidneyAttenuation_dB_per_cm_MHz": config.kidney_alpha,
        "PhantomAttenuation_dB_per_cm_MHz": config.phantom_alpha,
        "RFFrame": frame, "ROIFile": str(roi_path), "NumROIs": len(roi["left"]),
        "F0_MHz": frequency[f0], "F1_MHz": frequency[f1],
        "MBF_dB": nanmean(mbf), "SS_dB_per_MHz": nanmean(ss),
        "SI_dB": nanmean(si), "MeanR2": nanmean(r2),
        "GaussianESD_um": nanmean(esd),
        "GaussianEAC_dB_per_cm": nanmean(eac),
        "GaussianMeanR2": nanmean(gaussian_r2),
        "GaussianMeanChi2": nanmean(gaussian_chi2),
        "SampleID": sample_id, "Status": "completed", "ErrorMessage": "",
    }
    return result, cube, vectors, sample_axial, sample_lateral, roi, settings


def failure_result(row: pd.Series, config: Config, message: str, sample_id: str, frame: float = np.nan, roi_path: str = "") -> dict[str, Any]:
    result = {name: np.nan for name in MANIFEST_COLUMNS}
    result.update({
        "PatientID": row.PatientID, "VisitType": row.VisitType,
        "VisitNumber": row.VisitNumber, "ImageNumber": row.ImageNumber,
        "CICt": row.CICtLabel, "CaseName": row.CaseName,
        "RFFile": row.RFFile, "Phantom": row.Phantom,
        "TGCRemoval": config.apply_tgc_removal,
        "AttenuationCorrection": config.apply_attenuation_correction,
        "SampleAttenuationMode": config.sample_attenuation_mode,
        "SampleAttenuation_dB_per_cm_MHz": config.sample_alpha,
        "NonKidneyAttenuation_dB_per_cm_MHz": config.non_kidney_alpha,
        "KidneyAttenuation_dB_per_cm_MHz": config.kidney_alpha,
        "PhantomAttenuation_dB_per_cm_MHz": config.phantom_alpha,
        "RFFrame": frame, "ROIFile": roi_path, "SampleID": sample_id,
        "Status": "failed", "ErrorMessage": message,
    })
    return result


def checkpoint(manifest: list[dict[str, Any]], group_info: pd.DataFrame, output_dir: Path) -> None:
    frame = pd.DataFrame(manifest).reindex(columns=MANIFEST_COLUMNS)
    frame.to_csv(output_dir / "all_case_parameter_summary.csv", index=False)
    for group in group_info.itertuples():
        frame[frame.MapGroup == group.MapGroup].to_csv(Path(group.Folder) / "case_parameter_summary.csv", index=False)


def run_worker(config: Config, selected: pd.DataFrame, chunk_start: int, chunk_end: int, run_output: Path) -> None:
    selected = selected.copy()
    selected["GlobalCaseIndex"] = np.arange(1, len(selected) + 1)
    selected = selected.iloc[chunk_start - 1 : chunk_end].reset_index(drop=True)
    chunk_number = (chunk_start - 1) // config.chunk_size + 1
    output_dir = run_output / f"chunk_{chunk_number:03d}_cases_{chunk_start:04d}_{chunk_end:04d}"
    output_dir.mkdir(parents=True, exist_ok=True)
    selected, group_info = initialize_dimension_groups(selected, config, output_dir)
    selected.to_csv(output_dir / "selected_cases.csv", index=False)
    group_info.to_csv(output_dir / "map_dimension_groups.csv", index=False)
    reference_table = np.loadtxt(config.reference_bsc_path)
    processing = selected.sort_values(["MapGroup", "GroupCaseIndex"])
    manifest: list[dict[str, Any]] = []
    phantom_cache: dict[str, Any] = {}
    for case_number, (_, row) in enumerate(processing.iterrows(), start=1):
        print(
            f"[{case_number}/{len(processing)}] P{int(row.PatientID)} "
            f"{row.VisitType} {int(row.VisitNumber):02d}, Image {int(row.ImageNumber)}",
            flush=True,
        )
        sample_base = re.sub(r"\.raw$", "", str(row.RFFile), flags=re.I)
        try:
            result, cube, vectors, axial, lateral, roi, settings = process_case(row, config, reference_table, phantom_cache)
            if cube.shape[:2] != (int(row.MapRows), int(row.MapColumns)):
                raise ValueError(f"Computed map shape {cube.shape[:2]} != expected {(row.MapRows, row.MapColumns)}")
            result.update({
                "MapRows": row.MapRows, "MapColumns": row.MapColumns,
                "MapGroup": row.MapGroup, "GroupCaseIndex": row.GroupCaseIndex,
                "MapHDF5File": row.GroupHDF5File,
            })
            write_case_hdf5(
                Path(row.GroupHDF5File), int(row.GroupCaseIndex), cube,
                vectors, axial, lateral, roi, settings,
                {
                    "sample_id": result["SampleID"], "source_raw_file": str(row.RFFile),
                    "source_raw_path": str(config.raw_root / str(row.CaseName) / str(row.ImageFolder) / str(row.RFFile)),
                    "source_rf_frame": result["RFFrame"], "source_roi_file": result["ROIFile"],
                    "processing_status": "completed", "error_message": "",
                },
            )
        except Exception as exc:
            message = f"{type(exc).__name__}: {exc} [{traceback.format_exc(limit=1).strip()}]"
            print(f"Warning: skipping {row.RFFile}: {message}", file=sys.stderr, flush=True)
            result = failure_result(row, config, message, sample_base)
            result.update({
                "MapRows": row.MapRows, "MapColumns": row.MapColumns,
                "MapGroup": row.MapGroup, "GroupCaseIndex": row.GroupCaseIndex,
                "MapHDF5File": row.GroupHDF5File,
            })
            failed_map = np.full((int(row.MapRows), int(row.MapColumns), len(MAP_ORDER)), np.nan, dtype=np.float32)
            write_case_hdf5(
                Path(row.GroupHDF5File), int(row.GroupCaseIndex), failed_map,
                None, None, None, None, None,
                {
                    "sample_id": sample_base, "source_raw_file": str(row.RFFile),
                    "source_raw_path": str(config.raw_root / str(row.CaseName) / str(row.ImageFolder) / str(row.RFFile)),
                    "source_rf_frame": np.nan, "source_roi_file": "",
                    "processing_status": "failed", "error_message": message,
                },
            )
        manifest.append(result)
        if case_number % config.chunk_size == 0 or case_number == len(processing):
            checkpoint(manifest, group_info, output_dir)
    allowed = max(2, math.ceil(0.10 * len(manifest)))
    failed = sum(item["Status"] == "failed" for item in manifest)
    if failed > allowed:
        raise RuntimeError(
            f"{failed} of {len(manifest)} cases failed, exceeding the "
            f"allowed maximum {allowed}"
        )
    print(f"Worker finished: {output_dir}")


def merge_chunk_manifests(run_output: Path) -> None:
    paths = sorted(run_output.glob("chunk_*/all_case_parameter_summary.csv"))
    if not paths:
        raise RuntimeError(f"No chunk manifests found under {run_output}")
    pd.concat([pd.read_csv(path) for path in paths], ignore_index=True).to_csv(
        run_output / "all_case_parameter_summary.csv", index=False
    )


def coordinator(config: Config, selected: pd.DataFrame, args: argparse.Namespace) -> Path:
    base = ANALYSIS_DIR / config.output_folder_name
    base.mkdir(parents=True, exist_ok=True)
    if args.start_chunk == 1:
        stem = f"run_minus18dB_{datetime.now().strftime('%Y%m%d_%H%M%S_%f')[:-3]}"
        run_output = base / stem
        collision = 0
        while run_output.exists():
            collision += 1
            run_output = base / f"{stem}_{collision:03d}"
        run_output.mkdir()
    else:
        if not args.resume_run_folder:
            raise ValueError("--resume-run-folder is required when --start-chunk > 1")
        run_output = base / args.resume_run_folder
        if not run_output.is_dir():
            raise FileNotFoundError(f"Resume run folder not found: {run_output}")
    summary_path = run_output / "chunk_run_summary.csv"
    summary: list[dict[str, Any]] = []
    if args.start_chunk > 1 and summary_path.is_file():
        previous = pd.read_csv(summary_path)
        summary = previous[previous["Chunk"] < args.start_chunk].to_dict("records")
    n_chunks = math.ceil(len(selected) / config.chunk_size)
    if not 1 <= args.start_chunk <= n_chunks:
        raise ValueError(f"--start-chunk must be within 1-{n_chunks}")
    for chunk in range(args.start_chunk, n_chunks + 1):
        first = (chunk - 1) * config.chunk_size + 1
        last = min(chunk * config.chunk_size, len(selected))
        command = [
            sys.executable, str(SCRIPT_PATH), "--chunk-start", str(first),
            "--chunk-end", str(last), "--run-output", str(run_output),
        ]
        if args.case_rf_file:
            command.extend(["--case-rf-file", args.case_rf_file])
        print(f"Starting fresh Python process for chunk {chunk}/{n_chunks}, cases {first}-{last}", flush=True)
        started = time.time()
        completed = subprocess.run(command, check=False)
        folder = run_output / f"chunk_{chunk:03d}_cases_{first:04d}_{last:04d}"
        summary.append({
            "Chunk": chunk, "FirstCase": first, "LastCase": last,
            "ExitStatus": completed.returncode,
            "ElapsedSeconds": time.time() - started, "Folder": str(folder),
        })
        pd.DataFrame(summary).to_csv(summary_path, index=False)
        if completed.returncode:
            raise RuntimeError(f"Fresh Python process for chunk {chunk} failed with status {completed.returncode}")
    merge_chunk_manifests(run_output)
    return run_output


def self_test() -> None:
    rng = np.random.default_rng(20260724)
    rows, columns, nfft = 512, 12, 1024
    axial = np.linspace(0.0, 50.0, rows)
    lateral = np.arange(columns, dtype=float)
    rf = rng.normal(size=(rows, columns))
    roi = {
        "top": np.array([10.0, 20.0]), "bottom": np.array([20.0, 30.0]),
        "left": np.array([2.0, 5.0]), "right": np.array([6.0, 9.0]),
    }
    frequency, nps = phantom_nps(rf, axial, lateral, 15.0, roi, nfft)
    assert nps.shape == (nfft, 2) and frequency.shape == (nfft,)
    ratio = np.exp(rng.normal(scale=0.05, size=nps.shape))
    mbf, coefficients, r2, db = spectral_analysis(ratio, frequency, 80, 220)
    assert mbf.shape == (2,) and coefficients.shape == (2, 2) and db.shape == (141, 2)
    band = frequency[80:221]
    synthetic_bsc = np.column_stack(
        [
            band**4 * np.exp(-0.827 * (np.pi * band * diameter / 1540.0) ** 2)
            for diameter in (35.0, 90.0)
        ]
    )
    esd, eac, fit_r2, chi2 = gaussian_form_factor_fit(
        synthetic_bsc, band, 1540.0
    )
    assert np.array_equal(esd, [35.0, 90.0])
    assert np.isfinite(eac).all() and np.isfinite(fit_r2).all() and np.isfinite(chi2).all()
    image = create_estimator_map(rf.shape, axial, lateral, roi, np.array([1.0, 3.0]))
    assert image.shape == rf.shape and np.isfinite(image).any()
    corrected = apply_stft_attenuation(rf, 15.0, axial, np.zeros_like(rf, dtype=bool), 0.5, 0.5)
    assert corrected.shape == rf.shape and np.isfinite(corrected).all()
    with tempfile.TemporaryDirectory() as temporary:
        path = Path(temporary) / "test.h5"
        with h5py.File(path, "w") as handle:
            handle.create_dataset("all_parametric_maps", shape=(rows, columns, 5, 1), dtype="f4", fillvalue=np.nan)
            handle.create_group("cases")
        write_case_hdf5(path, 1, np.zeros((rows, columns, 5)), None, axial, lateral, roi, {}, {"processing_status": "completed"})
        with h5py.File(path) as handle:
            assert handle["all_parametric_maps"].shape == (rows, columns, 5, 1)
    print("Self-test passed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case-rf-file", default=os.environ.get("QUS_CASE_RF_FILE", "").strip())
    parser.add_argument("--start-chunk", type=int, default=1)
    parser.add_argument("--resume-run-folder", default="")
    parser.add_argument("--chunk-start", type=int, default=int(os.environ.get("QUS_CHUNK_START", "0") or 0), help=argparse.SUPPRESS)
    parser.add_argument("--chunk-end", type=int, default=int(os.environ.get("QUS_CHUNK_END", "0") or 0), help=argparse.SUPPRESS)
    parser.add_argument("--run-output", type=Path, default=Path(os.environ["QUS_OUTPUT_DIR"]) if os.environ.get("QUS_OUTPUT_DIR") else None, help=argparse.SUPPRESS)
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.self_test:
        self_test()
        return
    config = default_config()
    require_inputs(config)
    print(
        f"Cohort={config.cohort_mode}, TGC={config.apply_tgc_removal}, "
        f"attenuation={config.apply_attenuation_correction}, "
        f"sample mode={config.sample_attenuation_mode}, "
        f"kidney alpha={config.kidney_alpha:.3f}, "
        f"non-kidney alpha={config.non_kidney_alpha:.3f}, "
        f"phantom alpha={config.phantom_alpha:.3f}",
        flush=True,
    )
    selected = select_cases(config)
    if args.case_rf_file:
        selected = selected[selected.RFFile.eq(args.case_rf_file)].reset_index(drop=True)
        if len(selected) != 1:
            raise RuntimeError(f"--case-rf-file matched {len(selected)} cases instead of one")
    is_worker = args.chunk_start > 0 and args.chunk_end > 0 and args.run_output is not None
    if is_worker:
        if not (1 <= args.chunk_start <= args.chunk_end <= len(selected)):
            raise ValueError(f"Invalid worker range {args.chunk_start}-{args.chunk_end} for {len(selected)} cases")
        run_worker(config, selected, args.chunk_start, args.chunk_end, args.run_output.resolve())
    else:
        output = coordinator(config, selected, args)
        print(f"All chunks finished. Output: {output}")


if __name__ == "__main__":
    main()
