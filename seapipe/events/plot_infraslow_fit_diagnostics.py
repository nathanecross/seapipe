#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fit 2-Gaussian models to infraslow PSDs and plot fits with parameter bounds.

Example:
  python plot_infraslow_fit_diagnostics.py \
    --input /Users/ncro8394/Documents/projects/healthy/derivatives/datasets/clusters_summary \
    --output /Users/ncro8394/Documents/projects/healthy/derivatives/datasets/clusters_summary/fit_plots \
    --mu-bounds 0.0075 0.04 \
    --max-plots 50
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from typing import Iterable, Tuple

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def gaussian(x: np.ndarray, a: float, mu: float, sigma: float) -> np.ndarray:
    return a * np.exp(-((x - mu) ** 2) / (2.0 * sigma ** 2))


def multi_gaussian_offset(x: np.ndarray, offset: float, *params: float) -> np.ndarray:
    y = np.full_like(x, offset, dtype=float)
    for i in range(0, len(params), 3):
        a, mu, sigma = params[i:i + 3]
        y += gaussian(x, a, mu, sigma)
    return y


def fit_gauss2_with_bounds(freqs: np.ndarray, psd: np.ndarray,
                           mu_bounds: Tuple[float, float]) -> Tuple[np.ndarray, Tuple[float, ...]]:
    """Fit 2 Gaussians with constant offset and bounded mu."""
    lo, hi = mu_bounds
    init_mask = (freqs > lo) & (freqs <= hi)
    if init_mask.any():
        idx = np.nanargmax(psd[init_mask])
        mu1 = freqs[init_mask][idx]
        a1 = psd[init_mask][idx]
    else:
        idx = np.nanargmax(psd)
        mu1 = freqs[idx]
        a1 = psd[idx]

    # Seed second Gaussian inside bounds but offset from mu1.
    mu2 = mu1 + 0.01
    if mu2 > hi:
        mu2 = max(lo, mu1 - 0.01)
    if mu2 == mu1:
        mu2 = (lo + hi) / 2.0

    a2 = a1 * 0.5
    offset0 = float(np.nanmin(psd)) if not math.isnan(np.nanmin(psd)) else 0.0
    if offset0 < 0:
        offset0 = 0.0

    p0 = [offset0, a1, mu1, 0.01, a2, mu2, 0.01]
    bounds_low = [0, 0, lo, 1e-4, 0, lo, 1e-4]
    bounds_high = [np.inf, np.inf, hi, 0.1, np.inf, hi, 0.1]

    popt, _ = curve_fit(multi_gaussian_offset, freqs, psd,
                        p0=p0, bounds=(bounds_low, bounds_high), maxfev=12000)
    fit_curve = multi_gaussian_offset(freqs, *popt)
    return fit_curve, tuple(popt)


def parse_freq_header(header: Iterable[str]) -> np.ndarray:
    vals = []
    for h in header:
        try:
            vals.append(float(h))
        except ValueError:
            vals.append(np.nan)
    vals = np.array(vals, dtype=float)
    max_header = np.nanmax(vals)
    if max_header > 2.0:
        return vals * 0.001
    return vals


def iter_rows(path: str):
    with open(path, "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if not header or len(header) < 2:
            return
        freqs = parse_freq_header(header[1:])
        for row in reader:
            if not row or len(row) < 2:
                continue
            label = row[0]
            vals = []
            for v in row[1:]:
                try:
                    vals.append(float(v))
                except ValueError:
                    vals.append(np.nan)
            yield label, freqs, np.array(vals, dtype=float)


def plot_one(label: str, freqs: np.ndarray, psd: np.ndarray,
             mu_bounds: Tuple[float, float], outpath: str,
             plot_max_hz: float = 0.1) -> None:
    fit_mask = (freqs >= 0) & (freqs <= plot_max_hz)
    x = freqs[fit_mask]
    y = psd[fit_mask]

    fit_curve, params = fit_gauss2_with_bounds(x, y, mu_bounds)
    offset = params[0]
    comps = params[1:]

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    ax.plot(x, y, color="0.5", linewidth=1, label="PSD")
    ax.plot(x, fit_curve, color="r", linewidth=1.2, label="Fit")

    # Plot components
    for i in range(0, len(comps), 3):
        a, mu, sigma = comps[i:i + 3]
        ax.plot(x, offset + gaussian(x, a, mu, sigma), linestyle=":", color="r", linewidth=0.8)

    # Parameter boundaries
    ax.axvline(mu_bounds[0], color="k", linestyle=":", linewidth=0.8)
    ax.axvline(mu_bounds[1], color="k", linestyle=":", linewidth=0.8)

    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Power (a.u.)")
    ax.set_title(label)
    ax.set_xlim(0, plot_max_hz)
    ax.legend(loc="best", fontsize=8)

    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Directory containing Sigma_*_fluctuations_psd.csv")
    parser.add_argument("--output", required=True, help="Output directory for plots")
    parser.add_argument("--mu-bounds", nargs=2, type=float, default=[0.0075, 0.04])
    parser.add_argument("--max-plots", type=int, default=50, help="Max plots per file")
    parser.add_argument("--subs", default="", help="Comma-separated subject IDs to include")
    args = parser.parse_args()

    mu_bounds = (float(args.mu_bounds[0]), float(args.mu_bounds[1]))
    keep_subs = [s.strip() for s in args.subs.split(",") if s.strip()]

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    files = sorted(
        f for f in os.listdir(args.input)
        if f.startswith("Sigma_") and f.endswith("_fluctuations_psd.csv")
    )

    for fn in files:
        path = os.path.join(args.input, fn)
        count = 0
        for label, freqs, psd in iter_rows(path):
            if keep_subs and label not in keep_subs:
                continue
            outname = f"{os.path.splitext(fn)[0]}_{label}_fit.png"
            outpath = os.path.join(args.output, outname)
            try:
                plot_one(label, freqs, psd, mu_bounds, outpath)
            except Exception as exc:
                print(f"[WARN] {fn} {label}: {exc}")
                continue
            count += 1
            if count >= args.max_plots:
                break


if __name__ == "__main__":
    main()
