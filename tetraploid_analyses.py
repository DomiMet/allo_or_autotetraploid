#!/usr/bin/env python3

"""
tetraploid_analysis.py – Analyze allele and genotype frequencies of tetraploid SNP data
with optional binning of observed frequencies for smoother curves.

This script reads a VCF file of tetraploid samples and computes allele
frequencies and genotype frequency distributions for selected samples. It outputs
a plot of genotype frequency vs. allele frequency in SVG file format, showing:
- rasterized observed datapoints,
- raw and smoothed mean curves, 
- theoretical curves for autotetraploid (tetrasomic) and allotetraploid (disomic) models.
"""

import argparse
import pysam
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline

# --- Parse command-line arguments ---
parser = argparse.ArgumentParser(
    description="Analyze tetraploid allele and genotype frequencies.")
parser.add_argument(
    "vcf_file",
    help="Path to the multi-sample VCF file (can be .vcf or .vcf.gz)"
)
parser.add_argument(
    "-s", "--samples",
    nargs='+', required=True,
    help="One or more sample names (as listed in VCF header) to include in analysis."
)
parser.add_argument(
    "-o", "--output",
    default="genotype_vs_allele_frequency.svg",
    help="Output SVG file name"
)
parser.add_argument(
    "--binsize", type=float, default=0.05,
    help="Bin width for allele-frequency binning (default: 0.05)"
)
parser.add_argument(
    "--dpi", type=int, default=600,
    help="DPI for rasterization of scatter points and figure resolution (default: 600)"
)
# New options for smoothing spline
parser.add_argument(
    "--smooth", action="store_true",
    help="Overlay smoothed mean curves in addition to raw means."
)
parser.add_argument(
    "--smooth-factor", type=float, default=0.01,
    help="Smoothing factor (default: 0.01)"
)
args = parser.parse_args()

vcf_path = args.vcf_file
selected_samples = args.samples
output_svg = args.output
bin_size = args.binsize
dpi = args.dpi

# --- Load VCF and validate samples ---
vcf_in = pysam.VariantFile(vcf_path)
all_samples = list(vcf_in.header.samples)
for samp in selected_samples:
    if samp not in all_samples:
        raise SystemExit(f"Sample '{samp}' not found in VCF header. Available: {', '.join(all_samples)}")

# --- Collect allele and genotype frequencies ---
allele_freqs = []
genotype_freqs = {i: [] for i in range(5)}  # 0:AAAA,1:AAAa,...,4:aaaa
for rec in vcf_in.fetch():
    # filter for biallelic SNPs
    if len(rec.ref) != 1 or len(rec.alts) != 1 or len(rec.alts[0]) != 1:
        continue
    counts = {i: 0 for i in range(5)}
    alt_sum = 0
    n_samples = 0
    for samp in selected_samples:
        gt = rec.samples[samp].get('GT')
        if not gt or any(a is None for a in gt):
            continue
        alt_ct = sum(1 for a in gt if a > 0)
        counts[alt_ct] += 1
        alt_sum += alt_ct
        n_samples += 1
    if n_samples == 0:
        continue
    q = alt_sum / (4 * n_samples)
    # skip monomorphic sites
    if q == 0.0 or q == 1.0:
        continue
    allele_freqs.append(q)
    for i in range(5):
        genotype_freqs[i].append(counts[i] / n_samples)

# --- Theoretical expectations ---
q_vals = np.linspace(0, 1, 501)
auto = {i: [] for i in range(5)}
for q in q_vals:
    auto[0].append((1 - q)**4)
    auto[1].append(4 * q * (1 - q)**3)
    auto[2].append(6 * q**2 * (1 - q)**2)
    auto[3].append(4 * q**3 * (1 - q))
    auto[4].append(q**4)
for i in auto:
    auto[i] = np.array(auto[i])
allo = {i: [] for i in range(5)}
for q in q_vals:
    if q <= 0.5:
        P = [
            (1 - 2*q)**2 if 2*q <= 1 else 0,
            4*q*(1 - 2*q) if 2*q <= 1 else 0,
            (2*q)**2 if 2*q <= 1 else 1,
            0, 0
        ]
    else:
        p = 1 - q
        P = [0, 0,
             (2*p)**2 if 2*p <= 1 else 1,
             4*p*(2*q - 1) if 2*p <= 1 else 0,
             (2*q - 1)**2 if 2*p <= 1 else 0]
    for i, val in enumerate(P):
        allo[i].append(val)
for i in allo:
    allo[i] = np.array(allo[i])

# --- Plot static SVG at double size ---
plt.close('all')
# original was figsize=(10,6); doubling to (20,12)
fig, ax = plt.subplots(figsize=(20, 12), dpi=dpi)
# rasterized observed datapoints
colors = ['blue','cyan','magenta','orange','red']
labels = ['AAAA','AAAa','AAaa','Aaaa','aaaa']
for i, label in enumerate(labels):
    ax.scatter(
        allele_freqs,
        genotype_freqs[i],
        s=4,
        alpha=0.2,
        color=colors[i],
        rasterized=True,
        label=None,
        zorder=3
    )
# binned mean frequencies (extended to [0,1])
bins = np.arange(0, 1 + bin_size, bin_size)
centers = 0.5 * (bins[:-1] + bins[1:])
for i, label in enumerate(labels):
    q_arr = np.array(allele_freqs)
    g_arr = np.array(genotype_freqs[i])
    means = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        vals = g_arr[(q_arr >= lo) & (q_arr < hi)]
        means.append(vals.mean() if vals.size else np.nan)
    means = np.array(means)
    # extend to start at q=0 and end at q=1
    start_val = means[0] if not np.isnan(means[0]) else 0
    end_val = means[-1] if not np.isnan(means[-1]) else 0
    x_vals = np.concatenate(([0], centers, [1]))
    y_vals = np.concatenate(([start_val], means, [end_val]))
    # plot raw mean line extended
    ax.plot(
        x_vals,
        y_vals,
        color=colors[i],
        linewidth=2,
        marker='o',
        label=f"Mean {label} raw",
        zorder=2
    )
    # overlay smoothed line if requested
    if args.smooth:
        valid2 = ~np.isnan(y_vals)
        if valid2.sum() >= 4:
            spline = UnivariateSpline(x_vals[valid2], y_vals[valid2], s=args.smooth_factor)
            x_smooth = np.linspace(0, 1, 300)
            y_smooth = spline(x_smooth)
            ax.plot(
                x_smooth,
                y_smooth,
                color=colors[i],
                linewidth=2.5,
                label=f"Mean {label} smooth",
                zorder=1
            )
# theoretical curves as pure vector
for i, label in enumerate(labels):
    ax.plot(
        q_vals, auto[i],
        linestyle='--',
        linewidth=2,
        color=colors[i],
        label=f"Auto {label}",
        zorder=1
    )
    ax.plot(
        q_vals, allo[i],
        linestyle='-',
        linewidth=2,
        color=colors[i],
        label=f"Allo {label}",
        zorder=1
    )
# finalize axes and legend
ax.set_title('Genotype Frequency vs Allele Frequency')
ax.set_xlabel('Allele Frequency (q)')
ax.set_ylabel('Genotype Frequency (proportion)')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.legend(loc='best', fontsize='small', ncol=2)
plt.tight_layout()
# save SVG
fig.savefig(output_svg, format='svg', bbox_inches='tight', pad_inches=0.1)
print(f"SVG saved to {output_svg} (20x12 inches @ {dpi} DPI; rasterized datapoints + vector curves)")

