# Plot FSI and PWIA bin-centering factors for all pm settings.
#
# Layout:
#   rows    -> theta_nq bins
#   column 1 -> FSI BC factors
#   column 2 -> PWIA BC factors
#
# All available pm settings are overlaid in every panel.
#
# Usage:
#   python3 plot_all_pm_bc_by_angle.py <data_set> <sys_ext>
#
# Example:
#   python3 plot_all_pm_bc_by_angle.py 0 sys_ext
#
# Expected file names:
#   pm120_laget_bc_corr.txt
#   pmXXX_laget_bc_corr_setN.txt  for the other pm settings

import sys
import site
import os
import glob
import re

# Keep the system NumPy/Matplotlib compatibility used on ifarm.
_original_sys_path = sys.path[:]
_user_site = site.getusersitepackages()
sys.path[:] = [path for path in sys.path if path != _user_site]

import numpy as np
import matplotlib.pyplot as plt

sys.path[:] = _original_sys_path

from LT.datafile import dfile


def as_float_array(values):
    """Return a writable floating-point NumPy array."""
    return np.asarray(values, dtype=float).copy()


def replace_invalid(values):
    """Replace -1 and non-finite values with NaN."""
    values = as_float_array(values)
    values[values == -1] = np.nan
    values[~np.isfinite(values)] = np.nan
    return values


def extract_pm_set(filename):
    """Extract the pm setting from a BC-correction filename."""
    match = re.search(r"/pm(\d+)_laget_bc_corr", filename.replace("\\", "/"))
    return int(match.group(1)) if match else None


def discover_input_files(input_dir, data_set):
    """
    Locate the pm120 file and every file belonging to the requested data set.

    pm120 is treated specially because the original code uses no _setN suffix.
    """
    candidates = []

    pm120_file = os.path.join(input_dir, "pm120_laget_bc_corr.txt")
    if os.path.isfile(pm120_file):
        candidates.append(pm120_file)

    set_pattern = os.path.join(
        input_dir, f"pm*_laget_bc_corr_set{data_set}.txt"
    )
    candidates.extend(glob.glob(set_pattern))

    # Remove duplicates and sort numerically by pm.
    unique_files = {}
    for filename in candidates:
        pm_set = extract_pm_set(filename)
        if pm_set is not None:
            unique_files[pm_set] = filename

    return sorted(unique_files.items())


def read_bc_file(pm_set, filename):
    """Read the quantities needed from one BC-correction text file."""
    table = dfile(filename)

    return {
        "pm_set": pm_set,
        "filename": filename,
        "pm": replace_invalid(table["yb"]),
        "theta": replace_invalid(table["xb"]),
        "fsi": replace_invalid(table["bc_fact_fsi"]),
        "fsi_err": replace_invalid(table["bc_fact_fsi_err"]),
        "pwia": replace_invalid(table["bc_fact_pwia"]),
        "pwia_err": replace_invalid(table["bc_fact_pwia_err"]),
    }


def angle_mask(theta_values, angle, tolerance=0.1):
    """
    Match an angle safely even when the text file stores floating-point values.
    A tolerance of 0.1 degree is enough for nominal integer angle bins.
    """
    return np.isfinite(theta_values) & np.isclose(
        theta_values, angle, atol=tolerance, rtol=0.0
    )


def finite_bc_mask(x, y, yerr):
    """Select finite values suitable for errorbar plotting."""
    return (
        np.isfinite(x)
        & np.isfinite(y)
        & np.isfinite(yerr)
        & (yerr >= 0)
    )


if len(sys.argv) != 3:
    sys.exit(
        "Usage: python3 plot_all_pm_bc_by_angle.py "
        "<data_set> <sys_ext>"
    )

data_set = int(sys.argv[1])
sys_ext = sys.argv[2].rstrip("/")

if not os.path.isdir(sys_ext):
    sys.exit(f"Input directory does not exist: {sys_ext}")

input_files = discover_input_files(sys_ext, data_set)

if not input_files:
    sys.exit(
        f"No BC-correction files were found in '{sys_ext}' for set {data_set}.\n"
        "Expected pm120_laget_bc_corr.txt and/or "
        f"pm*_laget_bc_corr_set{data_set}.txt"
    )

print("\nFiles included:")
for pm_set, filename in input_files:
    print(f"  pm{pm_set}: {filename}")

all_data = []
for pm_set, filename in input_files:
    try:
        all_data.append(read_bc_file(pm_set, filename))
    except Exception as error:
        print(f"WARNING: skipping pm{pm_set}: {error}")

if not all_data:
    sys.exit("None of the discovered files could be read.")

# Nominal neutron recoil-angle bins.
theta_bins = [5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105]

# One row per angle and two columns: FSI and PWIA.
nrows = len(theta_bins)
fig, axes = plt.subplots(
    nrows=nrows,
    ncols=2,
    figsize=(14, 3.0 * nrows),
    sharex=False,
    sharey=True,
    squeeze=False,
)

# Use one consistent marker for each pm setting in both columns.
markers = ["o", "s", "^", "v", "D", "P", "X", "<", ">", "*", "h", "8"]

# Used to create one common legend for the full canvas.
legend_handles = []
legend_labels = []

# Collect all plotted x values so that both columns can use a sensible range.
all_plotted_pm = []

for row, angle in enumerate(theta_bins):
    ax_fsi = axes[row, 0]
    ax_pwia = axes[row, 1]

    plotted_in_row = False

    for index, dataset in enumerate(all_data):
        pm_set = dataset["pm_set"]
        marker = markers[index % len(markers)]

        mask = angle_mask(dataset["theta"], angle)

        if not np.any(mask):
            continue

        # x is the event/bin missing momentum stored in yb.
        x = dataset["pm"][mask]

        fsi = dataset["fsi"][mask]
        fsi_err = dataset["fsi_err"][mask]

        pwia = dataset["pwia"][mask]
        pwia_err = dataset["pwia_err"][mask]

        # Sort each pm sample by its actual missing-momentum coordinate.
        order = np.argsort(x)
        x = x[order]
        fsi = fsi[order]
        fsi_err = fsi_err[order]
        pwia = pwia[order]
        pwia_err = pwia_err[order]

        fsi_good = finite_bc_mask(x, fsi, fsi_err)
        pwia_good = finite_bc_mask(x, pwia, pwia_err)

        label = rf"$p_m={pm_set}$ MeV/c"

        if np.any(fsi_good):
            container = ax_fsi.errorbar(
                x[fsi_good],
                fsi[fsi_good],
                yerr=fsi_err[fsi_good],
                fmt=marker,
                linestyle="-",
                linewidth=1.0,
                markersize=4.5,
                capsize=2,
                label=label,
            )

            if label not in legend_labels:
                legend_handles.append(container)
                legend_labels.append(label)

            all_plotted_pm.extend(x[fsi_good].tolist())
            plotted_in_row = True

        if np.any(pwia_good):
            ax_pwia.errorbar(
                x[pwia_good],
                pwia[pwia_good],
                yerr=pwia_err[pwia_good],
                fmt=marker,
                linestyle="-",
                linewidth=1.0,
                markersize=4.5,
                capsize=2,
                label=label,
            )

            all_plotted_pm.extend(x[pwia_good].tolist())
            plotted_in_row = True

    for ax in (ax_fsi, ax_pwia):
        ax.axhline(1.0, linestyle="--", linewidth=1.0)
        ax.set_ylim(0.01, 2.5)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=9)

    ax_fsi.set_ylabel(
        rf"$\theta_{{nq}}={angle}\pm5^\circ$" + "\nBC factor",
        fontsize=10,
    )

    if not plotted_in_row:
        ax_fsi.text(
            0.5, 0.5, "No data",
            transform=ax_fsi.transAxes,
            ha="center", va="center"
        )
        ax_pwia.text(
            0.5, 0.5, "No data",
            transform=ax_pwia.transAxes,
            ha="center", va="center"
        )

# Column headings.
axes[0, 0].set_title("FSI bin-centering correction", fontsize=14)
axes[0, 1].set_title("PWIA bin-centering correction", fontsize=14)

# x labels on the bottom row.
axes[-1, 0].set_xlabel(r"$p_r$ (GeV/c)", fontsize=12)
axes[-1, 1].set_xlabel(r"$p_r$ (GeV/c)", fontsize=12)

# Apply the same x range to every panel, with a small margin.
if all_plotted_pm:
    x_min = np.nanmin(all_plotted_pm)
    x_max = np.nanmax(all_plotted_pm)
    margin = 0.03 * (x_max - x_min) if x_max > x_min else 0.05

    for ax in axes.flat:
        ax.set_xlim(x_min - margin, x_max + margin)

fig.suptitle(
    f"FSI and PWIA BC Factors: all available $p_m$ settings, set {data_set}",
    fontsize=17,
    y=0.995,
)

if legend_handles:
    fig.legend(
        legend_handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.975),
        ncol=min(6, len(legend_labels)),
        fontsize=9,
        frameon=True,
    )

fig.tight_layout(rect=[0.04, 0.025, 0.99, 0.95])

output_dir = f"{sys_ext}_plots_all_momenta"
os.makedirs(output_dir, exist_ok=True)

pdf_output = os.path.join(
    output_dir,
    f"all_pm_set{data_set}_FSI_PWIA_BC_by_angle.pdf",
)

png_output = os.path.join(
    output_dir,
    f"all_pm_set{data_set}_FSI_PWIA_BC_by_angle.png",
)

fig.savefig(pdf_output, bbox_inches="tight")
fig.savefig(png_output, dpi=200, bbox_inches="tight")
plt.close(fig)

print("\nCombined FSI/PWIA BC plot saved to:")
print(os.path.abspath(pdf_output))
print(os.path.abspath(png_output))
