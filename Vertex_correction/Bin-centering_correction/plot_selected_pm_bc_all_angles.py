# Overlay FSI and PWIA bin-centering factors for:
# pm = 120, 580, 800, and 900 MeV/c
#
# Each subplot corresponds to one theta_nq bin.
# For each pm setting:
#   FSI  = filled square
#   PWIA = open square
#
# Usage:
#   python3 plot_selected_pm_bc_overlay.py <data_set> <sys_ext>
#
# Example:
#   python3 plot_selected_pm_bc_overlay.py 0 sys_ext
#
# Expected filenames:
#   sys_ext/pm120_laget_bc_corr.txt
#   sys_ext/pm580_laget_bc_corr_set0.txt
#   sys_ext/pm800_laget_bc_corr_set0.txt
#   sys_ext/pm900_laget_bc_corr_set0.txt

import sys
import site
import os

# ------------------------------------------------------------
# Load compatible NumPy and Matplotlib on ifarm
# ------------------------------------------------------------
_original_sys_path = sys.path[:]
_user_site = site.getusersitepackages()
sys.path[:] = [p for p in sys.path if p != _user_site]

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path[:] = _original_sys_path

from LT.datafile import dfile


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
def to_float_array(values):
    """Convert LT.datafile values to a writable float NumPy array."""
    return np.asarray(values, dtype=float).copy()


def clean_array(values):
    """Convert invalid entries, including -1, to NaN."""
    values = to_float_array(values)
    values[values == -1] = np.nan
    values[~np.isfinite(values)] = np.nan
    return values


def get_filename(input_dir, pm_set, data_set):
    """
    Return the expected input filename.

    pm120 has no data-set suffix in the original file convention.
    """
    if pm_set == 120:
        return os.path.join(
            input_dir,
            "pm120_laget_bc_corr.txt"
        )

    return os.path.join(
        input_dir,
        f"pm{pm_set}_laget_bc_corr_set{data_set}.txt"
    )


def read_bc_file(filename, pm_set):
    """Read the BC-factor columns from one file."""
    f = dfile(filename)

    return {
        "pm_set": pm_set,
        "pm": clean_array(f["yb"]),
        "theta": clean_array(f["xb"]),
        "bc_fsi": clean_array(f["bc_fact_fsi"]),
        "bc_fsi_err": clean_array(f["bc_fact_fsi_err"]),
        "bc_pwia": clean_array(f["bc_fact_pwia"]),
        "bc_pwia_err": clean_array(f["bc_fact_pwia_err"]),
    }


def valid_errorbar_mask(x, y, yerr):
    """Return entries that can safely be passed to errorbar."""
    return (
        np.isfinite(x)
        & np.isfinite(y)
        & np.isfinite(yerr)
        & (yerr >= 0)
    )


# ------------------------------------------------------------
# Command-line input
# ------------------------------------------------------------
if len(sys.argv) != 3:
    sys.exit(
        "Usage: python3 plot_selected_pm_bc_overlay.py "
        "<data_set> <sys_ext>\n"
        "Example: python3 plot_selected_pm_bc_overlay.py 0 sys_ext"
    )

data_set = int(sys.argv[1])
syst_ext = sys.argv[2].rstrip("/")

if not os.path.isdir(syst_ext):
    sys.exit(f"Input directory does not exist: {syst_ext}")


# ------------------------------------------------------------
# Settings to include
# ------------------------------------------------------------
pm_sets = [120, 580, 800, 900]

theta_nq_bins = [
    5,15,25,35,
    45, 55,65,75,85,95,105
]


# ------------------------------------------------------------
# Read the four input files
# ------------------------------------------------------------
all_data = []

print("\nReading BC-correction files:")

for pm_set in pm_sets:
    filename = get_filename(syst_ext, pm_set, data_set)

    if not os.path.isfile(filename):
        print(f"WARNING: file not found for pm{pm_set}:")
        print(f"         {filename}")
        continue

    try:
        dataset = read_bc_file(filename, pm_set)
        all_data.append(dataset)
        print(f"  pm{pm_set}: {filename}")
    except Exception as error:
        print(f"WARNING: could not read pm{pm_set}: {error}")

if not all_data:
    sys.exit("No selected BC-correction files could be loaded.")


# ------------------------------------------------------------
# Create one subplot per theta_nq angle
# ------------------------------------------------------------
# ------------------------------------------------------------
# Create one subplot per theta_nq angle
# ------------------------------------------------------------
# Use a dynamic grid so every requested theta_nq angle gets a panel.
n_angles = len(theta_nq_bins)
ncols = 3
nrows = int(np.ceil(n_angles / ncols))

fig, axes = plt.subplots(
    nrows,
    ncols,
    figsize=(16, 3.8 * nrows),
    sharex=True,
    sharey=True,
    squeeze=False
)

axes = axes.ravel()

# One color per pm setting.
# No explicit color names are required; these come from Matplotlib's cycle.
# color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
# pm_colors = {
#     pm_set: color_cycle[i % len(color_cycle)]
#     for i, pm_set in enumerate(pm_sets)
# }

# # Small horizontal separation lets open and filled squares remain visible
# # when PWIA and FSI have exactly the same x coordinate.
# x_offset = {
#     "fsi": -0.0025,
#     "pwia": +0.0025,
# }

# all_x_values = []


# # ------------------------------------------------------------
# # Plot every pm setting in every theta_nq panel
# # ------------------------------------------------------------
# for panel_index, theta_value in enumerate(theta_nq_bins):

#     ax = axes[panel_index]
#     panel_has_data = False

#     for dataset in all_data:

#         pm_set = dataset["pm_set"]
#         color = pm_colors[pm_set]

#         # Use isclose because theta may be stored as floating point.
#         angle_mask = np.isclose(
#             dataset["theta"],
#             theta_value,
#             atol=0.1,
#             rtol=0.0
#         )

#         if not np.any(angle_mask):
#             continue

#         x = dataset["pm"][angle_mask]

#         bc_fsi = dataset["bc_fsi"][angle_mask]
#         bc_fsi_err = dataset["bc_fsi_err"][angle_mask]

#         bc_pwia = dataset["bc_pwia"][angle_mask]
#         bc_pwia_err = dataset["bc_pwia_err"][angle_mask]

#         # Sort by the actual missing-momentum coordinate.
#         order = np.argsort(x)

#         x = x[order]
#         bc_fsi = bc_fsi[order]
#         bc_fsi_err = bc_fsi_err[order]
#         bc_pwia = bc_pwia[order]
#         bc_pwia_err = bc_pwia_err[order]

#         fsi_mask = valid_errorbar_mask(
#             x,
#             bc_fsi,
#             bc_fsi_err
#         )

#         pwia_mask = valid_errorbar_mask(
#             x,
#             bc_pwia,
#             bc_pwia_err
#         )

#         # ----------------------------------------------------
#         # FSI: filled square
#         # ----------------------------------------------------
#         if np.any(fsi_mask):
#             ax.errorbar(
#                 x[fsi_mask] + x_offset["fsi"],
#                 bc_fsi[fsi_mask],
#                 yerr=bc_fsi_err[fsi_mask],
#                 fmt="s",
#                 linestyle="none",
#                 markersize=5,
#                 markerfacecolor=color,
#                 markeredgecolor=color,
#                 markeredgewidth=1.1,
#                 ecolor=color,
#                 capsize=2,
#                 elinewidth=1.0,
#                 zorder=3
#             )

#             all_x_values.extend(x[fsi_mask].tolist())
#             panel_has_data = True

#         # ----------------------------------------------------
#         # PWIA: open square
#         # ----------------------------------------------------
#         if np.any(pwia_mask):
#             ax.errorbar(
#                 x[pwia_mask] + x_offset["pwia"],
#                 bc_pwia[pwia_mask],
#                 yerr=bc_pwia_err[pwia_mask],
#                 fmt="s",
#                 linestyle="none",
#                 markersize=5,
#                 markerfacecolor="none",
#                 markeredgecolor=color,
#                 markeredgewidth=1.3,
#                 ecolor=color,
#                 capsize=2,
#                 elinewidth=1.0,
#                 zorder=4
#             )

#             all_x_values.extend(x[pwia_mask].tolist())
#             panel_has_data = True

#     # BC factor = 1 reference line
#     ax.axhline(
#         1.0,
#         linestyle="--",
#         linewidth=1.0,
#         color="black",
#         alpha=0.7
#     )

#     # Requested y-axis range
#     ax.set_ylim(0.001, 2.5)

#     ax.set_title(
#         rf"$\theta_{{nq}}={theta_value}\pm5^\circ$",
#         fontsize=12
#     )

#     ax.grid(True, alpha=0.3)

#     if not panel_has_data:
#         ax.text(
#             0.5,
#             0.5,
#             "No data",
#             transform=ax.transAxes,
#             ha="center",
#             va="center",
#             fontsize=10
#         )


# # ------------------------------------------------------------
# # Hide panels that do not correspond to a requested angle
# # ------------------------------------------------------------
# for unused_ax in axes[len(theta_nq_bins):]:
#     unused_ax.set_visible(False)


# # ------------------------------------------------------------
# # Set a common x-axis range
# # ------------------------------------------------------------
# if all_x_values:
#     x_min = np.nanmin(all_x_values)
#     x_max = np.nanmax(all_x_values)

#     if x_max > x_min:
#         margin = 0.04 * (x_max - x_min)
#     else:
#         margin = 0.05

#     for ax in axes[:len(theta_nq_bins)]:
#         ax.set_xlim(
#             x_min - margin,
#             x_max + margin
#         )


# # ------------------------------------------------------------
# # Common labels and title
# # ------------------------------------------------------------
# fig.suptitle(
#     (
#         "FSI and PWIA Bin-Centering Factors\n"
#         r"$p_m=120,\ 580,\ 800,\ 900$ MeV/c"
#     ),
#     fontsize=18,
#     y=0.98
# )

# fig.supxlabel(
#     r"$p_r$ (GeV/c)",
#     fontsize=15
# )

# fig.supylabel(
#     "Bin-Centering Factor",
#     fontsize=15
# )


# # ------------------------------------------------------------
# # Build legends
# # ------------------------------------------------------------

# # Legend 1: filled/open square meaning
# factor_handles = [
#     Line2D(
#         [0],
#         [0],
#         marker="s",
#         linestyle="none",
#         markersize=7,
#         markerfacecolor="black",
#         markeredgecolor="black",
#         label="FSI"
#     ),
#     Line2D(
#         [0],
#         [0],
#         marker="s",
#         linestyle="none",
#         markersize=7,
#         markerfacecolor="none",
#         markeredgecolor="black",
#         markeredgewidth=1.3,
#         label="PWIA"
#     )
# ]

# factor_legend = fig.legend(
#     handles=factor_handles,
#     loc="upper right",
#     bbox_to_anchor=(0.985, 0.965),
#     fontsize=11,
#     title="BC model",
#     title_fontsize=11
# )

# fig.add_artist(factor_legend)

# # Legend 2: color associated with each pm setting
# pm_handles = []

# for pm_set in pm_sets:
#     pm_handles.append(
#         Line2D(
#             [0],
#             [0],
#             marker="s",
#             linestyle="none",
#             markersize=7,
#             markerfacecolor=pm_colors[pm_set],
#             markeredgecolor=pm_colors[pm_set],
#             label=rf"$p_m={pm_set}$ MeV/c"
#         )
#     )

# fig.legend(
#     handles=pm_handles,
#     loc="upper right",
#     bbox_to_anchor=(0.985, 0.82),
#     fontsize=10,
#     title=r"$p_m$ setting",
#     title_fontsize=11
# )


# # ------------------------------------------------------------
# # Adjust spacing
# # ------------------------------------------------------------
# fig.tight_layout(
#     rect=[0.05, 0.05, 0.84, 0.94]
# )


# # ------------------------------------------------------------
# # Save PDF and PNG
# # ------------------------------------------------------------
# output_dir = f"{syst_ext}_plots_selected_pm"
# os.makedirs(output_dir, exist_ok=True)

# pdf_output = os.path.join(
#     output_dir,
#     (
#         f"pm120_580_800_900_set{data_set}_"
#         "FSI_PWIA_BC_all_angles.pdf"
#     )
# )

# png_output = os.path.join(
#     output_dir,
#     (
#         f"pm120_580_800_900_set{data_set}_"
#         "FSI_PWIA_BC_all_angles.png"
#     )
# )

# fig.savefig(
#     pdf_output,
#     bbox_inches="tight"
# )

# fig.savefig(
#     png_output,
#     dpi=200,
#     bbox_inches="tight"
# )

# plt.close(fig)


# # ------------------------------------------------------------
# # Print output locations
# # ------------------------------------------------------------
# print("\nPlots saved to:")
# print(os.path.abspath(pdf_output))
# print(os.path.abspath(png_output))




# ------------------------------------------------------------
# Create one separate 10 x 10 plot for each theta_nq angle
# ------------------------------------------------------------

# One color for each pm setting
color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

pm_colors = {
    pm_set: color_cycle[i % len(color_cycle)]
    for i, pm_set in enumerate(pm_sets)
}

# Slight horizontal offset so FSI and PWIA markers do not overlap
x_offset = {
    "fsi": -0.0025,
    "pwia": +0.0025,
}


# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------
output_dir = f"{syst_ext}_plots_selected_pm"
os.makedirs(output_dir, exist_ok=True)

saved_files = []


# ------------------------------------------------------------
# Make one figure for each theta_nq bin
# ------------------------------------------------------------
for theta_value in theta_nq_bins:

    # One independent 10 x 10 inch figure
    fig, ax = plt.subplots(figsize=(10, 10))

    panel_has_data = False
    angle_x_values = []

    for dataset in all_data:

        pm_set = dataset["pm_set"]
        color = pm_colors[pm_set]

        # Select entries corresponding to the current theta_nq angle
        angle_mask = np.isclose(
            dataset["theta"],
            theta_value,
            atol=0.1,
            rtol=0.0
        )

        if not np.any(angle_mask):
            continue

        x = dataset["pm"][angle_mask]

        bc_fsi = dataset["bc_fsi"][angle_mask]
        bc_fsi_err = dataset["bc_fsi_err"][angle_mask]

        bc_pwia = dataset["bc_pwia"][angle_mask]
        bc_pwia_err = dataset["bc_pwia_err"][angle_mask]

        # Sort according to the x coordinate
        order = np.argsort(x)

        x = x[order]
        bc_fsi = bc_fsi[order]
        bc_fsi_err = bc_fsi_err[order]
        bc_pwia = bc_pwia[order]
        bc_pwia_err = bc_pwia_err[order]

        fsi_mask = valid_errorbar_mask(
            x,
            bc_fsi,
            bc_fsi_err
        )

        pwia_mask = valid_errorbar_mask(
            x,
            bc_pwia,
            bc_pwia_err
        )

        # ----------------------------------------------------
        # FSI: filled square
        # ----------------------------------------------------
        if np.any(fsi_mask):

            ax.errorbar(
                x[fsi_mask] + x_offset["fsi"],
                bc_fsi[fsi_mask],
                yerr=bc_fsi_err[fsi_mask],
                fmt="s",
                linestyle="none",
                markersize=7,
                markerfacecolor=color,
                markeredgecolor=color,
                markeredgewidth=1.2,
                ecolor=color,
                capsize=3,
                elinewidth=1.2,
                zorder=3
            )

            angle_x_values.extend(x[fsi_mask].tolist())
            panel_has_data = True

        # ----------------------------------------------------
        # PWIA: open square
        # ----------------------------------------------------
        if np.any(pwia_mask):

            ax.errorbar(
                x[pwia_mask] + x_offset["pwia"],
                bc_pwia[pwia_mask],
                yerr=bc_pwia_err[pwia_mask],
                fmt="s",
                linestyle="none",
                markersize=7,
                markerfacecolor="none",
                markeredgecolor=color,
                markeredgewidth=1.5,
                ecolor=color,
                capsize=3,
                elinewidth=1.2,
                zorder=4
            )

            angle_x_values.extend(x[pwia_mask].tolist())
            panel_has_data = True


    # --------------------------------------------------------
    # Reference line at BC factor = 1
    # --------------------------------------------------------
    # ax.axhline(
    #     1.0,
    #     linestyle="--",
    #     linewidth=1.2,
    #     color="black",
    #     alpha=0.7
    # )

# --------------------------------------------------------
# Reference lines at 1, ±20%, and ±50%
# --------------------------------------------------------

# Central reference: BC factor = 1
    ax.axhline(
        y=1.0,
        linestyle=":",
        linewidth=1.5,
        color="black",
        alpha=0.9,
        label="BC factor = 1"
    )

    # ±20% lines: 0.8 and 1.2
    for y_value in [0.8, 1.2]:
        ax.axhline(
            y=y_value,
            linestyle=":",
            linewidth=1.2,
            color="black",
            alpha=0.6
        )

    # ±50% lines: 0.5 and 1.5
    for y_value in [0.5, 1.5]:
        ax.axhline(
            y=y_value,
            linestyle=":",
            linewidth=1.2,
            color="black",
            alpha=0.4
        )

    # --------------------------------------------------------
    # Axis formatting
    # --------------------------------------------------------
    ax.set_ylim(0.001, 2.5)

    if angle_x_values:

        x_min = np.nanmin(angle_x_values)
        x_max = np.nanmax(angle_x_values)

        if x_max > x_min:
            margin = 0.05 * (x_max - x_min)
        else:
            margin = 0.05

        ax.set_xlim(
            x_min - margin,
            x_max + margin
        )

    ax.set_xlabel(
        r"$p_r$ (GeV/c)",
        fontsize=16
    )

    ax.set_ylabel(
        "Bin-Centering Factor",
        fontsize=16
    )

    ax.set_title(
        (
            "FSI and PWIA Bin-Centering Factors\n"
            rf"$\theta_{{nq}}={theta_value}\pm5^\circ$"
        ),
        fontsize=18
    )

    ax.tick_params(
        axis="both",
        which="major",
        labelsize=14
    )

    ax.grid(
        True,
        alpha=0.3
    )


    # --------------------------------------------------------
    # Indicate when no data were found
    # --------------------------------------------------------
    if not panel_has_data:

        ax.text(
            0.5,
            0.5,
            "No data for this angle",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=14
        )


    # --------------------------------------------------------
    # Legend showing FSI and PWIA marker styles
    # --------------------------------------------------------
    factor_handles = [
        Line2D(
            [0],
            [0],
            marker="s",
            linestyle="none",
            markersize=8,
            markerfacecolor="black",
            markeredgecolor="black",
            label="FSI"
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            linestyle="none",
            markersize=8,
            markerfacecolor="none",
            markeredgecolor="black",
            markeredgewidth=1.5,
            label="PWIA"
        )
    ]


    # --------------------------------------------------------
    # Legend showing the pm-setting colors
    # --------------------------------------------------------
    pm_handles = []

    for pm_set in pm_sets:

        pm_handles.append(
            Line2D(
                [0],
                [0],
                marker="s",
                linestyle="none",
                markersize=8,
                markerfacecolor=pm_colors[pm_set],
                markeredgecolor=pm_colors[pm_set],
                label=rf"$p_m={pm_set}$ MeV/c"
            )
        )


    # First legend
    factor_legend = ax.legend(
        handles=factor_handles,
        loc="upper left",
        fontsize=12,
        title="BC model",
        title_fontsize=12
    )

    ax.add_artist(factor_legend)

    # Second legend
    ax.legend(
        handles=pm_handles,
        loc="upper right",
        fontsize=12,
        title=r"$p_m$ setting",
        title_fontsize=12
    )


    # --------------------------------------------------------
    # Save this individual angle plot
    # --------------------------------------------------------
    pdf_output = os.path.join(
        output_dir,
        (
            f"pm120_580_800_900_set{data_set}_"
            f"theta_nq_{theta_value}_FSI_PWIA_BC.pdf"
        )
    )

    png_output = os.path.join(
        output_dir,
        (
            f"pm120_580_800_900_set{data_set}_"
            f"theta_nq_{theta_value}_FSI_PWIA_BC.png"
        )
    )

    fig.tight_layout()

    fig.savefig(
        pdf_output,
        bbox_inches="tight"
    )

    fig.savefig(
        png_output,
        dpi=200,
        bbox_inches="tight"
    )

    saved_files.append(pdf_output)
    saved_files.append(png_output)

    plt.close(fig)


# ------------------------------------------------------------
# Print all output locations
# ------------------------------------------------------------
print("\nIndividual 10 x 10 plots saved to:")

for saved_file in saved_files:
    print(os.path.abspath(saved_file))