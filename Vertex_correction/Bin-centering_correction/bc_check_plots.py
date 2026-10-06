#This code checks the Bin-Centering corrections are done properly and make
#sense. 
#The code makes plots of: 
#1) BC corr. ratio for PWIA, FSI
#2) dataXsec before/after FSI BC. Corrected
#3) fsiXsec_avg_simc, fsiXsec_theory  (the two cross sections ratio give bc corr. factor)

import sys
import site

# The user site contains LT, but also NumPy 2.x.  The system Matplotlib on
# ifarm was compiled against NumPy 1.x, so load the compatible system NumPy
# and Matplotlib before restoring the user site and importing LT.
_original_sys_path = sys.path[:]
_user_site = site.getusersitepackages()
sys.path[:] = [path for path in sys.path if path != _user_site]

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tic

sys.path[:] = _original_sys_path

from LT.datafile import dfile
import LT.box as B
import os

def convert2NaN(arr=np.array([]), value=0):
    #method to convert a specified value in a array to nan (not a number)
    
    for i in enumerate(arr):
        if arr[i[0]]==value:
            arr[i[0]] = np.nan
    return arr


#User Input (Usage: python3 check_bc_corr.py 120 1 <sys_ext>)
if len(sys.argv) != 4:
    sys.exit('Usage: python3 check_bc_corr.py <pm_set> <data_set> <sys_ext>')

pm_set = int(sys.argv[1])
data_set = int(sys.argv[2])
syst_ext = sys.argv[3]  


#Make relevant directories to store plots
dir_name = syst_ext+"_plots_120"

#check if directory exists, else creates it.
if not os.path.exists(dir_name):
    os.makedirs(dir_name)

#Read data file
if pm_set==120:
    fname = syst_ext + '/pm%i_laget_bc_corr.txt' % (pm_set)
else:
    fname = syst_ext + '/pm%i_laget_bc_corr_set%i.txt' %(pm_set, data_set)

f = dfile(fname)

#Get Relevant headers to plot
pm = f['yb']    #missing momentum bin
thnq = f['xb']  #theta_nq bin 

#Get Bin-Centering Factors
bc_fact_fsi = f['bc_fact_fsi']
bc_fact_fsi_err = f['bc_fact_fsi_err']

bc_fact_pwia = f['bc_fact_pwia']
bc_fact_pwia_err = f['bc_fact_pwia_err']

#Get Data Xsec
fsiRC_dataXsec = f['fsiRC_dataXsec']               #only rad. corrected  data
fsiRC_dataXsec_err = f['fsiRC_dataXsec_err']
fsiRC_dataXsec_fsibc_corr = f['fsiRC_dataXsec_fsibc_corr']          #rad + b.c. corr.  data
fsiRC_dataXsec_fsibc_corr_err = f['fsiRC_dataXsec_fsibc_corr_err']

#Get Laget Theory Xsec (@ the avg. kinematics)
fsiXsec_theory = f['fsiXsec_theory']

# Reconstruct the SIMC average because calculate_bc_corr.py writes the
# bin-centering factor (theory / SIMC average), not separate SIMC columns.
with np.errstate(divide='ignore', invalid='ignore'):
    fsiXsec_SIMC_avg = np.where(
        bc_fact_fsi > 0, fsiXsec_theory / bc_fact_fsi, np.nan)
    fsiXsec_SIMC_avg_err = np.where(
        (bc_fact_fsi > 0) & (fsiXsec_theory > 0),
        bc_fact_fsi_err * fsiXsec_theory / bc_fact_fsi**2,
        np.nan)

#Convert to nan is value is -1:
convert2NaN(bc_fact_fsi, value=-1)
convert2NaN(bc_fact_fsi_err, value=-1)
convert2NaN(bc_fact_pwia, value=-1)
convert2NaN(bc_fact_pwia_err, value=-1)
convert2NaN(fsiRC_dataXsec, value=-1)
convert2NaN(fsiRC_dataXsec_err, value=-1)
convert2NaN(fsiRC_dataXsec_fsibc_corr, value=-1)
convert2NaN(fsiRC_dataXsec_fsibc_corr_err, value=-1)
convert2NaN(fsiXsec_SIMC_avg, value=-1)
convert2NaN(fsiXsec_SIMC_avg_err, value=-1)
convert2NaN(fsiXsec_theory, value=-1)

# Angles to plot
thnq_arr = [5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105]

# ============================================================
# Create three canvases, each with a 3 x 4 grid
# ============================================================

fig_data, axes_data = plt.subplots(
    3, 4,
    figsize=(18, 13),
    sharex=True,
    sharey=True
)

fig_bc, axes_bc = plt.subplots(
    3, 4,
    figsize=(18, 13),
    sharex=True,
    sharey=True
)

fig_model, axes_model = plt.subplots(
    3, 4,
    figsize=(18, 13),
    sharex=True,
    sharey=True
)

# Flatten the axes arrays so they can be indexed with one number
axes_data = axes_data.flatten()
axes_bc = axes_bc.flatten()
axes_model = axes_model.flatten()


# ============================================================
# Loop over neutron recoil angles
# ============================================================

for i, ithnq in enumerate(thnq_arr):

    print("thnq =", ithnq)

    mask = (thnq == ithnq)

    # Sort points by missing momentum so they appear in order
    order = np.argsort(pm[mask])

    pm_plot = pm[mask][order]

    fsiRC_plot = fsiRC_dataXsec[mask][order]
    fsiRC_err_plot = fsiRC_dataXsec_err[mask][order]

    fsiBC_plot = fsiRC_dataXsec_fsibc_corr[mask][order]
    fsiBC_err_plot = fsiRC_dataXsec_fsibc_corr_err[mask][order]

    bc_fsi_plot = bc_fact_fsi[mask][order]
    bc_fsi_err_plot = bc_fact_fsi_err[mask][order]

    bc_pwia_plot = bc_fact_pwia[mask][order]
    bc_pwia_err_plot = bc_fact_pwia_err[mask][order]

    simc_avg_plot = fsiXsec_SIMC_avg[mask][order]
    simc_avg_err_plot = fsiXsec_SIMC_avg_err[mask][order]

    theory_plot = fsiXsec_theory[mask][order]


    # ========================================================
    # 1. Data cross sections
    # ========================================================

    ax_data = axes_data[i]

    ax_data.errorbar(
        pm_plot,
        fsiRC_plot,
        yerr=fsiRC_err_plot,
        fmt="o",
        markersize=4,
        capsize=2,
        label="Rad. corrected"
    )

    ax_data.errorbar(
        pm_plot,
        fsiBC_plot,
        yerr=fsiBC_err_plot,
        fmt="s",
        markersize=4,
        capsize=2,
        label="Rad. + BC corrected"
    )

    ax_data.set_yscale("log")
    ax_data.set_ylim(1e-10, 1e-2)

    ax_data.set_title(
        rf"$\theta_{{nq}}={ithnq}\pm5^\circ$",
        fontsize=12
    )

    ax_data.grid(True, which="both", alpha=0.3)


    # ========================================================
    # 2. Bin-centering factors
    # ========================================================

    ax_bc = axes_bc[i]

    ax_bc.errorbar(
        pm_plot,
        bc_fsi_plot,
        yerr=bc_fsi_err_plot,
        fmt="o",
        markersize=4,
        capsize=2,
        label="FSI"
    )

    ax_bc.errorbar(
        pm_plot,
        bc_pwia_plot,
        yerr=bc_pwia_err_plot,
        fmt="s",
        markersize=4,
        capsize=2,
        label="PWIA"
    )

    # Change the BC-factor y-axis range here
    ax_bc.set_ylim(0.01, 1.8)

    # Reference line for BC factor = 1
    ax_bc.axhline(
        1.0,
        linestyle="--",
        linewidth=1
    )

    ax_bc.set_title(
        rf"$\theta_{{nq}}={ithnq}\pm5^\circ$",
        fontsize=12
    )

    ax_bc.grid(True, alpha=0.3)


    # ========================================================
    # 3. SIMC-average and theory cross sections
    # ========================================================

    ax_model = axes_model[i]

    ax_model.errorbar(
        pm_plot,
        simc_avg_plot,
        yerr=simc_avg_err_plot,
        fmt="o",
        markersize=4,
        capsize=2,
        label="Avg. SIMC Xsec"
    )

    ax_model.plot(
        pm_plot,
        theory_plot,
        marker="s",
        markersize=4,
        linestyle="none",
        label="Laget theory Xsec"
    )

    ax_model.set_yscale("log")
    ax_model.set_ylim(1e-10, 1e-2)

    ax_model.set_title(
        rf"$\theta_{{nq}}={ithnq}\pm5^\circ$",
        fontsize=12
    )

    ax_model.grid(True, which="both", alpha=0.3)


# ============================================================
# Remove unused 12th panel
# ============================================================

axes_data[-1].set_visible(False)
axes_bc[-1].set_visible(False)
axes_model[-1].set_visible(False)


# ============================================================
# Common labels and titles
# ============================================================

fig_data.suptitle(
    f"Data Cross Sections: pm{pm_set}, Set {data_set}",
    fontsize=18
)

fig_data.supxlabel(
    r"$p_{r}$ (GeV/c)",
    fontsize=15
)

fig_data.supylabel(
    r"Cross Section ($\mu$b / MeV sr$^{2}$)",
    fontsize=15
)


fig_bc.suptitle(
    f"Bin-Centering Correction Factors: pm{pm_set}, Set {data_set}",
    fontsize=18
)

fig_bc.supxlabel(
    r"$p_{r}$ (GeV/c)",
    fontsize=15
)

fig_bc.supylabel(
    "Bin-Centering Factor",
    fontsize=15
)


fig_model.suptitle(
    f"SIMC-Average and Theory Cross Sections: pm{pm_set}, Set {data_set}",
    fontsize=18
)

fig_model.supxlabel(
    r"$p_{r}$ (GeV/c)",
    fontsize=15
)

fig_model.supylabel(
    r"Cross Section ($\mu$b / MeV sr$^{2}$)",
    fontsize=15
)


# ============================================================
# Add one common legend to each canvas
# ============================================================

handles, labels = axes_data[0].get_legend_handles_labels()
fig_data.legend(
    handles,
    labels,
    loc="upper right",
    bbox_to_anchor=(0.98, 0.96),
    fontsize=11
)

handles, labels = axes_bc[0].get_legend_handles_labels()
fig_bc.legend(
    handles,
    labels,
    loc="upper right",
    bbox_to_anchor=(0.98, 0.96),
    fontsize=11
)

handles, labels = axes_model[0].get_legend_handles_labels()
fig_model.legend(
    handles,
    labels,
    loc="upper right",
    bbox_to_anchor=(0.98, 0.96),
    fontsize=11
)


# ============================================================
# Adjust spacing
# ============================================================

fig_data.tight_layout(rect=[0.03, 0.03, 0.98, 0.94])
fig_bc.tight_layout(rect=[0.03, 0.03, 0.98, 0.94])
fig_model.tight_layout(rect=[0.03, 0.03, 0.98, 0.94])


# ============================================================
# Save combined PDFs
# ============================================================

data_output = os.path.join(
    dir_name,
    f"pm{pm_set}set{data_set}_dataXsecBC_all_angles.pdf"
)

bc_output = os.path.join(
    dir_name,
    f"pm{pm_set}set{data_set}_BCfactor_all_angles.pdf"
)

model_output = os.path.join(
    dir_name,
    f"pm{pm_set}set{data_set}_avgXsec_all_angles.pdf"
)

fig_data.savefig(
    data_output,
    bbox_inches="tight"
)

fig_bc.savefig(
    bc_output,
    bbox_inches="tight"
)

fig_model.savefig(
    model_output,
    bbox_inches="tight"
)


# Close figures after saving
plt.close(fig_data)
plt.close(fig_bc)
plt.close(fig_model)


# Print exact save locations
print("\nCombined plots saved to:")

print(os.path.abspath(data_output))
print(os.path.abspath(bc_output))
print(os.path.abspath(model_output))