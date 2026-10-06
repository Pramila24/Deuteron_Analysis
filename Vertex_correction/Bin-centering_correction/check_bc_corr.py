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
sys_ext = sys.argv[3]  


dir_name = f"{sys_ext}_plots_pm{pm_set}"
os.makedirs(dir_name, exist_ok=True)

#check if directory exists, else creates it.
if not os.path.exists(dir_name):
    os.makedirs(dir_name)

#Read data file
if pm_set==120:
    fname = sys_ext + '/pm%i_laget_bc_corr.txt' % (pm_set)
else:
    fname = sys_ext + '/pm%i_laget_bc_corr_set%i.txt' %(pm_set, data_set)

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


thnq_arr = [5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105]

for i, ithnq in enumerate(thnq_arr):
    
    print('thnq=',ithnq)
         
    B.pl.clf()
    B.pl.figure(i)
    
    #Plot data Xsec
    B.plot_exp(pm[thnq==ithnq], fsiRC_dataXsec[thnq==ithnq], fsiRC_dataXsec_err[thnq==ithnq], marker='o', color='k', logy=True, label='Rad. Corrected')
    B.plot_exp(pm[thnq==ithnq], fsiRC_dataXsec_fsibc_corr[thnq==ithnq], fsiRC_dataXsec_fsibc_corr_err[thnq==ithnq], marker='o', color='r', logy=True, label='Rad. + BC. Corrected')
    
    B.pl.legend() 

    B.pl.title('Data Cross Sections, $\theta_{nq}=%i \pm 5$ deg'%(ithnq), fontsize=15)                          #Set common title
    B.pl.xlabel(r'$p_{r}$ (GeV/c)',  fontsize=12)                           #Set common x-label
    B.pl.ylabel(r'Data Cross Section $\mu b$ / MeV sr$^{2}$', fontsize=12)  #Set common y-label
    # Data cross sections — logarithmic y-axis
    B.pl.ylim(1e-10, 1e-2)

    B.pl.savefig(dir_name+'/pm%iset%i_dataXsecBC_thnq%i.pdf'%(pm_set, data_set,ithnq))
  
    #--------------PLOT BC. Corr. Factor------------------
    #-------------- PLOT BC CORRECTION FACTOR ------------------

    mask = (thnq == ithnq)

    fig, ax = plt.subplots(figsize=(8, 6))

    B.plot_exp(
    pm[mask],
    bc_fact_fsi[mask],
    bc_fact_fsi_err[mask],
    marker='o',
    color='k',
    logy=False,
    label='Bin-Centering Factor, FSI'
    )

    B.plot_exp(
    pm[mask],
    bc_fact_pwia[mask],
    bc_fact_pwia_err[mask],
    marker='o',
    color='r',
    logy=False,
    label='Bin-Centering Factor, PWIA'
    )

    ax.set_title(
    rf'Bin Centering Corr. Factor, $\theta_{{nq}}={ithnq} \pm 5$ deg',
    fontsize=15
    )

    ax.set_xlabel(r'$p_{r}$ (GeV/c)', fontsize=12)
    ax.set_ylabel('BC Factor', fontsize=12)

# Change the BC-factor y-axis range here
    ax.set_ylim(0.01, 1.8)

    ax.legend()
    ax.grid(True, alpha=0.3)

    output_file = (
    f'{dir_name}/pm{pm_set}set{data_set}_BCfactor_thnq{ithnq}.pdf'
    )

    fig.tight_layout()
    fig.savefig(output_file, bbox_inches='tight')
    plt.close(fig)

    print(f'Saved: {output_file}')
    print(f'BC-factor y limits: {ax.get_ylim()}')
    
    #Plot data Xsec
    B.plot_exp(pm[thnq==ithnq], fsiXsec_SIMC_avg[thnq==ithnq], fsiXsec_SIMC_avg_err[thnq==ithnq], marker='o', color='k', logy=True, label='Avg. SIMC Xsec')
    B.plot_exp(pm[thnq==ithnq], fsiXsec_theory[thnq==ithnq], marker='o', color='r', logy=True, label='Laget Theory Xsec')
    
    B.pl.legend() 

    B.pl.title(r'Model Cross Sections, $\theta_{nq}=%i \pm 5$ deg '%(ithnq), fontsize=15)                          #Set common title
    B.pl.xlabel(r'$p_{r}$ (GeV/c)',  fontsize=12)                           #Set common x-label
    B.pl.ylabel(r'Data Cross Section $\mu b$ / MeV sr$^{2}$', fontsize=12)  #Set common y-label
    # SIMC and theory cross sections — logarithmic y-axis
    B.pl.ylim(1e-10, 1e-2)
    B.pl.savefig(dir_name+'/pm%iset%i_avgXsec_thnq%i.pdf'%(pm_set, data_set,ithnq))
    
