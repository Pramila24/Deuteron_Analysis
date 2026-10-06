import numpy as np
import os
import sys                             
from sys import argv    
import getopt       
from LT.datafile import dfile

'''
Code that reads averaged Xsec(data and simc) and theory cross sections.  
The ratio (sig_theory_FSI(or PWIA)/sig_avg_FSI(or PWIA) is taken to determine the bin centering
correction factor, for either model, which enables one to do systematics on model dependency of bin-centering correction. 
The bin-centering factor is then applied to the FSI radiative corrected data ONLY. The bin-centering factor using FSI
is used to determine the momentum distributions
'''


#------------------------------------------------------------                   
# header information for the output file    
header = """                                            
# averaged kinematics results                                                                              
# averaged kinematic varibles used as input to calculate the averaged                                
# kinematics: Ei, omega, th_e, pf                                                                     
#
#---Useful Header Definitions---
#bc_fact_pwia(fsi):  bin-centering correction factor using PWIA or FSI model
#pwia(fsi)RC_dataXsec: average data Xsec radiative corrected using PWIA(FSI) model --> study rad. corr systematics  __
#fsiRC_dataXsec_pwiabc_corr: data Xsec radiative corrected using FSI model and Bin-Center corrected using PWIA model  | these can be used to study bin-centering corr. systematics from model dependency
#fsiRC_dataXsec_fsibc_corr: data Xsec radiative corrected using FSI model and Bin-Center corrected using FSI model  __| -->We write FSI bin-center corr. cross section to official table
#red_dataXsec: reduced data Xsec radiative corrected and bin-center corrected using FSI model (K_sig_cc1 evaluated at FSI avg. kin is used)
#red_pwiaXsec: reduced theoretical Xsec (K_sig_cc1 evaluated at PWIA avg. kin is used)
#red_fsiXsec: reduced theoretical Xsec (K_sig_cc1 evaluated at FSI avg. kin is used)                    
#i_b = 2D bin id number
#xb = th_nq value at bin center                                                                                                   
#yb = pmiss value at bin center        
#Averaged Kinemtics Definitions
#K_sig_cc1  [ub MeV^2 / sr^2]
#cross sections [ub sr^-2 MeV^-1]
#reduced cross sections   [ub sr^-2 MeV^-1] / [ub MeV^2 / sr^2] = [MeV]^-3                                                             
#Ei_avg : averaged beam energy [GeV]
#kf_avg : averaged final e- momentum [GeV/c]
#the_avg : averaged e- scattering angle [deg]
#nu_avg : averaged energy transfer [GeV]
#q_avg : averaged 3-momentum transfer [GeV/c]
#Q2_avg : averaged 4-momentum transfer [GeV/c]^2
#pf_avg : averaged final proton momentum [GeV/c]
#pm_avg : averaged missing momentum [GeV/c]
#th_pq_cm : averaged center of mass in-plane angle between proton and q-vector [deg]
#cphi_pq_avg : cos(phi_pq), where phi_pq is the averaged out of plane angle between (pf, q) or the scattering and reaction planes [deg]
#xbj_avg : avergaed x-Bjorken scale [unitless]

# current header line:  
#! i_b[i,0]/ i_x[i,1]/ i_y[i,2]/ xb[f,3]/ yb[f,4]/ bc_fact_pwia[f,5]/  bc_fact_pwia_err[f,6]/  bc_fact_fsi[f,7]/  bc_fact_fsi_err[f,8]/  pwiaRC_dataXsec[f,9]/  pwiaRC_dataXsec_err[f,10]/  fsiRC_dataXsec[f,11]/   fsiRC_dataXsec_err[f,12]/    fsiRC_dataXsec_pwiabc_corr[f,13]/   fsiRC_dataXsec_pwiabc_corr_err[f,14]/  fsiRC_dataXsec_fsibc_corr[f,15]/  fsiRC_dataXsec_fsibc_corr_err[f,16]/   pwiaXsec_theory[f,17]/   fsiXsec_theory[f,18]/   red_dataXsec[f,19]/   red_dataXsec_err[f,20]/   red_pwiaXsec[f,21]/   red_fsiXsec[f,22]/   pwia_Ksigcc1[f,23]/   fsi_Ksigcc1[f,24]/  Ei_avg[f,25]/  kf_avg[f,26]/  the_avg[f,27]/  nu_avg[f,28]/  q_avg[f,29]/  Q2_avg[f,30]/   pf_avg[f,31]/   pm_avg[f,32]/   th_pq_cm[f,33]/   cphi_pq_avg[f,34]/   xbj_avg[f,35]/                                                                                                                                                                          
"""

#User Input
pm_set = int(sys.argv[1])                                                                            
data_set = int(sys.argv[2])
sys_ext = sys.argv[3]
theory_xsec_dir = os.path.join(
    '/w/hallc-scshelf2102/c-deuteron/pokhrelp/bc_corr_new_cuts/cuts_study_Q2_3-4/Em_n_cointime_mod/bc_diagnostics/theory_xsec/theory_output_wshmscoll/',
    'pm%i_Q2_3-4' % pm_set,
)
                                                                                                            
print(argv)
         
#usage: /apps/python/2.7.12/bin/python.py calc_bc_corr.py 580  1  corr_bc_output          

# Check if output directory exists; create it if necessary.
if not os.path.exists(sys_ext):
    os.makedirs(sys_ext)

avg_xsec_dir = (
    "/w/hallc-scshelf2102/c-deuteron/pokhrelp/"
    "/bc_corr_new_cuts/cuts_study_Q2_3-4/Em_n_cointime_mod/bc_diagnostics/xsec_avg/averaged_xsec_wshmscoll"
)

if pm_set == 120:
    output_file = os.path.join(
        sys_ext,
        f"pm{pm_set}_laget_bc_corr.txt"
    )

    # Load theory cross section at averaged kinematics.
    ft = dfile(os.path.join(
        theory_xsec_dir,
        f"pm{pm_set}_laget_theory_set{data_set}.txt"
    ))

    # Load acceptance-averaged cross section.
    fa = dfile(os.path.join(
        avg_xsec_dir,
        f"pm{pm_set}_laget_set{data_set}.txt"
    ))

else:
    output_file = os.path.join(
        sys_ext,
        f"pm{pm_set}_laget_bc_corr_set{data_set}.txt"
    )

    # Load theory cross section at averaged kinematics.
    ft = dfile(os.path.join(
        theory_xsec_dir,
        f"pm{pm_set}_laget_theory_set{data_set}.txt"
    ))

    # Load acceptance-averaged cross section.
    fa = dfile(os.path.join(
        avg_xsec_dir,
        f"pm{pm_set}_laget_set{data_set}.txt"
    ))

o = open(output_file, "w")
o.write(header)

#Get Bin Information (Does no matter which file, as they have the same binning scheme)
ib_t = ft.get_data('i_b')             #2D Bin NUmber 
ix_t = ft.get_data('i_x')             #X-Bin Number (thnq bin number)
iy_t = ft.get_data('i_y')             #Y-Bin Number (Pm bin number)
thnq_t = ft.get_data('xb')            #bin-center value of theta_nq
pm_t = ft.get_data('yb')              #bin-center value of pmiss

#Get Theory Xsec 
pwiaXsec_theory = ft.get_data('pwiaXsec') 
fsiXsec_theory = ft.get_data('fsiXsec') 
pwia_Ksig_cc1 = ft.get_data('pwia_Ksig_cc1')    #deForest K*sig_cc1 factor for calculating the reduced cross section
fsi_Ksig_cc1 = ft.get_data('fsi_Ksig_cc1')
#Get Averaged kinematics (corresponding to Laget FSI model)
Ei_avg = ft.get_data('Ei_avg') / 1000.    #[GeV]
kf_avg = ft.get_data('kf_avg') / 1000.    #[GeV]
the_avg = ft.get_data('the_avg')          #[deg]
nu_avg = ft.get_data('omega_avg') / 1000. #[GeV]
q_avg = ft.get_data('q_avg') / 1000.      #[GeV]
Q2_avg = ft.get_data('Q2_avg') / (1.e6)     #[GeV]^2
pf_avg = ft.get_data('pf_avg') / 1000.    #[GeV]
pm_avg = ft.get_data('pm_avg') / 1000.    #[GeV]
th_pq_cm = ft.get_data('th_pq_cm')        #[deg]
cphi_pq_avg = ft.get_data('cphi_pq_avg')
xbj_avg = ft.get_data('xbj_avg')


#Get Average Xsec  
pwiaXsec_avg = fa.get_data('pwiaXsec')                             #PWIA Model Average Xsec (from SIMC)                              
pwiaXsec_avg_err = fa.get_data('pwiaXsec_err')                                                                
fsiXsec_avg = fa.get_data('fsiXsec')                               #FSI Model Average Xsec (from SIMC)
fsiXsec_avg_err = fa.get_data('fsiXsec_err')
pwiaRC_dataXsec_avg = fa.get_data('pwiaRC_dataXsec')               #Radiative Corrected Data Xsec Using PWIA model
pwiaRC_dataXsec_avg_err = fa.get_data('pwiaRC_dataXsec_err')
fsiRC_dataXsec_avg = fa.get_data('fsiRC_dataXsec')                 #Radiative Corrected Data Xsec Using FSI model
fsiRC_dataXsec_avg_err = fa.get_data('fsiRC_dataXsec_err')

# The theory file can contain a filtered subset of the full averaged-Xsec bin
# grid. Its global i_b numbering also differs, so match using the shared 2D
# bin coordinates rather than row position or i_b.
ix_avg = fa.get_data('i_x')
iy_avg = fa.get_data('i_y')
avg_index_by_bin = {
    (int(ix), int(iy)): i for i, (ix, iy) in enumerate(zip(ix_avg, iy_avg))
}

#Loop over ith 2Dbin
for i, ib in enumerate(ib_t):

    try:
        j = avg_index_by_bin[(int(ix_t[i]), int(iy_t[i]))]
    except KeyError:
        raise ValueError(
            'Theory bin (i_x=%i, i_y=%i) is missing from the averaged-Xsec file'
            % (ix_t[i], iy_t[i]))

    #--------Set initial values to -1 unless the bin passes the condition-----------
    bc_factor_fsi      = -1                  #bin-centering corr. factor (using FSI model)
    bc_factor_pwia     = -1                  #bin-centering corr. factor (using PWIA model)
    bc_factor_fsi_err  = -1
    bc_factor_pwia_err = -1

    #--------Declare Bin Centering Corrected Data Xsec (that has already been Radiative Crrected by using FSI Model)----------
    fsiRC_dataXsec_fsibc_corr      = -1            #bin-centered corr. data Xsec (using FSI model for bc corr.)
    fsiRC_dataXsec_fsibc_corr_err  = -1

    fsiRC_dataXsec_pwiabc_corr     = -1            #bin-centered corr. data Xsec (using PWIA model for bc corr.)
    fsiRC_dataXsec_pwiabc_corr_err = -1

    #Reduced Xsec (Intitial set to -1)
    red_dataXsec     = -1.                #reduced dataXsec is determined from FSI radiative corrected, FSI bin-center corrected data_Xsec_avg (so use Ksig_cc1 determined from FSI avg. kinematics)
    red_dataXsec_err = -1.
    red_fsiXsec      = -1.                #theoretical FSI redXsec is determined from FSI_Xsec_theory / FSI_Ksig_cc1
    red_pwiaXsec     = -1.                #theoretical PWIA redXsec is determined from PWIA_Xsec_theory / PWIA_Ksig_cc1

    #--------------------------------------Calculate Bin Centering Correction Factor and Reduced Xsec--------------------------------------
   
    if fsiXsec_theory[i]>0. and fsiXsec_avg[j]>0.:
        bc_factor_fsi = fsiXsec_theory[i] / fsiXsec_avg[j]                                                       #FSI Bin Centering Correction Factor
        bc_factor_fsi_err = fsiXsec_theory[i] / (fsiXsec_avg[j] * fsiXsec_avg[j]) * fsiXsec_avg_err[j]

        if fsiRC_dataXsec_avg[j]>0.:
            #Apply FSI Bin Centering Correction Factor to Data Xsec (to ONLY FSI radiative corrected data, as this is what will be used in the final result, since FSI is closest to data)
            fsiRC_dataXsec_fsibc_corr = fsiRC_dataXsec_avg[j] *  bc_factor_fsi                                   #FSI Bin Centering Correction Applied (for FSI rad. corr data Xsec)
            fsiRC_dataXsec_fsibc_corr_err = np.sqrt(bc_factor_fsi**2 * (fsiRC_dataXsec_avg_err[j]*fsiRC_dataXsec_avg_err[j])  +   (fsiRC_dataXsec_avg[j]*fsiRC_dataXsec_avg[j]) *bc_factor_fsi_err**2)   
     
            #Calculate the Experimental Reduced Xsec (Momentum Distributions given the right kinematic region, thnq:(35,45) deg and Q2>4)
            if fsi_Ksig_cc1[i]>0.:
                red_dataXsec = fsiRC_dataXsec_fsibc_corr / fsi_Ksig_cc1[i]
                red_dataXsec_err = fsiRC_dataXsec_fsibc_corr_err/fsi_Ksig_cc1[i]

    if pwiaXsec_theory[i]>0. and pwiaXsec_avg[j]>0.:
        bc_factor_pwia = pwiaXsec_theory[i] / pwiaXsec_avg[j]                                                   #PWIA Bin Centering Correction Factor
        bc_factor_pwia_err = pwiaXsec_theory[i] / (pwiaXsec_avg[j] * pwiaXsec_avg[j]) * pwiaXsec_avg_err[j]

        if fsiRC_dataXsec_avg[j]>0.:
            #Apply PWIA Bin Centering Correction Factor to Data Xsec
            fsiRC_dataXsec_pwiabc_corr = fsiRC_dataXsec_avg[j] *  bc_factor_pwia                               #PWIA Bin Centering Correction Applied (for FSI rad. corr data Xsec)
            fsiRC_dataXsec_pwiabc_corr_err = np.sqrt(bc_factor_pwia**2 * (fsiRC_dataXsec_avg_err[j]*fsiRC_dataXsec_avg_err[j])  +   (fsiRC_dataXsec_avg[j]*fsiRC_dataXsec_avg[j]) *bc_factor_pwia_err**2)   
            


    #Calculate the Theoretical Reduced Cross Sections (Momentum Distributions in PWIA)
    if fsiXsec_theory[i]>0. and fsi_Ksig_cc1[i]>0.:
        red_fsiXsec =  fsiXsec_theory[i] / fsi_Ksig_cc1[i]

    if pwiaXsec_theory[i]>0. and pwia_Ksig_cc1[i]>0.:
        red_pwiaXsec =  pwiaXsec_theory[i] / pwia_Ksig_cc1[i]

    l= "%i %i %i %f %f %f %f %f %f %.12e %.12e %.12e %.12e %.12e %12e %.12e  %.12e  %.12e  %.12e   %.12e   %.12e   %.12e    %.12e    %f     %f     %f     %f     %f     %f     %f     %f     %f     %f     %f     %f     %f \n"%(ib, ix_t[i], iy_t[i], thnq_t[i], pm_t[i], bc_factor_pwia, bc_factor_pwia_err, bc_factor_fsi, bc_factor_fsi_err, pwiaRC_dataXsec_avg[j], pwiaRC_dataXsec_avg_err[j], fsiRC_dataXsec_avg[j], fsiRC_dataXsec_avg_err[j],  fsiRC_dataXsec_pwiabc_corr,  fsiRC_dataXsec_pwiabc_corr_err,  fsiRC_dataXsec_fsibc_corr,  fsiRC_dataXsec_fsibc_corr_err,  pwiaXsec_theory[i],  fsiXsec_theory[i], red_dataXsec, red_dataXsec_err,  red_pwiaXsec, red_fsiXsec,  pwia_Ksig_cc1[i],  fsi_Ksig_cc1[i], Ei_avg[i], kf_avg[i], the_avg[i], nu_avg[i], q_avg[i], Q2_avg[i], pf_avg[i], pm_avg[i], th_pq_cm[i], cphi_pq_avg[i], xbj_avg[i])               
    o.write(l)

o.close()
                        
print('Done')
