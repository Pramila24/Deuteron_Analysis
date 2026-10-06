#!/usr/bin/env python3

# calculate averaged kinematics from SIMC analysis of 2D Histos (Pm vs. th_nq bins)

import sys
import os

#Set proper paths to import ROOT
sys.path.append('../../../../pyroot/')
sys.path.append('/apps/root/PRO/lib/')
sys.path.append('/apps/root/PRO/')

# do the root operations directly here
import ROOT as R

import numpy as np

hbarc = 197.327053
alpha = 1. / 137.0359895


def GMp(Q2, param=''):
   Q2 = Q2 * 1e-6
   Q = np.sqrt(Q2)
   mu_p = 2.793

   if param == 'bosted':
      return mu_p / (1. + 0.35 * Q + 2.44 * Q2 + 0.5 * Q**3 + 1.04 * Q2**2 + 0.34 * Q**5)

   if param == 'JRA':
      Q22 = Q2**2
      Q23 = Q2**3
      Q24 = Q2**4
      Q25 = Q2**5
      Q26 = Q2**6
      return mu_p / (1. + Q2 * 3.19 + Q22 * 1.355 + Q23 * 0.151 + Q24 * (-0.114E-01) + Q25 * 0.533E-03 + Q26 * (-0.900E-05))

   return mu_p / (1. + Q2 / 0.71)**2


def GEp(Q2, param=''):
   Q2 = Q2 * 1e-6
   Q = np.sqrt(Q2)

   if param == 'bosted':
      return 1. / (1. + 0.62 * Q + 0.68 * Q2 + 2.8 * Q**3 + 0.83 * Q2**2)

   if param == 'JRA':
      Q22 = Q2**2
      Q23 = Q2**3
      Q24 = Q2**4
      Q25 = Q2**5
      Q26 = Q2**6
      return 1. / (1. + Q2 * 3.226 + Q22 * 1.508 + Q23 * (-0.3773) + Q24 * 0.611 + Q25 * (-0.1853) + Q26 * 0.1596E-01)

   return 1. / (1. + Q2 / 0.71)**2


def sigMott(kf, th_e, Q2):
   if Q2 <= 0.:
      return -1.

   th_e = th_e * dtr
   sig = (2. * alpha * hbarc * kf * np.cos(th_e / 2.) / Q2)**2
   return sig * 1e4


def deForest(Ef, Q2, q, Pf, Pm, th_e, th_p, c_phi, thpq, sig_Mott, GE_p, GM_p):
   q2mu = -Q2
   q4mu = q2mu**2

   th_e = th_e * dtr
   gamma = thpq * dtr

   Epf = np.sqrt(MP * MP + Pf * Pf)
   Er = np.sqrt(MN * MN + Pm * Pm)
   Kfact = Pf * Epf

   pipf = 0.5 * (q**2 - (Pf**2 + Pm**2))
   rec = 1 - (Epf / Er) * (pipf / Pf**2)
   f_rec = 1. / rec

   Vc = q4mu / q**4
   Vt = -q2mu / (2. * q**2) + np.tan(th_e / 2.)**2
   Vi = -q2mu / q**2 * np.sqrt(-q2mu / q**2 + np.tan(th_e / 2.)**2) * c_phi
   Vs = -q2mu / q**2 * c_phi**2 + np.tan(th_e / 2.)**2

   Ebar_f = np.sqrt(Pm * Pm + MP * MP)
   om_bar = Epf - Ebar_f
   q2mu_bar = q**2 - om_bar**2
   EbarE = Ebar_f * Epf

   tau = Q2 / (4. * MP**2)
   F1 = (GE_p + tau * GM_p) / (1. + tau)
   kF2 = (GM_p - GE_p) / (1. + tau)

   Wc = 1. / (4. * EbarE) * ((Ebar_f + Epf)**2 * (F1**2 + (q2mu_bar / (4. * MP**2)) * kF2**2) - q**2 * (F1 + kF2)**2)
   Wt = q2mu_bar / (2. * EbarE) * (F1 + kF2)**2
   Ws = Pf**2 * np.sin(gamma)**2 / EbarE * (F1**2 + (q2mu_bar / (4. * MP**2)) * kF2**2)
   Wi = -Pf * np.sin(gamma) * (Ebar_f + Epf) / EbarE * (F1**2 + q2mu_bar / (4. * MP**2) * kF2**2)

   sig_eN = sig_Mott * (Vc * Wc + Vt * Wt + Vs * Ws + Vi * Wi)
   de_Forest = Kfact * f_rec * sig_eN

   return Kfact, f_rec, sig_eN, de_Forest


class HistoDataArrays:
   pass


def get_histo_data_arrays(hist):
   if not hist:
      raise RuntimeError('Missing input histogram')

   xaxis = hist.GetXaxis()
   yaxis = hist.GetYaxis()
   data = HistoDataArrays()
   data.nx = hist.GetNbinsX()
   data.ny = hist.GetNbinsY()
   data.dx = xaxis.GetBinWidth(1)
   data.dy = yaxis.GetBinWidth(1)
   data.xmin = xaxis.GetXmin()
   data.ymin = yaxis.GetXmin()
   data.i = []
   data.ix = []
   data.iy = []
   data.xb = []
   data.yb = []
   data.cont = []

   serial_bin = 0
   for ix in range(1, data.nx + 1):
      for iy in range(1, data.ny + 1):
         serial_bin += 1
         content = hist.GetBinContent(ix, iy)
         data.i.append(serial_bin)
         data.ix.append(ix)
         data.iy.append(iy)
         data.xb.append(xaxis.GetBinCenter(ix))
         data.yb.append(yaxis.GetBinCenter(iy))
         data.cont.append(content if content != 0.0 else -1)

   return data


# some constants
dtr = np.pi/180.

#MeV
MP = 938.272
MN = 939.566
#MD = 1875.6127
MD = 1875.61
me = 0.51099 

#------------------------------------------------------------
# header information for the output file
header = \
"""
# averaged kinematics results
# averaged kinematic varibles used as input to calculate the averaged 
# kinematics: Ei, omega, th_e, pf
#
# variables with _mc attached are from histograms not calculated
# alpha is the spectatror (neutron) alpha
#\\ xb = th_nq
#\\ yb = pm
# current header line:
#! i_b[i,0]/ i_x[i,1]/ i_y[i,2]/ xb[f,3]/ yb[f,4]/ Ei[f,5]/ kf[f,6]/ th_e[f,7]/ omega_mc[f,8]/ omega[f,9]/ Q2[f,10]/ Q2_calc[f,11]/ q_mc[f,12]/ q_lab[f,13]/ Ep_calc[f,14]/ pf[f,15]/ pm_mc[f,16]/ pm[f,17]/ En_calc[f,18]/ beta_cm[f,19]/ gamma_cm[f,20]/ PfPar_q[f,21]/ PfPerp_q[f,22]/ theta_pq[f,23]/ theta_pq_calc[f,24]/ PfPar_cm[f,25]/ th_pq_cm[f,26]/ th_nq_mc[f,27]/ th_nq_calc[f,28]/  cos_phi[f,29]/  sin_phi[f,30]/  alpha_c[f,31]/  GEp[f,32]/   GMp[f,33]/   sigMott[f,34]/   Ksig_cc1[f,35]/  nx[i,36]/ ny[i,37]/ cont[f,38]/        
"""
#------------------------------------------------------------
#print argv
#usage: python calc_average_kin.py 120 fsi 1 avgkin_output [JRA|dipole|bosted]

#User INput
pm_set = int(sys.argv[1])
model = sys.argv[2]
data_set = int(sys.argv[3])
sys_ext = sys.argv[4]   #systematics directory name extension  
ff_choice = sys.argv[5] if len(sys.argv) > 5 else 'JRA'

ff_param = ff_choice
if ff_choice.lower() in ('dipole', 'default'):
   ff_param = ''
   ff_label = 'dipole'
else:
   ff_label = ff_choice

if ff_label not in ('JRA', 'bosted', 'dipole'):
   print('ERROR: form-factor choice must be JRA, bosted, or dipole')
   sys.exit(1)

#Create Directory to put output if it does not exist
dir_name="./%s" % (sys_ext)

#check if directory exists, else creates it.
if not os.path.exists(dir_name):
   os.makedirs(dir_name)

#create output file to write avg kin
if pm_set == 120:
   output_file = '%s/pm%i_%s_norad_%s_avgkin.txt'%(sys_ext, pm_set, model, ff_label)
else:
   output_file = '%s/pm%i_%s_norad_%s_avgkin_set%i.txt'%(sys_ext, pm_set, model, ff_label, data_set)


o = open(output_file,'w')

#Open root file to read avg kin histos
model_tag = model
if not model_tag.startswith('jml'):
   model_tag = 'jml%s' % model_tag

root_file = '../../cross_sections/pm%i_Q2_3-4/d2_pm%i_%s_norad_xsec_output_wshmscoll_cut_diagnostics_gvill_input_file_Q2_3-4.root'%(pm_set, pm_set, model_tag)

if not os.path.exists(root_file):
   print('ERROR: ROOT input file not found: ', root_file)
   sys.exit(1)


# open ROOTfile
rf = R.TFile(root_file)

# start with 2D yield histo, Fill(Pm, thnq, FullWeight)
all = get_histo_data_arrays(rf.H_Pm_vs_thnq_v) 
# write the necessary header parameters
o.write('# histogram parameters \n')
o.write('# form factor parametrization = {0:}\n'.format(ff_label))
o.write('#\\ dx = {0:}\n'.format(repr(all.dx)))
o.write('#\\ dy = {0:}\n'.format(repr(all.dy)))
o.write('#\\ nx = {0:}\n'.format(repr(all.nx)))
o.write('#\\ ny = {0:}\n'.format(repr(all.ny)))
o.write('#\\ xmin = {0:}\n'.format(repr(all.xmin)))
o.write('#\\ ymin = {0:}\n'.format(repr(all.ymin)))
# write header
o.write(header)

#Get 2D Histogram Bin Info (Avg. kin)
bin_info_Ei        = get_histo_data_arrays(rf.H_Ein_2Davg)          #inc. beam energy [GeV]
bin_info_kf        = get_histo_data_arrays(rf.H_kf_2Davg)           #final e- momentum [GeV]
bin_info_the       = get_histo_data_arrays(rf.H_theta_elec_2Davg)   #final e- angle [deg]
bin_info_Pf        = get_histo_data_arrays(rf.H_Pf_2Davg)           #final p momentum [GeV]
bin_info_thp       = get_histo_data_arrays(rf.H_theta_prot_2Davg)   #final p angle [deg]
bin_info_q         = get_histo_data_arrays(rf.H_q_2Davg)            # |q| momentum transfer
bin_info_thq       = get_histo_data_arrays(rf.H_theta_q_2Davg)      # q-angle with +z beam
bin_info_Q2        = get_histo_data_arrays(rf.H_Q2_2Davg)           # Q2 4-momentum transfer
bin_info_nu        = get_histo_data_arrays(rf.H_omega_2Davg)        # omega, energy transfer
bin_info_xbj       = get_histo_data_arrays(rf.H_xbj_2Davg)          # Xbj, Bjorken
bin_info_Pm        = get_histo_data_arrays(rf.H_Pm_2Davg)           # Missing Momentum
bin_info_thpq      = get_histo_data_arrays(rf.H_theta_pq_2Davg)     # theta_pq [deg]
bin_info_thnq      = get_histo_data_arrays(rf.H_theta_nq_2Davg)     # theta_nq [deg]
bin_info_cphi_pq   = get_histo_data_arrays(rf.H_cphi_pq_2Davg)      # cos(phi_pq) (-1,1)
bin_info_sphi_pq   = get_histo_data_arrays(rf.H_sphi_pq_2Davg)      # sin(phi_pq) (-1,1)


#Loop over bin number (xbin, ybin)->(th_nq_bin, Pm_bin)
for i,acont in enumerate(all.cont):
   
   # get bin values
   i_bin = all.i[i]
   i_xbin = all.ix[i]
   i_ybin = all.iy[i]
   thnq_b = all.xb[i]
   pm_b = all.yb[i]
   if (acont == -1):
      # skip zero content bins
      #continue
      print('acont = ',acont)
   else:
      
      # convert rad to deg and GeV to MeV 
      Ei        = bin_info_Ei.cont[i]*1000.      
      kf        = bin_info_kf.cont[i]*1000.
      the       = bin_info_the.cont[i]
      Pf        = bin_info_Pf.cont[i]*1000.
      thp       = bin_info_thp.cont[i]
      q         = bin_info_q.cont[i]*1000.
      thq       = bin_info_thq.cont[i]
      Q2        = bin_info_Q2.cont[i]*1.e6
      nu        = bin_info_nu.cont[i]*1000.
      xbj       = bin_info_xbj.cont[i]
      Pm        = bin_info_Pm.cont[i]*1000.
      thpq      = bin_info_thpq.cont[i]
      thnq      = bin_info_thnq.cont[i]
      cphi_pq   = bin_info_cphi_pq.cont[i]
      sphi_pq   = bin_info_sphi_pq.cont[i]


      # calculate electron kinematics from measured, averaged quantities
      Ef = np.sqrt(kf*kf + me*me)  
      nu_calc = Ei - Ef

      Q2_calc = 4.*Ei*Ef*np.sin(the*dtr/2.)**2
      q_calc = np.sqrt(Q2_calc + nu_calc*nu_calc)    # |q| in the lab frame
      if q_calc==0.:
         #unphysical, skip
         continue
         # calculate hadron kinematics
      Ep = np.sqrt( MP**2 + Pf**2)
      # calculated missing momentum
      Pm_calc2 = (nu_calc+MD-Ep)**2 - MN**2
      if (Pm_calc2 < 0.):
         print('calculated pm**2 < 0. ', Pm_calc2, ' use Pm_avg : ', Pm)
         Pm_calc = Pm   #set it to the average Pm from 2D histo
      else:
         Pm_calc = np.sqrt ( Pm_calc2 )
      En_calc = np.sqrt(MN**2 + Pm_calc**2);

      # center of mass motion
      beta_cm = q_calc/(MD+nu_calc)
      gamma_cm = 1./np.sqrt(1. - beta_cm**2)

      # Momentum Components for Proton (in q-frame)
      Pf_par = ( Pf**2 + q_calc**2 - Pm_calc**2)/ (2.*q_calc)
      Pf_perp2 = Pf**2 - Pf_par**2
      if (Pf_perp2 < 0.):
         #print 'calculated Pf_perp**2<0. : ', Pf_perp2,' --->   estimate it using theta_pq :', thpq
         #print 'Pf_par = ', Pf_par, ', Pf_perp(pf) = ', Pf*np.sin(dtr*thpq), ', Pf_perp(pm) = ', Pm*np.sin(dtr*thnq)   
         Pf_perp = Pf*np.sin(dtr*thpq)
         th_pq_calc = thpq
      else:
         Pf_perp = np.sqrt(Pf_perp2)
         cthpq = Pf_par/Pf     #Cos(theta_pq)
         th_pq_calc = np.arccos(cthpq)/dtr             
      Pf_par_cm = gamma_cm*Pf_par - gamma_cm*beta_cm*Ep   #parallel component of proton in cm

      # proton angle in the cm
      thp_calc_cm = 0.

      if Pf_par_cm == 0. :
         thp_calc_cm = np.pi
      if Pf_par_cm > 0. :
         thp_calc_cm = np.arctan(Pf_perp/Pf_par_cm)
      if Pf_par_cm < 0. :
         thp_calc_cm = np.pi+np.arctan(Pf_perp/Pf_par_cm)

      theta_pq_cm = thp_calc_cm/dtr

      # calculate angles using calculated Pmiss
      denom = q_calc**2 + Pm_calc**2 - Pf**2
      num = (2.*q_calc*Pm_calc)

      cth_nq = -2.    #Cos(theta_nq)
      theta_nq_calc = -1.
      if num > 0. : 
         cth_nq = denom/num
         theta_nq_calc = 0.
      if abs(cth_nq) <=1.:
         theta_nq_calc = np.arccos(cth_nq)/dtr;
      # calculate alpha
      pz_n = Pm_calc*np.cos(theta_nq_calc*dtr)
      p_n_minus = En_calc - pz_n
      alpha_calc = p_n_minus/MN

      #Calculate the deForest Cross Section Factor  K * sig_cc1,  where K = Pf * Ep , units: ub * MeV^2 / sr^2
      sig_Mott =  sigMott(kf, the, Q2_calc)   #ub / sr
      GE_p = GEp(Q2_calc, ff_param)
      GM_p = GMp(Q2_calc, ff_param)
      Kfact, f_rec, sig_eN, de_Forest = deForest(Ef, Q2_calc, q_calc, Pf, Pm_calc, the, thp, cphi_pq, th_pq_calc, sig_Mott, GE_p, GM_p)
      
      #print('ix=',i_xbin,' iy=',i_ybin,' pm=',Pm_calc,' Kfact=',Kfact,' f_rec=',f_rec,' sig_eN=',sig_eN)
      # write output file


      l = "%i %i %i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %i %i %f\n"%( \
                                                                                                                                  # 0 2d bin number
                                                                                                                                  i_bin, \
                                                                                                                                  # 1 
                                                                                                                                  i_xbin, \
                                                                                                                                  # 2
                                                                                                                                  i_ybin, \
                                                                                                                                  # 3 central thnq_bin
                                                                                                                                  thnq_b, \
                                                                                                                                  # 4 central pm_bin
                                                                                                                                  pm_b, \
                                                                                                                                  # 5 avg. beam energy 
                                                                                                                                  Ei, \
                                                                                                                                  # 6 avg. e- momentum
                                                                                                                                  kf, \
                                                                                                                                  # 7 avg. e- angle
                                                                                                                                  the, \
                                                                                                                                  # 8 MC average energy transfer
                                                                                                                                  nu, \
                                                                                                                                  # 9 calc. average energy transer
                                                                                                                                  nu_calc, \
                                                                                                                                  # 10 MC average 4-Momentum tranfer
                                                                                                                                  Q2, \
                                                                                                                                  # 11 calc. average 4-Momentum transfer
                                                                                                                                  Q2_calc, \
                                                                                                                                  # 12 MC average |q| 3-momentum tranfer
                                                                                                                                  q, \
                                                                                                                                  # 13 calc. average |q| 3-momentum transfer
                                                                                                                                  q_calc, \
                                                                                                                                  # 14 calc. average final proton energy (assume proton mass)
                                                                                                                                  Ep, \
                                                                                                                                  # 15 MC average final proton momentum
                                                                                                                                  Pf, \
                                                                                                                                  # 16 MC average missing momentum
                                                                                                                                  Pm, \
                                                                                                                                  # 17 calc. average Missing momentum  (assume deuteron mass)
                                                                                                                                  Pm_calc, \
                                                                                                                                  # 18
                                                                                                                                  En_calc, \
                                                                                                                                  # 19
                                                                                                                                  beta_cm, \
                                                                                                                                  # 20
                                                                                                                                  gamma_cm, \
                                                                                                                                  # 21
                                                                                                                                  Pf_par, \
                                                                                                                                  # 22
                                                                                                                                  Pf_perp, \
                                                                                                                                  # 23
                                                                                                                                  thpq, \
                                                                                                                                  # 24
                                                                                                                                  th_pq_calc, \
                                                                                                                                  # 25
                                                                                                                                  Pf_par_cm, \
                                                                                                                                  # 26
                                                                                                                                  theta_pq_cm, \
                                                                                                                                  # 27
                                                                                                                                  thnq, \
                                                                                                                                  # 28
                                                                                                                                  theta_nq_calc, \
                                                                                                                                  # 29
                                                                                                                                  cphi_pq, \
                                                                                                                                  # 30
                                                                                                                                  sphi_pq, \
                                                                                                                                  # 31
                                                                                                                                  alpha_calc, \
                                                                                                                                  # 32
                                                                                                                                  GE_p, \
                                                                                                                                  # 33
                                                                                                                                  GM_p, \
                                                                                                                                  # 34
                                                                                                                                  sig_Mott, \
                                                                                                                                  # 35
                                                                                                                                  de_Forest, \
                                                                                                                                  # 36
                                                                                                                                  all.nx, \
                                                                                                                                  # 37
                                                                                                                                  all.ny, \
                                                                                                                                  # 38
                                                                                                                                  all.cont[i])
                                                                          
      o.write(l)
o.close()
