#Code to extract the average DATA and SIMC cross sections from the 2D Pm vs. theta_nq Histos

import argparse
import sys
from pathlib import Path

import bin_info as BI

#Set proper paths to import ROOT
sys.path.append('../../../../pyroot/')
sys.path.append('/apps/root/PRO/lib/')
sys.path.append('/apps/root/PRO/')

# do the root operations directly here
import ROOT as R

import numpy as np

# from math import *
# import everything for handling 2d histograms
# import bin_info as BI
# use a corrected version (this ignores the overflow bins)
# import bin_info1 as BI
#Usage:: python3 calculate_xsec.py 900 2 averaged_xsec/

# some constants
dtr = np.pi/180.

#MeV
MP = 938.272
MN = 939.566
#MD = 1875.6127
MD = 1875.61
me = 0.51099; 


#------------------------------------------------------------
# header information for the output file
header = \
"""
# averaged kinematics results
# averaged kinematic varibles used as input to calculate the averaged 
# kinematics: Ei, omega, th_e, pf
#
#\\ xb = th_nq
#\\ yb = pm
# current header line:
#! i_b[i,0]/ i_x[i,1]/ i_y[i,2]/ xb[f,3]/ yb[f,4]/ pwiaXsec[f,5]/  pwiaXsec_err[f,6]/  fsiXsec[f,7]/   fsiXsec_err[f,8]/  pwiaRC_dataXsec[f,9]/  pwiaRC_dataXsec_err[f,10]/  fsiRC_dataXsec[f,11]/  fsiRC_dataXsec_err[f,12]/ 
"""
#------------------------------------------------------------


def parse_args():
   parser = argparse.ArgumentParser(
      description="Extract average DATA and SIMC cross sections from ROOT histograms."
   )
   parser.add_argument("pm_set", type=int, help="missing-momentum setting, e.g. 580")
   parser.add_argument("data_set", type=int, help="data-set number, e.g. 1")
   parser.add_argument("avg_kin_ext", help="systematic-variation directory, e.g. syst_ext")
   parser.add_argument(
      "--output-dir",
      type=Path,
      help="output directory (default: avg_kin_ext relative to this script)",
   )
   return parser.parse_args()


def read_histograms(root_path, histogram_names):
   if not root_path.is_file():
      raise FileNotFoundError("ROOT input file not found: {}".format(root_path))

   root_file = R.TFile.Open(str(root_path), "READ")
   if not root_file or root_file.IsZombie():
      raise OSError("Could not open ROOT input file: {}".format(root_path))

   try:
      result = []
      for name in histogram_names:
         histogram = root_file.Get(name)
         if not histogram:
            raise KeyError("Histogram {!r} not found in {}".format(name, root_path))
         result.append(BI.get_histo_data_arrays(histogram))
      return result
   finally:
      root_file.Close()


def main():
   args = parse_args()
   pm_set = args.pm_set
   data_set = args.data_set
   avg_kin_ext = args.avg_kin_ext

   script_dir = Path(__file__).resolve().parent
   output_dir = args.output_dir.resolve() if args.output_dir else script_dir / avg_kin_ext
   output_dir.mkdir(parents=True, exist_ok=True)

   if pm_set == 120:
      output_file = output_dir / "pm{}_laget_set{}.txt".format(pm_set, data_set)
   else:
      output_file = output_dir / "pm{}_laget_set{}.txt".format(pm_set, data_set)

   analysis_dir = script_dir.parent.parent
   pm_dir = analysis_dir / "cross_sections" / "pm{}_Q2_3-4".format(pm_set)

   bin_file = pm_dir / (
      "d2_pm{}_jmlpwia_norad_"
      "xsec_output_wshmscoll_cut_diagnostics_gvill_input_file_Q2_3-4.root".format(pm_set)
   )

   cross_section_file = pm_dir / (
      "Radiative_corr_cs_combined_pm{}_"
      "wshmscut_gvill_input_file_Q2_3-4.root".format(pm_set)
   )
   
   root_file_pwia = cross_section_file
   root_file_fsi = cross_section_file

   binning, = read_histograms(bin_file, ["H_Pm_vs_thnq_v"])
   bin_info_fsiRC_dataXsec, bin_info_fsiXsec = read_histograms(
      root_file_fsi, ["h_Xsec_fsi_data_rad_corrected_PS1", "h_cross_section_simc_noradfsi_PS1"]
   )
   bin_info_pwiaRC_dataXsec, bin_info_pwiaXsec = read_histograms(
      root_file_pwia, ["h_Xsec_pwia_data_rad_corrected_PS1", "h_cross_section_simc_noradpwia_PS1"]
   )

   with output_file.open("w") as o:
      o.write('# histogram parameters \n')
      o.write('#\\ dx = {0:}\n'.format(repr(binning.dx)))
      o.write('#\\ dy = {0:}\n'.format(repr(binning.dy)))
      o.write('#\\ nx = {0:}\n'.format(repr(binning.nx)))
      o.write('#\\ ny = {0:}\n'.format(repr(binning.ny)))
      o.write('#\\ xmin = {0:}\n'.format(repr(binning.xmin)))
      o.write('#\\ ymin = {0:}\n'.format(repr(binning.ymin)))
      o.write(header)


#Retrieve the 2D bin content for each of the files
      cont1 = bin_info_fsiRC_dataXsec.cont
      cont2 = bin_info_fsiXsec.cont
      cont3 = bin_info_pwiaRC_dataXsec.cont
      cont4 = bin_info_pwiaXsec.cont

      ibin = bin_info_fsiRC_dataXsec.i
      i_xbin = bin_info_fsiRC_dataXsec.ix
      i_ybin = bin_info_fsiRC_dataXsec.iy
      thnq_b = bin_info_fsiRC_dataXsec.xb
      pm_b = bin_info_fsiRC_dataXsec.yb

#Loop over bin number (xbin, ybin)->(th_nq_bin, Pm_bin)
      for i, ib in enumerate(ibin):

   
   #Check Bin Content
         if(cont1[i] == 0.):
            fsiRC_dataXsec = -1.
            fsiRC_dataXsec_err = -1.
         else:
            fsiRC_dataXsec = bin_info_fsiRC_dataXsec.cont[i]
            fsiRC_dataXsec_err = bin_info_fsiRC_dataXsec.dcont[i]

         if(cont2[i] == 0.):
            fsiXsec = -1.
            fsiXsec_err = -1.
         else:
            fsiXsec = bin_info_fsiXsec.cont[i]
            fsiXsec_err = bin_info_fsiXsec.dcont[i]

         if(cont3[i] == 0.):
            pwiaRC_dataXsec = -1.
            pwiaRC_dataXsec_err = -1.
         else:
            pwiaRC_dataXsec = bin_info_pwiaRC_dataXsec.cont[i]
            pwiaRC_dataXsec_err = bin_info_pwiaRC_dataXsec.dcont[i]

         if(cont4[i] == 0.):
            pwiaXsec = -1.
            pwiaXsec_err = -1.
         else:
            pwiaXsec = bin_info_pwiaXsec.cont[i]
            pwiaXsec_err = bin_info_pwiaXsec.dcont[i]
      

         line = "%i %i %i %f %f %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n" % (
            ib, i_xbin[i], i_ybin[i], thnq_b[i], pm_b[i], pwiaXsec,
            pwiaXsec_err, fsiXsec, fsiXsec_err, pwiaRC_dataXsec,
            pwiaRC_dataXsec_err, fsiRC_dataXsec, fsiRC_dataXsec_err
         )
         o.write(line)

   print("Wrote {}".format(output_file))


if __name__ == "__main__":
   main()
