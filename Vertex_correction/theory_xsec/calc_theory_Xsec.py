#Calculate the Laget Theory Cross Section at the Averaged Kinematics

# select of reading binary grid data
use_binary = 1

# import laget module
if use_binary:
    import Laget_Xsec_fp_1 as LX
else:
    import Laget_Xsec_fp_1 as LX
        
import numpy as np
import sys
import os
import glob
from sys import argv


class SimpleDataFile:
    def __init__(self, path):
        self.path = path
        self.comments = []
        self.keys = []
        self.key_types = {}
        self.data = []
        self._read()

    def _read(self):
        with open(self.path) as f:
            for line in f:
                stripped = line.strip()
                if not stripped:
                    self.comments.append(line.rstrip())
                    continue
                if stripped.startswith('#!'):
                    self.comments.append(line.rstrip())
                    self._read_keys(stripped[2:])
                    continue
                if stripped.startswith('#'):
                    self.comments.append(line.rstrip())
                    continue

                values = stripped.split()
                if not self.keys:
                    continue
                row = {}
                for key, value in zip(self.keys, values):
                    if self.key_types.get(key) == 'i':
                        row[key] = int(float(value))
                    else:
                        row[key] = float(value)
                self.data.append(row)

    def _read_keys(self, key_line):
        self.keys = []
        self.key_types = {}
        for item in key_line.split('/'):
            item = item.strip()
            if not item or '[' not in item:
                continue
            key = item.split('[', 1)[0].strip()
            key_type = item.split('[', 1)[1].split(',', 1)[0].strip()
            self.keys.append(key)
            self.key_types[key] = key_type

    def add_key(self, key, key_type='f'):
        if key not in self.keys:
            self.keys.append(key)
            self.key_types[key] = key_type
            for row in self.data:
                row[key] = 0

    def save(self, path):
        with open(path, 'w') as out:
            for line in self.comments:
                if not line.startswith('#!'):
                    out.write(line + '\n')
            key_parts = []
            for i, key in enumerate(self.keys):
                key_parts.append('%s[%s,%i]' % (key, self.key_types.get(key, 'f'), i))
            out.write('#! ' + '/ '.join(key_parts) + '/\n')
            for row in self.data:
                values = []
                for key in self.keys:
                    if self.key_types.get(key) == 'i':
                        values.append('%i' % int(row.get(key, 0)))
                    else:
                        values.append('%.12e' % float(row.get(key, 0.0)))
                out.write(' '.join(values) + '\n')


def get_file(path):
    return SimpleDataFile(path)


def get_data(data_file, key):
    return np.array([row[key] for row in data_file.data])

#------------------------------------------------------------
# header information for the output file
header = \
"""
# Laget Cross Section Results at the Averaged Kinematics
# averaged kinematic varibles used as input to calculate the cross section 
# kinematics: Ei, Q2, omega, th_pq_cm, phi_pq
#
#\\ xb = th_nq
#\\ yb = pm
# current header line:
#! i_b[i,0]/ i_x[i,1]/ i_y[i,2]/ xb[f,3]/ yb[f,4]/  pwiaXsec[f,5]/  pwiaGEp[f,6]/   pwiaGMp[f,7]/   pwia_sigMott[f,8]/   pwia_Ksig_cc1[f,9]/    
"""
#------------------------------------------------------------
#create output file to write cross section @ avg kinematics

#usage: ipython calc_theory_Xsec.py 580 1 theory_output [JRA|dipole|bosted]

#User Input
if len(sys.argv) < 4:
    print(
        "Usage: python3 calc_theory_Xsec.py "
        "<pm_set> <data_set> <output_dir> [JRA|dipole|bosted]"
    )
    sys.exit(1)
pm_set = int(sys.argv[1])
data_set = int(sys.argv[2])
theory_output = sys.argv[3]
ff_label = sys.argv[4] if len(sys.argv) > 4 else 'JRA'
avgkin_input_dir = '/w/hallc-scshelf2102/c-deuteron/pokhrelp/bc_corr_new_cuts/cuts_study_Q2_3-4/Em_n_cointime_mod/bc_diagnostics/avg_kinematics/avg_kinematics_output_wshmscoll'
pm_dir = 'pm%i_Q2_3-4' % pm_set
theory_dir = os.path.join(theory_output, pm_dir)

#check if directory exists, else creates it.
if not os.path.exists(theory_dir):
    os.makedirs(theory_dir)

print(argv)

def avgkin_file(model):
    names = [
        'pm%i_%s_norad_%s_avgkin_set%i.txt' % (pm_set, model, ff_label, data_set),
        'pm%i_%s_norad_%s_avgkin.txt' % (pm_set, model, ff_label),
        'pm%i_%s_norad_avgkin_set%i.txt' % (pm_set, model, data_set),
        'pm%i_%s_norad_avgkin.txt' % (pm_set, model),
    ]

    candidates = []
    for name in names:
        candidates.append(os.path.join(avgkin_input_dir, pm_dir, name))
        candidates.append(os.path.join(avgkin_input_dir, name))
        candidates.append(os.path.join(theory_output, pm_dir, name))
        candidates.append(os.path.join(theory_output, name))
        candidates.append(os.path.join('..', 'avg_kinematics', theory_output, pm_dir, name))
        candidates.append(os.path.join('..', 'avg_kinematics', theory_output, name))

    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate

    # Some production runs used different set numbers for PWIA and FSI.
    # Fall back only when exactly one matching set file exists.
    fallback_patterns = [
        os.path.join(avgkin_input_dir, pm_dir,
                     'pm%i_%s_norad_%s_avgkin_set*.txt' % (pm_set, model, ff_label)),
        os.path.join(avgkin_input_dir,
                     'pm%i_%s_norad_%s_avgkin_set*.txt' % (pm_set, model, ff_label)),
    ]
    fallback_matches = sorted(set(
        match for pattern in fallback_patterns for match in glob.glob(pattern)
    ))
    if len(fallback_matches) == 1:
        print('WARNING: requested set%i %s file was not found; using %s' %
              (data_set, model, fallback_matches[0]))
        return fallback_matches[0]
    if len(fallback_matches) > 1:
        print('ERROR: multiple averaged kinematics files match %s; '
              'cannot choose a set safely:' % model)
        for match in fallback_matches:
            print('  ', match)
        sys.exit(1)

    print('ERROR: could not find averaged kinematics file for %s.' % model)
    print('Tried:')
    for candidate in candidates:
        print('  ', candidate)
    sys.exit(1)

def theory_file():
    return os.path.join(
        theory_dir,
        'pm%i_laget_theory_set%i.txt' % (pm_set, data_set)
    )

output_file = theory_file()

o = open(output_file,'w')
# write header
o.write(header)

# some constants
dtr = np.pi/180.
#MeV
MP = 938.272
MN = 939.566
#MD = 1875.6127
MD = 1875.61
me = 0.51099 

fm2ub = 1.e4   #fm^2 to ub

# initalize arrays: there is a link in the current directory
# assume it is in the current directory
deut_data_dir = './'

# select what kind of calculation and set flags
#do_fsi = 0
#do_pwia = 1

# select linear interpolation
lin_interp = 1
# do not save the grid as binary file
save_grid = 0

# intitialize the code
#if use_binary:
#    LX.init_laget('./', do_fsi, lin_interp, save_grid, use_binary)
#else:
#    LX.init_laget('./', do_fsi, lin_interp)


#======================================DO PWIA=============================================================

# read the averaged kinematics file
f = get_file(avgkin_file('pwia'))

#Read Headers to be written to the output Xsec file
i_b = get_data(f, 'i_b')      #ith (Pm, th_nq) bin number
i_x = get_data(f, 'i_x')      #ith th_nq bin number
i_y = get_data(f, 'i_y')      #ith Pm bin number
xb  = get_data(f, 'xb')       #theta_nq central bin value [deg]
yb  = get_data(f, 'yb')       #Pmiss cental bin value [GeV]
cont = get_data(f, 'cont')    #2D bin content

#===============================================================================
#=========  Laget_Xsec_fp_1.f,  get_sigma_laget() Input Parameters =============
#===============================================================================
#units are in MeV and deg [convert to radians, as Laget code takes MeV and radians]

th_cm = get_data(f, 'th_pq_cm')*dtr     #center of mass in-plane angle between proton and q-vector [rad]
omega = get_data(f,  'omega')           #averaged calculated energy transfer [MeV]
q = get_data(f, 'q_lab')                #averaged magnitude of 3-momentum transfer [MeV]
Q2 = get_data(f, 'Q2_calc')             #averaged calculated 4-Momentum Transfer [MeV^2]
cphi = get_data(f, 'cos_phi')           #averaged Cos(phi_pq) : out-of-plane angle between proton and q-vector, determined from SIMC  vertex quantities 
sphi = get_data(f, 'sin_phi')
phi_pq_carlos = np.arccos(np.clip(cphi, -1.0, 1.0))                    #phi, out-of-plane angle between the proton and q-vector (angle between reaction-scattering plane) [rad]
phi_pq_right= np.atan2(sphi,cphi)
Ei = get_data(f, 'Ei')                  #averaged Incident Beam Energy [MeV]
GEp = get_data(f, 'GEp')                #Proton Electric Form Factor at the avg. kinematics
GMp = get_data(f, 'GMp')                #Proton Magnetic Form Factor at the avg. kinematics
sigMott = get_data(f, 'sigMott')        #Mott cross section at the avg. kinematics
Ksig_cc1 = get_data(f, 'Ksig_cc1')      #Kinematic Factor K * deForest cc1 cross section at the avg. kinematics

# DO PWIA
if use_binary:
    LX.init_laget('./', 0, lin_interp, save_grid, use_binary)
else:
    LX.init_laget('./', 0, lin_interp)

sigma = []
for i, e0_i in enumerate(Ei):
    # check kinematics
    #print("i = ",i," e0_i= ",e0_i, " Q2= ",Q2[i]," nu= ", omega[i]," th_cm=",th_cm[i], " phi=",phi[i])
    #print "p_rec = ", LX.p_recoil(Q2[i], omega[i], th_cm[i])
    
    if cont[i]==0.0:
        sig         = -1.
        GEp[i]      = -1.
        GMp[i]      = -1.
        sigMott[i]  = -1.
        Ksig_cc1[i] = -1.
        #sigma.append(sig)
    else:
        # calculate cross sections (fm^2 sr^-2 mev^-1)
        sig = LX.get_sigma_laget(e0_i,Q2[i],omega[i],th_cm[i],phi_pq_right[i]) * fm2ub    #convert to (ub sr^-2 MeV^-1)
        #sigma.append(sig)
        #print('sig = ',sig)
        # Write to File
    l = '%i %i %i %.6f %.6f %.12e  %.6f  %.6f  %.6f  %.6f\n'%(i_b[i], i_x[i], i_y[i], xb[i], yb[i], sig, GEp[i], GMp[i], sigMott[i], Ksig_cc1[i])
    o.write(l)
o.close()

#======================================DO FSI================================================

# read the averaged kinematics file, and output file to write cross sections
f = get_file(avgkin_file('fsi'))

fname = theory_file()
output_file = get_file(fname)

#Add Key to the output_file header
output_file.add_key('fsiXsec', 'f')
output_file.add_key('fsiGEp', 'f')
output_file.add_key('fsiGMp', 'f')
output_file.add_key('fsi_sigMott', 'f')
output_file.add_key('fsi_Ksig_cc1', 'f')
#Add additional avg. kinematics to output_file header
output_file.add_key('Ei_avg')     
output_file.add_key('th_pq_cm')   
output_file.add_key('omega_avg')  
output_file.add_key('q_avg')      
output_file.add_key('Q2_avg')     
output_file.add_key('cphi_pq_avg')
output_file.add_key('sphi_pq_avg')
output_file.add_key('phi_pq_carlos')
output_file.add_key('phi_pq_right')
output_file.add_key('pf_avg')     
output_file.add_key('pm_avg')     
output_file.add_key('kf_avg')     
output_file.add_key('the_avg')    
output_file.add_key('xbj_avg')    

# Default FSI-derived fields to -1. PWIA and FSI can contain different bin
# subsets, so unmatched PWIA bins must not retain the add_key() default of 0.
fsi_output_keys = [
    'fsiXsec', 'fsiGEp', 'fsiGMp', 'fsi_sigMott', 'fsi_Ksig_cc1',
    'Ei_avg', 'th_pq_cm', 'omega_avg', 'q_avg', 'Q2_avg',
    'cphi_pq_avg','sphi_pq_avg','phi_pq_carlos','phi_pq_right', 'pf_avg', 'pm_avg', 'kf_avg', 'the_avg', 'xbj_avg',
]
for row in output_file.data:
    for key in fsi_output_keys:
        row[key] = -1.0

# Match FSI rows to PWIA output rows using the shared 2D bin coordinates.
# The global i_b values and row counts can differ between model samples.
output_index_by_bin = {
    (int(row['i_x']), int(row['i_y'])): i
    for i, row in enumerate(output_file.data)
}
fsi_i_x = get_data(f, 'i_x')
fsi_i_y = get_data(f, 'i_y')
cont = get_data(f, 'cont')

#===============================================================================
#=========  Laget_Xsec_fp_1.f,  get_sigma_laget() Input Parameters =============
#===============================================================================
#units are in MeV and deg [convert to radians, as Laget code takes MeV and radians]
#All these quantities are averaged from SIMC or they have been calculated using the assumed masses
th_cm = get_data(f, 'th_pq_cm')*dtr     #center of mass in-plane angle between proton and q-vector [rad]
omega = get_data(f,  'omega')           #energy transfer [MeV]
q = get_data(f, 'q_lab')                #magnitude of 3-momentum transfer [MeV]
Q2 = get_data(f, 'Q2_calc')             #4-Momentum Transfer [MeV^2]
cphi = get_data(f, 'cos_phi')           #averaged Cos(phi_pq) : out-of-plane angle between proton and q-vector, determined from SIMC  vertex quantities 
sphi = get_data(f, 'sin_phi')
phi_pq_carlos = np.arccos(np.clip(cphi, -1.0, 1.0))                    #phi, out-of-plane angle between the proton and q-vector (angle between reaction-scattering plane) [rad]
phi_pq_right= np.atan2(sphi,cphi)                    #phi, out-of-plane angle between the proton and q-vector (angle between reaction-scattering plane) [rad]
Ei = get_data(f, 'Ei')                  #Incident Beam Energy [MeV]
GEp = get_data(f, 'GEp')                #Proton Electric Form Factor at the avg. kinematics
GMp = get_data(f, 'GMp')                #Proton Magnetic Form Factor at the avg. kinematics
sigMott = get_data(f, 'sigMott')        #Mott cross section at the avg. kinematics
Ksig_cc1 = get_data(f, 'Ksig_cc1')      #deForest cc1 cross section at the avg. kinematics
#Additional averaged kinematics to add
Pf = get_data(f, 'pf')                  #averaged proton final momentum
pm_avg = get_data(f, 'pm')              #averaged calculated missing momentum assuming deuteron proton and neutron mass
kf = get_data(f, 'kf')                  #avergaed final e- momentum
th_e = get_data(f, 'th_e')               #averaged in-plane electron angle
xbj = Q2 / (2.*MP*omega)                  #averaged x-Bkorker determined from calculated averaged Q2 and omega

# DO FSI
if use_binary:
    LX.init_laget('./', 1, lin_interp, save_grid, use_binary)
else:
    LX.init_laget('./', 1, lin_interp)

#sigma = []
for i, e0_i in enumerate(Ei):
    # check kinematics
    #print("i = ",i," e0_i= ",e0_i, " Q2= ",Q2[i]," nu= ", omega[i]," th_cm=",th_cm[i], " phi=",phi[i])
    #print "p_rec = ", LX.p_recoil(Q2[i], omega[i], th_cm[i])
    output_i = output_index_by_bin.get((int(fsi_i_x[i]), int(fsi_i_y[i])))
    if output_i is None:
        print('WARNING: skipping FSI-only bin (i_x=%i, i_y=%i)' %
              (fsi_i_x[i], fsi_i_y[i]))
        continue
    output_row = output_file.data[output_i]

    if cont[i] != 0.0:
        # calculate cross sections (returns Xsec in fm^2 sr^-2 MeV^-1 (SIMC is in microbarn.  1 ub = 1e-4 fm^2)
        sig = LX.get_sigma_laget(e0_i,Q2[i],omega[i],th_cm[i],phi_pq_right[i]) * fm2ub    #convert to (ub sr^-2 MeV^-1)
        output_row['fsiXsec'] = float("%.12e"%(sig))                    #ub sr^-2 MeV^-1
        output_row['fsiGEp'] = float("%.6f"%(GEp[i]))
        output_row['fsiGMp'] = float("%.6f"%(GMp[i]))
        output_row['fsi_sigMott'] = float("%.6f"%(sigMott[i]))           #ub / sr
        output_row['fsi_Ksig_cc1'] = float("%.6f"%(Ksig_cc1[i]))         # ub MeV^2 / sr^2
        #Write averaged kinematics quantities 
        output_row['Ei_avg']     = float("%.6f"%(Ei[i]))                    
        output_row['th_pq_cm'] = float("%.6f" % (th_cm[i] / dtr))     
        output_row['omega_avg']  = float("%.6f"%(omega[i]))     
        output_row['q_avg']      = float("%.6f"%(q[i]))
        output_row['Q2_avg']     = float("%.6f"%(Q2[i]))
        output_row['cphi_pq_avg']= float("%.6f"%(cphi[i]))
        output_row['sphi_pq_avg']= float("%.6f"%(sphi[i]))
        output_row['phi_pq_carlos'] = float("%.6f" % phi_pq_carlos[i])
        output_row['phi_pq_right'] = float("%.6f" % phi_pq_right[i])
        output_row['pf_avg']     = float("%.6f"%(Pf[i]))
        output_row['pm_avg']     = float("%.6f"%(pm_avg[i]))
        output_row['kf_avg']     = float("%.6f"%(kf[i]))
        output_row['the_avg']    = float("%.6f"%(th_e[i]))
        output_row['xbj_avg']    = float("%.6f"%(xbj[i]))
        
# Write to File
output_file.save(fname)
