#!/usr/bin python
import numpy as np
from matplotlib.patches import Ellipse
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import copy
from matplotlib.colors import SymLogNorm
from math import copysign

#import functions defining the VP rheologies
from VP_rheology_functions import *

# import settings
from VP_rheology_settings import *

# Using LaTeX in figures
plt.rc('text', usetex=True)
plt.rc('font', family='sans')

# create fake deformation data
# data = create_data(random=False,i=2e-5,j=0,plot=False,sym=True,s=201)
# data = create_data(random=True,i=2e-1,j=0,plot=False,sym=True,s=201)
data = create_data(random=True,i=7e-4,j=0,plot=False,sym=True,s=204)

# compute simpler and additionnal variables
comp_sim_sr(data)

data['rheos']=[]

######################
# CHOICE OF RHEOLOGIES
######################

# Examples of choice of rheologies
# NFR: normal flow rule
# NNFR: non normal flow rule
# MC: Mohr-Coulomb
# Ell. : Elliptical/Ellipse

# Ellipse NFR 'name':{'rheo_t':'ell', 'e':2.0, 'kt':0.5},
# Ellipse NFR with rotation 'name':{'rheo_t':'ell_rot', 'e':2.0, 'kt':0.5},
# Ellipse NNFR v2: 'e2r4t':{'rheo_t':'ellt', 'e':2.0, 'efr':4.0, 'kt':0.00},
# Ellipse NNFR: 'e2r4':{'rheo_t':'ell', 'e':2.0, 'efr':4.0, 'kt':0.00},
# MC yield curve shear only: 'mc.7s':{'rheo_t':'mcs', 'mu':0.7, 'kt':0.05},
# MC yield curve Ell. plastic potential: 'mc.7e4':{'rheo_t':'mce', 'mu':0.7, 'kt':0.05, 'e':4.0},
# MC yield curve Ell. plastic potential v2:'mc.7eG4':{'rheo_t':'mceG', 'mu':0.7, 'kt':0.05, 'e':4.0},
# Teadrop yield curve NFR: 'td':{'rheo_t':'td', 'kt':0.05},
# Parabolic lens yield curve NFR: 'pl':{'rheo_t':'pl', 'kt':0.05},
# MC yield curve, TD plastic potential: 'mc7td':{'rheo_t':'mctd', 'mu':0.7, 'kt':0.05},
# MC yield curve, PL plastic potential: 'mc7pl':{'rheo_t':'mcpl', 'mu':0.7, 'kt':0.05},
# Ell. yield curve, TD plastic potential: 'e2td':{'rheo_t':'etd', 'kt':0.05, 'e':2.},
# Ell. yield curve, PL plastic potential: 'e2pl':{'rheo_t':'epl', 'kt':0.05, 'e':2.},
# mu(I) rheology with dilatation: 'muID':{'rheo_t':'muID'},

#same kt for all rheologies, write a different one to override
kt = 0.05

rheo_to_viz = {
    # 'ell2':{'rheo_t':'ell', 'e':2.0, 'kt':kt},
    # 'ell2':{'rheo_t':'ell', 'e':1.0, 'kt':kt},
    # 'ell1p4':{'rheo_t':'ell', 'e':1.4, 'kt':kt},
    # 'ell4':{'rheo_t':'ell', 'e':4.0, 'kt':kt, 'plot_inv':True},
    # 'e2r4t':{'rheo_t':'ellt', 'e':1.0, 'efr':2.0, 'kt':kt},
    # 'e2r4':{'rheo_t':'ell', 'e':2.0, 'efr':4.0, 'kt':kt},
    # 'e2r1':{'rheo_t':'ell', 'e':2.0, 'efr':1.0, 'kt':kt},
    'mc.7s':{'rheo_t':'mcs', 'mu':0.7, 'kt':kt},
    'mc.7e1p4':{'rheo_t':'mce', 'mu':0.7, 'kt':kt, 'e':1.4, 'mu_c':4},
    'mc.7e1p4G':{'rheo_t':'mceG', 'mu':0.7, 'kt':kt, 'e':1.4, 'mu_c':4},
    # 'td':{'rheo_t':'td', 'kt':kt},
    # 'pl':{'rheo_t':'pl', 'kt':kt},
    'mc7td':{'rheo_t':'mctd', 'mu':0.7, 'kt':kt, 'mu_c':4},
    'mc7pl':{'rheo_t':'mcpl', 'mu':0.7, 'kt':kt, 'mu_c':4},
    # 'e2td':{'rheo_t':'etd', 'kt':kt, 'e':2.},
    # 'e2pl':{'rheo_t':'epl', 'kt':kt, 'e':2.},
    # 'muID':{'rheo_t':d'muID', 'db':False},
}

#######################
# COMPUTING VISCOSITIES
#######################

compute_visc(data=data, rheos=rheo_to_viz)

####################
# COMPUTING STRESSES
####################

compute_stress(data=data)

# for rheo_n in data['rheos']:
    # print(data[rheo_n].keys())

###############
# RUNNING TESTS
###############

compute_tests(data=data)

###################
# PLOTTING STRESSES
###################

plot_stress(data=data)

plot_FR(data=data)

plot_prAng(data=data)

###########
# PLOT SHOW
###########

plt.show()