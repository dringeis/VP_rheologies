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
data = create_data(random=True,i=1e-1,j=0,plot=False,sym=True,s=201)

# compute simpler and additionnal variables
comp_sim_sr(data)

data['rheos']=[]

######################
# CHOICE OF RHEOLOGIES
######################

# Choice of rheologies
# Ellipse 'name':{'rheo_t':'ell', 'e':2.0, 'kt':0.5},
# Ellipse with rotation 'name':{'rheo_t':'ell_rot', 'e':2.0, 'kt':0.5},

rheo_to_viz = {
    'ell2':{'rheo_t':'ell', 'e':2.0, 'kt':0.05},
    'ell4':{'rheo_t':'ell', 'e':4.0, 'kt':0.05},
    'elr2':{'rheo_t':'ell', 'e':2.0, 'efr':4.0, 'kt':0.05},
}

#######################
# COMPUTING VISCOSITIES
#######################

compute_visc(data=data, rheos=rheo_to_viz)

####################
# COMPUTING STRESSES
####################

compute_stress(data=data)

for rheo_n in data['rheos']:
    print(data[rheo_n].keys())

###################
# PLOTTING STRESSES
###################

plot_stress(data=data)

#################################
# TEST FOR THE EFFECT OF ROTATION
#################################

# if ell and ell_rot :

#     sI_diff=sI_ell_rot-sI_ell
#     sII_diff=sII_ell_rot-sII_ell

#     print('sI_diff mean', np.nanmean(sI_diff))
#     print('sII_diff mean', np.nanmean(sII_diff))

#     print('sI_diff median', np.nanmedian(sI_diff))
#     print('sII_diff median', np.nanmedian(sII_diff))

#     print('sI std diff ', np.nanstd(sI_ell_rot)-np.nanstd(sI_ell))
#     print('sII std diff ', np.nanstd(sII_ell_rot)-np.nanstd(sII_ell))

#     print('sI_diff mean sI>-0.5', np.nanmean(sI_diff[sI_ell>-0.5]))
#     print('sII_diff mean sI<-0.5', np.nanmean(sII_diff[sI_ell<-0.5]))

#     plt.figure('diff sI')
#     plt.pcolormesh(sI_diff)
#     plt.axis('equal')
#     plt.colorbar()

#     plt.figure('diff sII')
#     plt.pcolormesh(sII_diff)
#     plt.axis('equal')
#     plt.colorbar()

#     A1=[sI_ell.flatten(),sI_ell_rot.flatten()]
#     B1=[sII_ell.flatten(),sII_ell_rot.flatten()]

#     plt.figure('diff yield curve')
#     plt.plot(A1,B1)
#     plt.axis('equal')

#     A2=[sI_ell.flatten(),sII_ell.flatten()]
#     B2=[sI_ell_rot.flatten(),sII_ell_rot.flatten()]

#     plt.figure('arrows diff yield curve')
#     for c in range(len(A2[0])):
#         # print(c)
#         # print(A2[0][c])
#         plt.arrow(A2[0][c], A2[1][c], B2[0][c] - A2[0][c], B2[1][c] - A2[1][c],head_width=0.01, length_includes_head=True,ec='k',fc='r')
#     plt.axis('equal')

###########
# PLOT SHOW
###########

plt.show()