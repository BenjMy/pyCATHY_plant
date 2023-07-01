#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 11:17:58 2023

@author: ben

The three-dimensional finite element mesh (6.4 × 6.4 × 4.0 m^3) consists of 9282 nodes 
and 52500 tetrahedral elements, divided into 25 horizontal layers with increasing 
thickness with depth (from 0.1 m to 0.5 m). 

However, our focus is only on an internal domain, with dimensions equal to those of
the volume investigated by ERT (1.2 × 1.2 × 1.2 m^3). 

The material, homogeneous and isotropic throughout the mesh, is described 
according to the parameters indicated in [6] for a sandy loam material,
 both in saturated and unsaturated conditions. 
 
The water table is located at a depth of 3.0 m. 

The boundary conditions consist of zero flux at the nodes above the water table (Neumann condition) 
and increasing pressure load with depth in the saturated zone (Dirichlet condition). 

Irrigation is imposed as an atmospheric boundary condition on a small portion of
 the surface face, with a flow rate of 4 lh^-1 for 5 hd^-1 for 15 days, then suspended for 5 days,
 and finally resumed for 10 days. 
 
Precipitation is considered negligible. The distribution of the pressure load for 
the initial conditions is obtained through a simulation of irrigation only for a duration of 15 days. 
The total simulation time is therefore 15 days, and the time step is kept constant at 120 s.

"""

#%% Compilation

# https://howtoinstall.co/en/m2c
#https://lists.mcs.anl.gov/pipermail/petsc-users/2010-November/007299.html

#%%
import pyCATHY
from pyCATHY import CATHY
import pandas as pd
import numpy as np
import os 

import utils_Busato
#%%
Busato_a_TREE = CATHY(dirName='../examples/Mary_pyCATHY/',
                      prj_name='49_2016_04_28_a_TREE') # + 'test')


#%%
Busato_a_TREE.run_preprocessor(verbose=True)

#%%
Busato_a_TREE.run_processor(verbose=True)


# #%%
# grid3d = Busato_a_TREE.read_outputs('grid3d')
# nnod3d = int(grid3d['nnod3'])

# x,y,z = [grid3d['mesh3d_nodes'][:,0],
#          grid3d['mesh3d_nodes'][:,1],
#          grid3d['mesh3d_nodes'][:,2],
#          ]
# nn = np.arange(1,len(x),1)



# areanodo = np.loadtxt(os.path.join(Busato_a_TREE.workdir, 
#                                    Busato_a_TREE.project_name, 
#                                    'output', 
#                                    'areanod'
#                                    )
#                       )
# nn2d=areanodo[:,0]
# areanod2d=areanodo[:,1]

# # areanod2d = area_e_nodo[:,2]



# # The three-dimensional finite element mesh (6.4 × 6.4 × 4.0 m^3) consists of 9282 nodes 
# # and 52500 tetrahedral elements, divided into 25 horizontal layers with increasing 
# # thickness with depth (from 0.1 m to 0.5 m). 

# #%%

# # area_e_nodo = np.loadtxt('./pyCATHY/area_e_nodo')

# # areanod2d = area_e_nodo[:,2]
# # nn2d = area_e_nodo[:,1]


# inputs_atmbc = {
#                 'xmin': 2.8,
#                 'xmax': 3.2,
#                 'ymin': 2.8,
#                 'ymax': 3.4,
#                 'flux': 1.11111E-06, #(in m^3/s, [L^3/T])
#                 'giorni_irrigazione': 15,
#                 'giorni_stop': 5,
#                 'giorni_irrigazione_2': 10,               
#                 }



# #%% Update atmbc

# # Irrigation is imposed as an atmospheric boundary condition on a small portion of
# # the surface face, with a flow rate of 4 lh^-1 for 5 hd^-1 for 15 days, then suspended for 5 days,
# # and finally resumed for 10 days. 

# HSPATM=0
# IETO=1
# zero=0.0
# time=25200


# #%% Case with suspension

# utils_Busato.atmbc_with_suspension(HSPATM,IETO,
#                                   zero,
#                                   inputs_atmbc,
#                                   nnod3d,nn,nn2d,
#                                   areanod2d,
#                                   x,y,z,
#                                   time,
#                           )

                
# #%% Update input plant
# # - Plant parm
# # - Plant meteo
# # - Plant growth
# # - Plant salt

# #%% Plant meteo


# time = np.loadtxt('./pyCATHY/create_input/time.txt') #contains the time steps at which the other meteorological parameter values are available 
# temp = np.loadtxt('./pyCATHY/create_input/Temperatura.txt') #contains the temperature values at the corresponding time steps in time.txt.
# rh = np.loadtxt('./pyCATHY/create_input/RH.txt') #contains the relative humidity values at the corresponding time steps in time.txt.
# par = np.loadtxt('./pyCATHY/create_input/PAR.txt') #contains the photosynthetically active radiation values at the corresponding time steps in time.txt.
# zen = np.loadtxt('./pyCATHY/create_input/ZEN.txt') #contains the zenith angle values at the corresponding time steps in the file.
# lai = np.loadtxt('./pyCATHY/create_input/LAI.txt') #contains the leaf area index values at the corresponding time steps in the file time.txt (read if NMETEOPLANT=5 in the file plant_parm).




# # Write plant_meteo file
# with open('plant_meteo', 'w') as file:
#     for i in range(len(time)):
#         file.write(f'{time[i]} TIME\n')
#         file.write(f'{temp[i]} {rh[i]} {par[i]} {zen[i]} {lai[i]}\n')


    
    
                            





      
      


