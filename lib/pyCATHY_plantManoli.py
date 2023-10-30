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
from pyCATHY.importers import cathy_inputs as CTin

#%%
CATHY_PLANT_1D = CATHY(dirName='../examples/Mary_pyCATHY/',
                      prj_name='CATHY_PLANT_1D',
                      version='G. Manoli',
                      )
                      #version='plant') # + 'test')


#%%
CATHY_PLANT_1D.run_preprocessor(verbose=True)

#%%

# file2read = os.path.join(CATHY_PLANT_1D.workdir,
#                           CATHY_PLANT_1D.project_name,
#                           'input',
#                           'plant_parm'
#                           )
# # plant_parm_dict = CTin.read_plant_parm(file2read)
 
# # CATHY_PLANT_1D.read_inputs('parm')


# file2read = os.path.join(CATHY_PLANT_1D.workdir,
#                           CATHY_PLANT_1D.project_name,
#                           'input',
#                           'parm'
#                           )

# CTin.read_parm(file2read,
#                version='G. Manoli'
#                )

#%%
# NEW PARM file in v0.0.1_
# - TRAFLAG
# - VELREC
# - OMEGA is integer !
# -  IPRT VTKF NPRT (TIMPRT(I),I=1,NPRT) instead of IPRT VTKF (TIMPRT(I),I=1,NPRT)
# NEW DEM ?
CATHY_PLANT_1D.run_processor(verbose=True)

                            




#%%

plant_leaf_explanations = {
    'TIME': 'Time step of CATHY',
    'TRASP': 'Total transpiration calculated by PSIR. It is expressed in m/s, so to obtain the value in m/s, it should be divided by ACANO',
    'QTOT': 'Total transpiration calculated as the sum of nodal values QPLANT (as a check, TRASP should be equal to QTOT)',
    'QTOT_HR': 'Total hydraulic redistribution flow (negative QPLANT values, from plant to soil)',
    'QTOT_TRASP': 'Root water uptake flow (positive QPLANT values, from soil to plant). The total transpiration is the sum of hydraulic redistribution and root water uptake; these two components can be considered separately',
    'PSILEAF': 'Leaf water potential. If this variable reaches the minimum value indicated by PSILMAX, then the plant is not functioning correctly',
    'PSIR': 'Water potential at the base of the plant trunk',
    'PSILMEAN': 'Mean of PSILEAF over the previous 24 hours',
    'GSTOMA(15)': "Stomatal conductance. In this case, it prints the value at layer 15 since the 'big leaf' is divided into layers (30 in this case)",
    'FC': 'Carbon flux',
    'LASTOMA': 'λ in [3] and [2]'
}


# Create a dictionary from the list
plant_leaf_dict = dict(plant_leaf_explanations)

# Display the dictionary
print(plant_leaf_dict)

file2read = os.path.join(CATHY_PLANT_1D.workdir,
                          CATHY_PLANT_1D.project_name,
                          'output',
                          'plant_leaf'
                          )

plant_leaf_df = pd.read_csv(file2read,delim_whitespace=True,header=None)
plant_leaf_df.columns = plant_leaf_explanations.keys()

plant_leaf_df.set_index('TIME', inplace=True)

plant_leaf_df.plot(y='TRASP')

plant_leaf_df.plot(y='QTOT_HR')

plant_leaf_df.plot(y='PSILEAF')






      
      


