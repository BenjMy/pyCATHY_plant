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
CATHY_PLANT_1D = CATHY(dirName='../examples/Mary_pyCATHY/',
                      prj_name='CATHY_PLANT_1D',
                      )
                      #version='plant') # + 'test')



CATHY_PLANT_1D.run_preprocessor(verbose=True)

    
    
                            





      
      


