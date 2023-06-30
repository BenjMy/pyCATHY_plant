To describe the presence and activity of the plant, the `plant_parm` file is used. This file contains the following input parameters and values:

- `NPLANT`: Number of plants in the model.
- `NPTYPE`: Number of plant types in the model.
- `ACANOFLAG`: Set to 0 to indicate if the projected canopy area is given in absolute value or relative to the nodal area.
- `PLANT_PRINT`: Set to 0 to specify the type of output.
- `NMETEODATA`: Set to 4 if the Leaf Area Index (LAI) is constant, or 5 if it is time-dependent and read from another input file.
- `RHOW`: Water density (constant value).
- `COATM`, `CCO2STAR`, and `CCO2ATM`: Concentrations of O2 in the atmosphere, reference CO2 concentration, and CO2 concentration in the environment, respectively.
- `SSTOMA`, `ASTOMA`, `TOLLNR`, `ITMAXNR`, `PSILEAF0`, `PSTEP`, and `NPMED`: Parameters for plant growth and water stress calculations.

Each row corresponds to a plant, with the first number indicating the node and the second indicating the plant type. The number of rows equals the value of NPLANT.

Additional parameters include:

- XMAXP, YMAXP, and ZMAXP: Semi-lengths describing the root extension along the x, y, and z axes, respectively.
- VRUX, VRUY, VRUZ, VRUX1, VRUY1, and VRUZ1: Unused parameters.
- LAI: Constant value of leaf area index.
- GXYLEM_MAX, C_GXYLEM, and D_GXYLEM: Maximum xylem conductance and coefficients for vulnerability curve.
- VCMAX, KC25, KO25, and COMP25: Parameters related to carbon fixation and compensation point.
- LA_MAX, LA_BETA, and LA_PSILMAX: Parameters for maximum marginal water use efficiency.
- PSILMAX: Threshold value to determine if the plant is functioning properly based on leaf water potential.
- LIMIT: Apparent quantum yield.
- HCANO: Canopy height.
- DATA_LEAF, ACANO, AXYLEM, GROOT_STAR, and GSOIL_STAR: Parameters related to canopy, xylem area, and root conductance.
- SALT_TOX: Unused parameter.
- NDATA_RDF, ZRDF, and RDFVAL: Parameters for root density function.

For columns from XMAXP to SALT_TOX, there will be as many columns as plant types, where the first column represents type 1 plants, the second column represents type 2 plants, and so on.

The second input file `plant_meteo` allows you to specify meteorological data for the model. It is created using the out_meteo.exe executable, and the required input files are:

- time.txt: Contains time steps at which meteorological parameters are available.
- Temperatura.txt: Contains temperature values corresponding to the time steps in time.txt.
- RH.txt: Contains relative humidity values corresponding to the time steps in time.txt.
- PAR.txt: Contains photosynthetically active radiation values corresponding to the time steps in time.txt.
- ZEN.txt: Contains zenith angle values corresponding to the time steps in time.txt.
- LAI.txt: Contains leaf area index values corresponding to the time steps in time.txt (read if NMETEOPLANT=5 in plant_parm file).

The `plant_growth` and `plant_salt` files are also used for some cases.

