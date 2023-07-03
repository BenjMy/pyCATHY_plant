# L. Busato example description

**The exemple is thourougly described in L. Busato PhD thesis manuscript: Non-invasive monitoring and numerical modeling of the Soil-Plant continuum (2017).** and accessible at: http://paduaresearch.cab.unipd.it/10255/1/Laura_Busato_tesi.pdf 

## Objective

This is a synthetic case study in order to obtain a dataset suitable to test the reconstruction of plant transpiration. 
This synthetic case study consists of two hydrological models: the former is an infiltration model, while the latter describes the activity of a plant. 
These two models are identical, except for the presence of the tree.


## 1st model
The first model is aimed at representing the variation of spatial and temporal soil moisture content patterns as a consequence of the sole infiltration of irrigation water.

### Mesh

- 3d finite element mesh representing the investigated domain, over which the system of equation describing both surface and subsurface processes is solved. 
In our case, it consists of 9282 nodes and 52500 elements, subdivided into 25 horizontal layers with thickness increasing with depth (i.e. 16 layers 0.1 m thick, 7 layers 0.2 m thick, and 2 layers 0.5 m thick),
- Total depth is equal to 4.0 m, while length and width are both equal to 6.4 m. In spite of the dimensions of our mesh.
- We are actually interested only on a smaller subdomain (called “inner domain”), assumed to recall the volume investigated during the ERT surveys. In particular, it is defined by a
superficial square with side equal to 1.2 m (red polygon in Fig. 6.1) and extended till a depth of 1.2 m. Nevertheless, modeling a bigger domain is mandatory, in order to avoid the influence of boundary effects on the inner zone.

### Soil

- We characterized the soil by means of the parameters provided by Leij et al. [1996];
- We opted for a loamy sand soil, homogeneous and isotropic; 
- Te water retention curve is described thanks to van Genuchten and Nielsen [1985]. 
- The water table is located 3.0 m below ground level. 


### Boundary conditions

- On all nodes of the lateral boundary of the vadose zone (i.e. above the water table), as well as on the bottom face of the mesh, we imposed a Neumann condition
of no flow, while
- On the boundary of the saturated zone, pressure head is assumed to linearly increase with depth from a value of ψ = 0 m in correspondence of the water table (i.e. Dirichlet boundary condition) (Fig. 6.1). 
- Irrigation is imposed as an atmospheric boundary condition and takes place on a small portion of the upper surface (i.e. 39 nodes covering an area of about 0.18 m 2 , Fig. 6.1), so as to resemble
the drip irrigation occurring at the field site. We assumed an irrigation of 4 l h −1 (i.e. 22.1 mm h −1 ) for 5 h d −1 taking place for 15 days, then suspended for 5 days, and then performed for other 10 days, while rain precipitation is not considered.
- Water ponding is neglected. 
- The initial pressure head distribution (i.e. the initial condition) is nonuniform and has been computed thanks to another infiltration model, where irrigation is simulated for 15 days. This expedient is necessary since the infiltration and the plant models need to share the same initial conditions in order to be comparable, and derives from the impossibility to “activate” the plant
after a certain amount of time from the beginning of the simulation. 



- Therefore, the total simulation time is equal to 45 days. 
- Rime adaptation is intentionally avoided, in order to maintain the time step size constant (and equal to 120 s) throughout the whole simulation.
 
## 2nd model

- The plant model is identical to the infiltration model except for the presence of the tree. Therefore, mesh, irrigation schedule, boundary conditions, initial conditions, and soil parameters are those described before.

### Plant parameters

- The leaf area index is assumed constant over time, i.e. the plant does not grow nor the canopy is trimmed. The parameters can be either measured, i.e. from real orange trees,
assumed, or taken from literature; 
- The root system covers an area of 1.2 × 1.2 m 2 centered on the tree trunk and reaches a maximum depth of 0.4 m below groundlevel, in accordance with Kotur and Keshava Murthy [1998]. For the root length density vertical profile we assumed an exponential distribution [Volpe et al., 2013];
- As already described above, plant transpiration is mainly driven by external stressors (i.e. atmosphere and weather conditions). Therefore, also for this synthetic model, it is necessary to take into account meteorological information to determine the atmospheric forcing. In this case we took advantage of a real dataset measured from the meteorological station at the Bulgherano field site (see chap. 3), which consisted of photosynthetically active radiation, relative humidity, and temperature, all measured at a height of 4 m, except for PAR, which was measured at 8 m;
- The acquisitions took place on a hourly basis from 25 th September 2013 to 24 th October 2013;
- On the basis of these data, we considered precipitation negligible.


# G. Manoli Busato example description

**The exemple is thourougly described in G. Manoli PhD thesis manuscript: Contribution to modeling of soil-plant-atmosphere interactions and coupled hydro-geophysical data assimilation (2013).** and accessible at: https://hdl.handle.net/11577/3423556

RWU: Root Water Uptake
ABL: Atmospheric boundary layer 

- The study site (Fig. 3.1) is a 21 ha field located at the southern margins of the Lagoon
of Venice, North-East of Italy
- The field was cropped with maize (Zea mais L.)

## Objective
The model is applied to the field site to understand the impact of land elevation, soil heterogeneities, and seawater contamination on land productivity. 


### Model calibration
The model is calibrated on a single plant and then applied at the plot scale to obtain a sort of upscaled version of the plant model.
Our interest lies in the simulation of the **long term crop productivity**, the temporal dynamics of salt concentration in the vadose zone is neglected and assumed equal to the measured soil salinity 

### Mesh

- A cubic portion of soil of dimensions 5×5×5m (subsequently called the plot) is discretized in the x and y-directions with a 0.2 m spacing (fine grid) allowing a detailed description of
132 plants, i.e. the typical number of plants seeded in a 5×5 m plot.
- The upscaled model is then obtained by calibration of the model parameters on the same plot discretized with a coarse grid (2.5m spacing in x and y) where a single plant is used to represent the behavior of the whole fine-scale plot.

### Boundary conditions

- Field data (rainfall, temperature, relative humidity, radiation) are used to calculate the atmospheric forcing.
- Input evaporation is considered as a potential rate, and actual evaporation is evaluated based on system state condition allowing the switching between Neumann and Dirichlet boundary
conditions [Camporese et al., 2010].
- Water table levels are specified for the Northern and Southern boundaries of the model domain according to the observed water level in the irrigation channels. 
- No flow boundary conditions are set on the other edges of the domain. 



