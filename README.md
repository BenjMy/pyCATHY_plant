# pyCATHY plant model DA 🌱💧

This documentation provides an overview of the project code and presents the encouraging results obtained through the implementation of Data Assimilation (DA) techniques to enhance the performance of the CATHY hydrological model solver and invert model parameters.

## Introduction
The objective of this project is to improve the accuracy and reliability of the CATHY model by integrating DA methodologies. In this implementation, we have utilized [pyCATHY](https://github.com/BenjMy/pycathy_wrapper), a Python wrapper for CATHY, to externally incorporate DA functionality and streamline file updates through Python routines.

**Successful Assimilation with pyCATHY**
DA using pyCATHY has shown promising results for: 
- ERT (Electrical Resistivity Tomography),
- soil water content, and tensiometers.
  
1 article based on the outcomes of this work (in prep.)

## Expanding Assimilation Scope 🪝

We are considering the incorporation of plant data observations such as **leaf water potential** or **stem flow** as additional assimilation observations. This expansion has the potential to provide valuable insights and contribute to a more comprehensive understanding of the hydrological model.


### The CATHY plant model

While CATHY V1 uses the Feddes 1D RWU sink term approach, **CATHY plant model** uses a modeling framework that combines a 1D description of stem water flow, leaf-level photosynthesis, and transpiration with a **3D representation of soil-root water exchanges**.
- [plant inputs parameters](plant_inputs)
- [plant outputs files](plant_outputs)

Source files are different from the [CATHY V1](https://bitbucket.org/cathy1_0/cathy/src/master/). 
- [source files](/src/)


- Work from Laura Busato - see [Busato](/examples/Busato/)
- Work from Gabriele - see [Manoli](/examples/1_Gabriele_Piante_NON_modificato/)

### DA Procedure within CATHY plant model

Please refer to the the [online doc](https://benjmy.github.io/pycathy_wrapper/) for detailed instructions on running the code and implementing data assimilation using pyCATHY.

The procedure to prepare is logistically the same with a few changes:
- Input preparation 
- Mapping operator (to write)

## Bugs 🐛

- See [Compilation issue](https://github.com/BenjMy/pyCATHY_plant/issues/1#issue-1782683943)
- See [Processor issue]()

## Authors 

- L. Busato
- B. Mary


