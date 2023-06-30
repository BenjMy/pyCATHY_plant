# pyCATHY plant model DA

This documentation provides an overview of the project code and presents the encouraging results obtained through the implementation of Data Assimilation (DA) techniques to enhance the performance of the CATHY hydrological model solver and invert model parameters.

## Introduction
The objective of this project is to improve the accuracy and reliability of the CATHY model by integrating DA methodologies. In this implementation, we have utilized [pyCATHY](https://github.com/BenjMy/pycathy_wrapper), a Python wrapper for CATHY, to externally incorporate DA functionality and streamline file updates through Python routines.

**Successful Assimilation with pyCATHY**
DA using pyCATHY has shown promising results for: 
- ERT (Electrical Resistivity Tomography),
- soil water content, and tensiometers.
  
1 article based on the outcomes of this work (in prep.)

## Expanding Assimilation Scope
We are considering the incorporation of plant data observations such as **leaf water potential** or **stem flow** as additional assimilation observations. This expansion has the potential to provide valuable insights and contribute to a more comprehensive understanding of the hydrological model.


### The CATHY plant model

Source files are different from the [CATHY V1](https://bitbucket.org/cathy1_0/cathy/src/master/). 

- Work from Laura Busato - see Busato
- Work from Gabriele - see Manoli

### DA Procedure within CATHY plant model

Please refer to the accompanying resources for detailed instructions on running the code and implementing data assimilation using pyCATHY.

The procedure to prepare is logistically the same with a few changes:
- Input preparation
- Mapping operator

### Bugs

- See ?
- See 
