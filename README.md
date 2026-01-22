# BayesApp
BayesApp calculated the pair distance distribution function (PDDF or p(r) function) for SAXS or SANS data
version 2.1

## How to install/run

#### Run throught the web application 
BayesApp can be run through the webapp, which also offers a graphical user interface (GUI) through our [GenApp web application](https://somo.chem.utk.edu/bayesapp/).    
For developers: The source-code for the web-app is available on the dedicated [GitHub page](https://github.com/ehb54/GenApp-BayesApp).    

#### Run locally on your own computer, step 1: compile the fortran code
To run locally (no GUI), the Fortran code, `bift.f` needs to be compiled, e.g. with (linux/mac):
```
gfortran bift.f -march=native -O2 -o bift
```
on windows, the executable is assumed to have .exe extension, so, e.g.:
```
gfortran bift.f -O2 -o bift.exe
```
The name has to be bift for bayesapp.py to be able to find it
Requirements is a gfortran compiler. Compilation is, in our experience, straightforward on Linux (Ubuntu), but troublesome on MacOS/Windows. Unfortunately, we do not have resources to help with compilation.   

#### Run locally on your own computer, step 2: run the python program
When compilation is done, the program can be run as a standard `python` (python3) program: 
```
python bayesapp.py -f datafile.dat
````
with -f (datafile) being the only required input.

##### requirements
* `bayesapp_helpfunctions.py` and the executable `bift` or `bift.exe` (compiled version of `bift.f`, see step 1) must be in the same folder. 
* `python` requirements (can be installed with `pip`): `numpy`, `scipy` and `matplotlib`. 

#### options/flags
too see all options, run this command 
```
python bayesapp.py -h
````

## How to cite
If you use BayesApp (web version or locally), please cite our most recent publication:    
[Larsen and Pedersen, 2021](https://doi.org/10.1107/S1600576721006877)

Optionally, also cite the publication describing the core algorithm (Bayesian Indirect Fourier Transformation, BIFT): 
[Hansen, 2000] (https://doi.org/10.1107/S0021889800012930)

If you use or report the `number of good parameters` as a measure for the information content in SAXS/SANS data, please cite: 
[Vestergaard and Hansen, 2006](https://doi.org/10.1107/S0021889806035291)

If you rescale your errorbars using BayesApp, please cite: 
[Larsen and Pedersen, 2021](https://doi.org/10.1107/S1600576721006877)

## Developers/maintainers
BayesApp was originally written by Steen Hansen, University of Copenhagen.    
The program has been further developed by and is currently maintained by Andreas Haahr Larsen, University of Copenhagen.   
The GUI was developed in collaboration with Emre Brookes, University of Montana.    


