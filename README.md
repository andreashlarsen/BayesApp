# BayesApp
BayesApp calculated the pair distance distribution function (PDDF or p(r) function) for SAXS or SANS data
version 2.0

## How to install/run

#### Run throught the web application 
BayesApp can be run through the webapp, which also offers a graphical user interface (GUI) through our [GenApp web application](https://somo.chem.utk.edu/bayesapp/).    
For developers: The source-code for the web-app is available on the dedicated [GitHub page](https://github.com/ehb54/GenApp-BayesApp).    

#### Run locally on your own computer, step 1: compile the fortran code
To run locally (no GUI), the Fortran code, `bift.f` needs to be compiled, e.g with:
```
gfortran bift.f -march=native -O2 -o bift
```
The name has to be bift for bayesapp.py to be able to find it
Requirements is a gfortran compiler. Compilation is, in our experience, straightforward on Linux (Ubuntu), but troublesome on MacOS/Windows. Unfortunately, we do not have resources to help with compilation.   

#### Run locally on your own computer, step 2: run the python program
When compilation is done, the program can be run as a standard `python` (python3) program: 
```
python bayesapp.py -f datafile.dat
````
with -f (datafile) being the only required input. 
`bayesapp_helpfunctions.py` and the executable `bift` (compiled version of `bift.f`, see step 1) should be in the same folder. 
`python` requirements (can be installed with `pip`) are standard libraries of scientific computing `numpy` and `scipy` as well as `matplotlib` for plotting. 

too see all options, run the command: 
```
python bayesapp.py -h
````

## How to cite
if you use BayesApp (web version or locally), please cite the most recent publication: 

optionally also cite the publication describing the core algorithm (Bayesian Indirect Fourier Transformation, BIFT): 
[Hansen, 2000] (https://doi.org/10.1107/S0021889800012930)

if you use or report the `number of good parameters` as a measure for the information content in SAXS/SANS data, please cite: 
[Vestergaard and Hansen, 2006](https://doi.org/10.1107/S0021889806035291)

if you rescale your errorbars using BayesApp, please cite: 
[Larsen and Pedersen](https://doi.org/10.1107/S1600576721006877)

## Developers/maintainers
BayesApp was originally written by Steen Hansen, University of Copenhagen.    
The program has been further developed by and is currently maintained by Andreas Haahr Larsen, University of Copenhagen.   
The GUI was developed in collaboration with Emre Brookes, University of Montana.    


