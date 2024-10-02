# Tumor-Growth
Tumor Growth 2023 internship. The code is largely inspired from "Learning hydrodynamic equations for active matter from particle simulations and experiments" PNAS, 2023, Vol. 120.

---
The three major steps in the learning framework are: 

1. Data interpolation
2. Spectral projection and smoothing 
3. PDE learning 

The codes provided need to be run sequentially to execute the above steps. The details are provided below. The starting point is the data provided in `data/defectors06_07_2023 13_3*.txt`. 

---
## 1. Data interpolation [MATLAB] 

- Run `coarse_grain_data.m` through MATLAB. This takes in the particle data from `data/quinke_particle_data_S2.mat` and saves the coarse-grained data at `data/coarse_grained_field(y,x,t).mat`. The generated image `data/coarse_grained_frame.png` shows the coarse-graining result at a particular frame. 

- With `Nx=100, Ny=50` (specified inside `coarse_grain_data.m`), the code should run for ~10 mins.

---
## 2. Spectral representation [Julia] 

- If you do not have Julia versions from 1.3 - 1.7, download it from: https://julialang.org/downloads/. Make sure the command `julia` works from the command line.

- Navigate to the project folder through the command line. 

- Saving spectral coefficients of the coarse-grained data: 
    - Run `julia save_representation_heights.jl`. This should run for a few mins.

- Saving smooth data and time-derivative after imposing cut-offs on the coefficients:
    - Run `julia save_smooth_data_heights.jl`. This should run for ~ 5 mins.

---
## 3. PDE Learning [Python 3]

We use Python >=3.6.

- Navigate to the project folder through the command line. 

- Installing required packages:

Conda[Recommended]: 
- If you do not have anaconda, download it from here: https://www.anaconda.com/products/individual. After the installation, make sure the command <conda list> works from your command line. 
- Create a new conda environment from the provided .yml file: run 
    `conda env create --name pdel --file=environment.yml`

- Activate the environment with `conda activate pdel`. 

OR 

PIP:
- If you already have a Python3 installation that you would like to use, do `pip3 install -r requirements.txt`.

### **Learning density equation:**
- Run `python3 learn_density_equation.py`. This should run for a few mins. 
- After the script runs, the folder `data/PDElearn/density_equation/` will have the learning results. 
- The file `pde_stability.txt` lists all the learned PDEs and the files `pde_stability_*.txt` will have the individual PDEs. 

### **Learning velocity equation:**  
- Do `python3 learn_velocity_equation.py`. This should run for a few mins.
- After the script runs, the folder `data/PDElearn/velocity_equation/` would have the learning results. 
- The file `pde_stability.txt` lists all the learned PDEs and the files `pde_stability_*.txt` will have the individual PDEs. 
