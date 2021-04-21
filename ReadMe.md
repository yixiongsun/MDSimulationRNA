# MDSimulationRNA
COMP 559 Final Project

## Description
This project is a simple Molecular dynamics simulation of RNA, coded in python and rendered with PyMol. Here is the youtube link for the video description:

## Installation
Python3 is required to run the simulation program. We recommend using a virtual environment for managing python packages.

To install the required packages, run 
```$ pip install -r requirements.txt```
This will install Biopython, Numpy and Taichi.

To render animations for the simulation, PyMol is required and a free trial can be downloaded here: https://pymol.org/2/

## Run the simulation
To run the simulation, we need RNA crystal structure files from the Protein Data Bank. There are a couple sample files found in the pdbs folder. You can find more online at https://www.rcsb.org/

The command to run the simulation is

```$ python run_simulation.py [pdb filename] [pdb chain] [time step] [total steps] [save step size] [gpu/cpu]```

where ```[]``` are arguments passed to the command. 

- ```[pdb filename]``` is the path to the pdb file.
- ```[pdb chain]``` is the desired chain in the crystal structures. These chains are specific to the structure and must be manually specified.
- ```[time step]``` is the time step of the simulation. Good time steps are around 0.01 which is around 0.5 fs.
- ```[total steps]``` is the total number of steps of the simulation.
- ```[save step size]``` defines at which time steps to save the positions of the particles. For example, 10 would mean saving the position of the simulation every 10 steps.
- ```[gpu/cpu]``` defines whether the simulation should be ran with the gpu or cpu

An example usage of the command

```$ python run_simulation.py pdbs/1ivs.cif C 0.01 2000 10 gpu```

## Contact
If you have any questions, please contact me at yixiong.sun@mail.mcgill.ca or create an issue in this repo.