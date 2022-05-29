# Modelling & Visualization Exam 2022 - Three Chemical Fields
See the attached PDF with the corresponding password. The solutions here are working correct code for the exam, including
the notorious correlation function calculation. See [here](https://www.youtube.com/watch?v=tjKbJIN0Coo) for a video of it in action in 3D.

![alt text](https://github.com/callous4567/UoE-Projects/blob/master/SimAndVis/MVis_2022_Exam/Julia%20SharedArrays/sweep_1160.png)

I've included three options. These include my default Python solution (2D with Numba, as standard) and also two Julia
options- one with multiprocessing for plots and simulation in 3D (using SharedArrays) and one that just works single-threaded. The Julia
renditions tend to be extremely fast with simulation (around 10-30 times faster) but questionable with plotting (I've only tested in 3D, and
Makie meshscatters are extremely inefficient- indeed most of the multiprocessing is put into handling plotting.) 

## Exam
Usual Python/Numba 2D Exam Solutions

## Julia
Single-cored. Can not recommend- plotting is extremely slow with Makie. 

## Julia SharedArrays  
Has multiprocessing for plotting + simulation (recommend using a 4:1 split of plotting:simulation.)
