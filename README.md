# Project5

In this project we simulate the time evolution of a wavefunction in a double-slit setup. 

To compile and link Main.cpp: 

g++  Main.cpp -larmadillo -o Main.exe

to run: 

./Main.exe h Dt T x_c sig_x p_x y_c sig_y p_y v0 slits

h = length-step, type double. 
Dt = time-step, type double. 
T = total time-period, type double. 
x_c = center in x-direction, type double. 
sig_x = width in x-direction, type double. 
p_x = momentum in x-direction, type double. 
y_c = center in y-direction, type double. 
sig_y = width in y-direction, type double. 
p_y = momentum in y-direction, type double. 
v0 = magnitude of the potential at the wall, type double. 
slits = number of slits, type int. Must be 1, 2 or 3. (For no wall, set v0=0.)


The Python codes has been run in Spyder using the "Run" command. To run wall.py, first run Main.cpp with lines 209 and 210 uncomented. To run the other Python codes, first run Main.cpp with appropriate input. An armadillo cube with the time-evolving wavefunction will be saved as Us.bin. Change the name of this file so it matches the names in the Python files. 
