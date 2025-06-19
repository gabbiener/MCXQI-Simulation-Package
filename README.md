# MCXQI Simuloation Package
Simulation software for X-ray Quantum Imaging using Monte-Carlo method
## Introduction
MCXQI package is using ray tracing technique. As a ray tracing technique it is modular, thus, composed of different modules perform different functions.
The package include the following modules:
1. make_source - generation of rays starting from a chosen position in space. each ray has a coordinate and direction in 3D.\
                 The module needs the following parameters: Position, Wavelength, Number of rays and other optional parameters.
2. SPDC        - simulates a non linear crystal with a certain dimentions, positoin and orientation in 3D. The module simulates the\
                 transfomration of each incoming ray to two out going rays with randomly chosen directions and wavelengths. The module\
                 makes sure the incoming energy and momentum is conserved, i.e. equal to the sum of out going energy and momentum.
3. design_object - uploads an 2d image and converts it to an object with variable absorption values according to the values in the image.
4. propagate_to - changes each ray position based on each ray direction and a distance between the original point and the destination point.
5. aperture - simulates a circular openning where rays can go through and outside the opening the rays stop from propagating.
6. print_rays - outputs the ray information into a text file.
7. absorption_object - a circular object with a certain absorption value inside the circle and absorption 1 outside of the circle.

## Installation
Requirements: C compiler and makefile software.
1. copy the source code.
2. compile all modules using make and the input file Makefile. Open a terminal and type\
   (base) file_location ~ % $${\color{red}make~all}$$
3. To compile individual modules open terminal and type:\
   (base) file_location ~ % $${\color{red}gcc~-lm~-Wall}$$ $${\color{red}module\underline{ }name.c}$$ $${\color{red}Required\underline{ }File1.o~Required\underline{ }File2.o}$$ $${\color{red}-o}$$ $${\color{red}Exec\underline{ }File\underline{ } Name}$$

## Running
To simulate we use a pipeline method of linux and run an entire optical setup. For example running a setup were X-ray light is reflected by a non-linear crystal and travels a certain distance using the terminal will look like:\
(base) file_location ~ % $${\color{red}\text{./make\\_source --position='[x,y,z]' --wavelength=w -n 10|./SPDC --P0='[x0,y0,z0]' --P1='[x1,y1,z1]' --P2='[x2,y2,z2]' -m '[h,k,l]'}}$$ $${\color{red}\text{-p 'angle'|./propagate\\_to --P0='[x0,y0,z0]' --P1='[x1,y1,z1]' --P2='[x2,y2,z2]'|./print\\_rays}}$$
