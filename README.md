# MCXQI Simulation Package
A Monte Carlo-based ray-tracing software for simulating X-ray Quantum Imaging\
\
Authors: Pakhshan Espoukeh*, Gabriel Biener*, Andrew Aquila**, Peter Schwander*

\*  University of Wisconsin – Milwaukee, Milwaukee, United States\
\** National Accelerator Laboratory, Menlo Park, United States

## Introduction
MCXQI is a modular Monte Carlo ray-tracing package which includes a set of dedicated modules which perform specific functions on the rays along their opticaL path.

The implementation of the software makes use of the UNIX pipeline mechanism. This allows for high flexibility in simulating a variety of optical setups. 

The package includes among many other modules:
1. _make_source_ - Generates rays emitted from an X-ray source of a given position, direction, size, divergence and bandwidth.\
                 The module needs the following parameters: Position, Wavelength, Number of rays and other optional parameters.
2. _SPDC_        - Simulates the process of Spontaneous Parametric Down-Conversion (SPDC) from a pump ray hitting a nonlinear crystal of a given size, position and orientation in space. The module simulates the\
                 stochastic conversion of a pump ray into two outgoing rays with directions and wavelengths conserving
                 energy and momentum taking account the cross section of SPDC.
3. _design_object_ - Simulates the effect of an object with spatially varying absorption from values specified by a two-dimensional image.
4. _propagate_to_ - Updates the position and path length of each ray based on the ray's direction and the point of intersection with any given plane in space.
5. _aperture_ - Simulates the effect of a circular aperture where rays can pass through inside the aperture and are blocked otherwise.
6. _print_rays_ -  Lists all ray information in text format.
7. _absorption_object_ - Simulates the effect of a round object with a given absorption inside the object and no absorption outside.

## Installation
Requirements: GCC compiler
1. Make a copy of the source code from the repository.
2. Compile all modules using make and the input file Makefile. Open a terminal and type on the command line\
   $ make all
<!-- 3. To compile individual modules open terminal and type:\
   (base) file_location ~ % $${\color{red}gcc~-lm~-Wall}$$ $${\color{red}module\underline{ }name.c}$$ $${\color{red}Required\underline{ }File1.o~Required\underline{ }File2.o}$$ $${\color{red}-o}$$ $${\color{red}Exec\underline{ }File\underline{ } Name}$$
   -->

## Execution
As an example, the simulation of an optical setup with an X-ray source hitting a non-linear crystal and observing the SPDC on a given plane would be executed as follows:

$ _make_source_ --position='[x,y,z]' --wavelength=w -n 10 | _SPDC_ --P0='[x0,y0,z0]' --P1='[x1,y1,z1]' --P2='[x2,y2,z2]' -m '[h,k,l]' -p 'angle' | _propagate_to_ --P0='[x0,y0,z0]' --P1='[x1,y1,z1]' --P2='[x2,y2,z2]' | _print_rays_

