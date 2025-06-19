# MCXQI Simuloation Package
Simulation software for X-ray Quantum Imaging using Monte-Carlo method
## Introduction
MCXQI package is using ray tracing technique. As a ray tracing technique it is modular, thus, composed of different modules perform different functions.
The package include the following modules:
1. make_source - generation of rays starting from a chosen position in space. each ray has a coordinate and direction in 3D.
                 The module needs the following parameters: Position, Wavelength, Number of rays and other optional parameters.
2. SPDC        - simulates a non linear crystal with a certain dimentions, positoin and orientation in 3D. The module simulates the
                 transfomration of each incoming ray to two out going rays with randomly chosen directions and wavelengths. The module
                 makes sure the incoming energy and momentum is conserved, i.e. equal to the sum of out going energy and momentum.
3. design_object - uploads an 2d image and converts it to an object with variable absorption values according to the values in the image. 
