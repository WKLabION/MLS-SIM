Code
\Supplementary Software\MLS-SIM Reconstruction

Prerequisites:
PC with Windows 10 operating system
At least 16GB RAM
MATLAB version R2021a or higher
CUDA Library compatible with MATLAB
A CUDA compatible graphic card and 11GB graphic memory

Installation guide:
Install all the software in the prerequisites in order. The typical install time on a "normal" desktop computer is about 30 minutes.

Expected output and runtime:
A super-resolution image reconstructed from the examplar raw image. The typical runtime is about 1 min with a NVIDIA V100 GPU.

main.m
This script reconstruct super-resolution image from raw image captured by MLS-SIM. This script could run directly.

example\
This folder contains an example raw image for reconstruction.

params\
This folder contains example PSFs and parameters of the imaging system which could be used for reconstruction of the example data.

subfunctions\
This folder contains subfunctions used in main.m.


