# ParticleTracker
Multi-frame predictive particle-tracking algorithm for 2-D systems implemented using MATLAB. This code was written by Nicholas T. Ouellette (Stanford) and Douglas H. Kelley (University of Rochester) and is available in its original form here: https://web.stanford.edu/~nto/software_tracking.shtml

This particle-tracking algorithm was developed for Lagrangian particle tracking for flow visualization, but is broadly applicable to a variety of systems and has been used to follow the motion of animals, insects, and sediments. Particle trajectories are also processed by convolution with a smoothing and differentiating kernel to calculate velocities and accelerations (Mordant et al., 2004; Ouellette, 2006).

Since particle tracking experiments often produce a large quantity of images, this repository also includes scripts for quickly reading and viewing large uncompressed TIFF stacks or sequences of TIFF files. This TIFF reading and viewing code is adapted from Joseph M. Stujenske's *Matlab_FastTiffReadWrite* repository: https://github.com/jmstujenske/Matlab_FastTiffReadWrite

## Tutorial
The *TrackExample.m* file is a script that walks the user through basic functions of the predictive particle-tracking algorithm using an example movie of particles in a 2-D flow. 

The *PlotTracksExample.m* file is a script that goes through additional data processing, and particle and track visualization functionalities of code in this repository, including the *TiffTrackViewer* tool. This tool is a viewer with intuitive controls for playing TIFFs as a movie with overlaid particle tracks, as well as visualizing individual track statistics and trajectories.

![output](particle_track_mov.gif)

## Citation
This particle tracking code is based on the following publications: 

Ouellette, N.T., Xu, H. & Bodenschatz, E. "A quantitative study of three-dimensional Lagrangian particle tracking algorithms," Exp Fluids 40, 301-–313 (2006). https://doi.org/10.1007/s00348-005-0068-7

Kelley, D. H. and Ouellette, N. T. "Using particle tracking to measure flow instabilities in an undergraduate laboratory experiment," Am. J. Phys. 79, 267--273 (2011). https://doi.org/10.1119/1.3536647

## References
Mordant, N., Crawford, A.M., and Bodenschatz, E. "Experimental Lagrangian acceleration probability density function measurement." Physica D 193, 245--251 (2004). https://doi.org/10.1016/j.physd.2004.01.041

Ouellette N.T. "Probing the Statistical Structure of Turbulence with Measurements of Tracer Particle Tracks." [PhD Dissertation, Cornell University] (2006).
