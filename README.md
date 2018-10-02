# neuron_registration
cellular-resolution slice-to-volume registration of neuronal imaging data

This project provides python and MATLAB implementations of the image registration algorithm described in Stoneking et al. 2018.
It enables images of thin brain slices to be registered against volume images of the same tissue that were previously acquired in vivo.
The method is designed to handle thin (effectively 2D) slices, missing data and nonlinear deformation. 

![single_point_registration_figure](figures/single_point_registration_figure.png?raw=true "Registration at a single point")
