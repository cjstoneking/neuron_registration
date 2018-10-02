# neuron_registration

### register images of tissue slices to images 

This project provides python and MATLAB implementations of the image registration algorithm described in Stoneking et al. 2018.
The algorithm enables images of thin histological slices of tissue to be registered against volume images of the same tissue that were previously acquired *in vivo* - i.e. it infers a transform that warps the slices to the volume. data from histology to be integrated with data obtained *in vivo*. One possible application would be to integrate in vivo functional imaging of neural activity with ex vivo antibody staining or in situ sequencing. 
The method is designed to handle thin (effectively 2D) slices, missing data, varying image modalities and nonlinear deformation. In particular, the ability to cope with missing data allows this method to be used on thin tissue slices which usually do not perfectly preserve the structures present in the original tissue.

![single_point_registration_figure](figures/single_point_registration_figure.png?raw=true "Registration at a single point")
