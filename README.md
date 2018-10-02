# neuron_registration

### register images of tissue slices to images 

This project provides python and MATLAB implementations of the image registration algorithm described in Stoneking et al. (2018).

The algorithm addresses an increasingly common problem in neuroscience: many questions require the experimenter to first perform functional imaging of neural activity *in vivo*, and subsequently to investigate neuroanatomy and/or cellular properties *ex vivo*. To fully exploit the resulting data, it is necessary to identify the same individual neurons in both datasets, i.e. to register images of brain slices acquired *ex vivo* to volume images acquired *in vivo*, with cellular resolution. However, this can be difficult when the slices involved are thin, because they often do not perfectly preserve the structures present in the original tissue.
In response, this algorithm was designed to handle input with a substantial degree of missing data, as well as image modalities that differ between *in vivo* and *ex vivo*, and nonlinear deformation of the slice relative to the volume. 

![single_point_registration_figure](figures/single_point_registration_figure.png?raw=true "Registration at a single point")

Example registration at a single point: local features in a tissue slice (A) are used to identify the best match in the corresponding volume image (B). (C) shows the location of this match within the larger image volume.




