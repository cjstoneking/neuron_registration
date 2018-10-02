# neuron_registration

### register images of thin brain slices to images of the intact brain

This project provides python and MATLAB implementations of the image registration algorithm described in Stoneking et al. (2018).

The algorithm addresses an increasingly common problem in neuroscience: many questions require the experimenter to first perform functional imaging of neural activity *in vivo*, and subsequently to investigate neuroanatomy and/or cellular properties *ex vivo*. To fully exploit the resulting data, it is necessary to identify the same individual neurons in both datasets, i.e. to register images of brain slices acquired *ex vivo* to volume images acquired *in vivo*, with cellular resolution. However, this can be difficult when the slices involved are thin, because they often do not perfectly preserve the structures present in the original tissue.
In response, this algorithm was designed to handle input with a substantial degree of missing data, as well as image modalities that differ between *in vivo* and *ex vivo*, and nonlinear deformation of the slice relative to the volume. 

![single_point_registration_figure](figures/single_point_registration_figure.png?raw=true "Registration at a single point")

Example registration at a single point: local features in a tissue slice (A) are used to identify the best match in the corresponding volume image (B). (C) shows the location of this match within the larger image volume.


Different interpolation methods are provided to enable the full registration transform to be estimated from multiple single-point correspondences. For cases with significant nonlinear deformation of the slices, Gaussian process interpolation is used to estimate a transform that is smooth but has fidelity to the data. In other cases, the slices will approximately occupy different planes within the volume, and a method is provided to find the best-fitting planes.


![interpolation_figure](figures/interpolation_figure.png?raw=true "Interpolation result")

Example reconstruction of a slice image from a volume image, based on interpolation from multiple single-point correspondences. (A) shows the original slice image, red circles mark points that were registered to the 3D volume. (B) shows a reconstruction of the slice from the 3D volume. Differences in fluorescence in the slice seem to be genuine features of the data, most likely due to a time gap between the acquisition of the two datasets.
