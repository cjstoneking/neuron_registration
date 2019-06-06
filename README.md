# neuron_registration

### register images of thin brain slices to images of the intact brain

This project provides python and MATLAB implementations of the image registration algorithm described in Stoneking et al. (2019).

The algorithm addresses an increasingly common problem in neuroscience: many questions require the experimenter to first perform functional imaging of neural activity *in vivo*, and subsequently to investigate neuroanatomy and/or cellular properties *ex vivo*. To fully exploit the resulting data, it is necessary to identify the same individual neurons in both datasets. This can be achieved by registering images of brain slices acquired *ex vivo* to volume images acquired *in vivo*, i.e. identifying matching features in these sets of images and thus estimating a transform that maps the slices to the volume. Although registration is a common image analysis problem, it can be challenging in this case because the slices involved may be quite thin, and often do not perfectly preserve the structures present in the original tissue.
In response, this algorithm was designed perform slice-to-volume image registration on input with a substantial degree of missing data. It is robust enough to handle image modalities that differ between *in vivo* and *ex vivo*, and flexible enough to infer nonlinear deformation of the slice relative to the volume. 

![single_point_registration_figure](figures/single_point_registration_figure.png?raw=true "Registration at a single point")

Example registration at a single point: local features in an image of a tissue slice (A) are used to identify the best match in the corresponding volume image (B). (C) shows the location of this match within the larger image volume. The signal in these images is from labeled blood vessels.

The registration process involves first performing a number of single-point registrations, in which a single point in the slice is chosen and its best match in the volume is identified. This is achieved by maximizing local cross-correlation, using a Bayesian optimization approach to minimize computation time. An example single-point registration is shown above.
In the second stage, a global registration between the slice and the volume is obtained by interpolating between the single-point registrations. We provide two different interpolation options: *best-fitting plane*, which is appropriate if the slice has little nonlinear deformation relative to the volume, and *Gaussian process*, which can capture nonlinear deformation by estimating a smooth nonlinear transform.


![interpolation_figure](figures/interpolation_figure.png?raw=true "Interpolation result")

Example reconstruction of a slice image from a volume image, based on interpolation from multiple single-point correspondences. (A) shows the original slice image, red circles mark points that were registered to the 3D volume. (B) shows a reconstruction of the slice from the 3D volume, based on the best-fitting plane. The signal in these images is from labeled neurons. Differences in fluorescence in the slice are due to differences in the signal from individual cells and seem to be genuine features of the data, most likely due to a time gap between the acquisition of the two datasets.
