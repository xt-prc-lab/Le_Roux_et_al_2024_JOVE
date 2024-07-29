# Description:

This repository contains the code used in the paper "Equibiaxial stretching device for high magnification live-cell confocal fluorescence microscopy", by Anabel-Lise Le Roux et al.

It allows to calculate the 2D stretch of an image, with respect to an unstretched reference image. It makes use of the following experimental data: 

  * A 2D microscopy image of fluorescent markers embedded on the surface of an stretched PDMS membrane.
  * Another 2D microscopy image of the relaxed (unstretched) PDMS membrane.

This repository is organized in the following directories:

  * [Stretch_Code](https://github.com/xt-prc-lab/Le_Roux_et_al_2024_JOVE/tree/main/Stretch_Code): contains the code to calculate the 2D stretch field of a membrane.
  * [Examples](https://github.com/xt-prc-lab/Le_Roux_et_al_2024_JOVE/tree/main/Examples): contains example experiments to be analyzed.

# Prerequisites:

The software contained in this repository is written in Matlab. In order to run, it needs a valid installation of:

 * [Matlab](https://www.mathworks.com/products/matlab.html): tested with Matlab versions R2019a, R2020b and R2022a.

It has only been tested under Linux (Gentoo Linux and Ubuntu 20.04) and Windows 10, but it should work in any other operating system with minor modifications.

# Dependencies:

This code makes use of the following software dependencies, that must be downloaded and placed in the corresponding folder:

  * [inpaint_nans.m](https://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans): in-paints over nans in a 2-D array.
  * [natsort.m](https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort) and [natsortfiles.m](https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort): natural-order sort of filenames or filepaths.
  * [removeOutliers_2D.m](https://github.com/FranckLab/FIDIC/blob/master/removeOutliers_2D.m): remove outliers from 2D PIV data.
