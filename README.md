# Color Homography Theorem
Demostration of color homography Theorem with experiments in color correction, color feature matching, and color indexing, published in TPAMI 2017.

Images of co-planar points in 3-dimensional space taken from different camera positions are a homography apart. Homographies are at the heart of geometric methods in computer vision and are used in geometric camera calibration, 3D reconstruction, stereo vision and image mosaicking among other tasks. In this paper we show the surprising result that homographies are the apposite tool for relating image colors of the same scene when the capture conditions ‒ illumination color, shading and device ‒ change.

[[Paper](http://homepages.inf.ed.ac.uk/rbf/PAPERS/cshgong17.pdf)]
[[Supplementary Material](http://www2.cmp.uea.ac.uk/~ybb15eau/projects/color_homography/supl_review.pdf)]

## References
Please cite the following papers if you find the code useful.
1. Finlayson, G.D., Gong, H. (**Co-first author**) and Fisher, R., 2017. Color Homography: theory and applications. IEEE Transactions on Pattern Analysis and Machine Intelligence.
2. Finlayson, G.D., Gong, H. and Fisher, R., 2016. Color Homography Color Correction, Color and Imaging Conference.
3. Gong, H. and Finlayson, G.D., 2017. Root-Polynomial Color Homography Color Correction, AIC Congress.

# Setup
The code was tested on Ubuntu 16.04 with MATLAB 2016b.

## Download datasets
1. Download files
* [24-patch color checker with non-uniform shading](http://www2.cmp.uea.ac.uk/~ybb15eau/db/HG_ColourChecker.zip)
* ALOI object recognition dataset (subset) [col](http://aloi.science.uva.nl/tars/aloi_col.tar) [ill](http://aloi.science.uva.nl/tars/aloi_ill.tar) [view](http://aloi.science.uva.nl/tars/aloi_view.tar)

2. Make a directory named 'data' and extract the above packages to the 'data' dir. We will have the following structures:
```
./data/HG_HG_ColourChecker
./data/aloi/col/
./data/aloi/ill/
./data/aloi/view/
```
We need to rename some extracted directories. For instance, the extracted 'aloi_col' directory becomes a sub-directory './data/aloi/col/'.
## Compile ASIFT
```
$ cd ./homo_sim/ASIFT
$ make
```

# Color structure matching
Features of chromaticity histogram can be used for estimating a color homography which does a color mapping.
## Color structure matching by ASIFT
see `./homo_sim/H_test.m` which reproduces Fig.3 of the paper.

<img src="http://www2.cmp.uea.ac.uk/~ybb15eau/single-color.jpg" width="500">

## Color indexing test
We implemented 3 other color indexing methods:
* Swain's histogram intersection (`./homo_sim/swain_index.m`)
* Gevers's m1m2m3 color feature (`./homo_sim/gevers_index.m`)
* Finlayson's comprehensive normalization (`./homo_sim/cn_index.m`)

The color homography color indexing is implemented in `./homo_sim/homo_index.m`. They share the similar interface. For example, to get the result for swain's indexing method:
```
>> MP = swain_index('aloi/col');
```
where `aloi/col` is the color temperature test directory, `MP`'s rows refer to object ID and its columns refer to capture condition. See the [ALOI dataset website](http://aloi.science.uva.nl/) for the details.

# Color Correction
Given pixel-to-pixel correspondances, color homography is applied for shading-invariant color mapping.

## Shading invariant color mapping
see `./homo_sim/dismatch_im.m` which reproduces Fig.5 of the paper.

## Color correction test
We implemented 2 other color correction methods:
* Linear least-squares (`homo_cal/lscal.m`)
* Alternating least-squares (`homo_cal/alshomocal.m`)
* Root-polynomial least-squares (`homo_cal/rpcal.m`)

We also show 2 color homography approaches:
* 2D color homography with RANSAC (`homo_cal/ransachomocal_luv.m`) [1,2]
* 2nd order root-polynomial homography (`homo_cal/alsRPcal.m`) [3]

*Please note that the results here are slightly different from the paper's due to a preivous evaluation bug. However, the relative difference is still similar.*

To evaluate for all methods, in MATLAB execute
```
>> cal_real
```

Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East Anglia.
