function M = alshomocal(rgb,xyz)
%ALSHOMOCAL estimates a color homography matrix for shading-independant
%   color correction. This version is without outlier detection.
%
%   ALSHOMOCAL(RGB,XYZ) returns the color homography color correction
%   matrix. 
%
%   Parameters:
%   RGB: Input RGB with shadings
%   XYZ: Input XYZ ground truth

%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.

%   References:
%   Finlayson, G.D., Gong, H. and Fisher, R., 2016. Color Homography, 
%   Progress in Colour Studies (PICS).

addpath('../homo_solver');

M = uea_H_from_x_als(rgb',xyz',50,1e-20);

rmpath('../homo_solver');
