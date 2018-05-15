function M = alsRPcal(rgb,xyz,csz)
%ALSRPCAL estimates a root-polynomial color homography matrix for
%   shading-independant color correction.
%
%   Parameters:
%   RGB: Input RGB with shadings
%   XYZ: Input XYZ ground truth
%   CSZ: Color checker size (e.g. [4,6])

%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.

%   References:
%   Gong, H. and Finlayson, G.D., 2017. Root-Polynomial Color Homography
%   Color Correction, AIC Congress.

addpath('../homo_solver');
% use 10 iteration and DCT shading constraint
[M,~,D] = uea_H_from_x_alsr(rgb,xyz,csz,10,'DCT');

% some plotting codes
%{
figure;
I_RGB = reshape(rgb,[csz,3]);
imagesc(I_RGB);axis equal; axis off; 

I_D = reshape(diag(D),csz);
figure;
imagesc(I_D,[0,1.5]); axis equal; axis off; colormap(gray);
%}

rmpath('../homo_solver');

end
