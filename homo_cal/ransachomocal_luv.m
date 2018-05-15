function M = ransachomocal_luv(rgb,xyz,white,rgb_u)
%RANSACHOMOCAL_LUV estimates a 'good-luck' color homography matrix for
%   shading-independant color correction.
%
%   RANSACHOMOCAL_LUV(RGB,XYZ,WHITE,RGB_U) returns the 'best-luck' color 
%   homography color correction matrix. 
%
%   Parameters:
%   RGB: Input RGB with shadings
%   XYZ: Input XYZ ground truth
%   WHITE: white point RGB for evaluation
%   RGB_U: shading removed RGB for evaluation

%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.

%   References:
%   Finlayson, G.D., Gong, H. and Fisher, R., 2017. Color Homography: theory
%   and applications. IEEE Transactions on Pattern Analysis and Machine
%   Intelligence.

%   Note:
%   Due to the randomness of RANSAC, we show the best theoretically 
%   'good-luck' estimation result. In practice, RANSAC does not guarantee
%   the optimum result. A more robust scheme to control the randomness is
%   the future work.

addpath('../homo_solver');

% best error
b_err = Inf;
n_trail = 100;

% try 100 times
for i = 1:n_trail
    t_M = uea_alsransac_luv(rgb',xyz',white',0.2);
    xyz_est = uea_homocvt(rgb_u,t_M);
    t_err = mean(luv_err(xyz_est,xyz));

    if t_err<b_err
        % return the good CC matrix
        M = t_M; b_err = t_err;
    end
end

rmpath('../homo_solver');

end

function err = luv_err(xyz_est,xyz_std)
    % normalize by a white patch's green intensity
    XYZ_est = xyz_est./xyz_est(4,2);
    
    % LUV error
    luv_est = HGxyz2luv(XYZ_est,xyz_std(4,:));
    luv_ref = HGxyz2luv(xyz_std,xyz_std(4,:)); % reference LUV
    err = deltaE1976(luv_ref,luv_est);
end

