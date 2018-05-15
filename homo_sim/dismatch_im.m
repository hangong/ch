%This scripts generates similar color homography plots in Fig.3 of the PAMI
% paper.
%
%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.

%   References:
%   Finlayson, G.D., Gong, H. and Fisher, R., 2017. Color Homography: theory
%   and applications. IEEE Transactions on Pattern Analysis and Machine
%   Intelligence.
clear all;

addpath('../homo_solver');
addpath('../utility');

tpath = 'H_test_im/';

gamma = 1;

% original linear image
oim = imread([tpath,'698_r45.png']);
rim = imread([tpath,'698_i110.png']);

%rim = imread([tpath,'444_r45.png']);
%oim = imread([tpath,'444_i110.png']);

%rim = imread([tpath,'881_i250.png']);
%oim = imread([tpath,'881_i110.png']);

%rim = imread([tpath,'3_i250.png']);
%oim = imread([tpath,'3_i110.png']);

oim = im2double(oim).^gamma;
rim = im2double(rim).^gamma;

res = 320;
% compute chromaticity distributions for both images
H = H_from_feat(oim,rim,res,true,4);

rmpath('../homo_solver');
rmpath('../utility');
