function H_test(fn)
%H_TEST verifies the color homography color mapping
%
%   H_TEST(FN) shows the results of color homography mapping shown in 
%   Fig. 5 in the PAMI paper. 
%
%   Parameters:
%   FN: image number (see the directory 'H_test_im').

%   Examples:
%   H_test(7);

%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.

%   References:
%   Finlayson, G.D., Gong, H. and Fisher, R., 2017. Color Homography: theory
%   and applications. IEEE Transactions on Pattern Analysis and Machine
%   Intelligence.

addpath('../homo_solver');

%% configuration
dpath = 'H_test_im/';

% load query image
rim = imread([dpath,num2str(fn),'_i110.png']);
rim = im2double(rim);
pr = reshape(rim,[],3)';

qim = imread([dpath,num2str(fn),'_l6c1.png']);
qim = im2double(qim);
pq = reshape(qim,[],3)';

[H,~,D] = uea_H_from_x_als(pr,pq,50);
M = pq/pr;

pc11 = H*pr*D;
pc12 = H*pr;
pc2 = M*pr;
cim11 = reshape(pc11',size(qim));
cim12 = reshape(pc12',size(qim));
cim2 = reshape(pc2',size(qim));

addpath('../utility');
ch1 = chrodist(pr,128);
ch2 = chrodist(pq,128);
figure; ch = imShowPair(ch1,ch2);
rmpath('../utility');

figure;
subplot(2,3,1); imshow(rim); title('A');
subplot(2,3,2); imshow(qim); title('B');
subplot(2,3,3); imshow(cim11); title('A2B w shading');
subplot(2,3,4); imshow(cim12); title('A2B w/o shading');
subplot(2,3,5); imshow(cim2); title('A2B by Least-Squares');

% save figures
%{
imwrite(rim,['/tmp/',num2str(fn),'_A.jpg']);
imwrite(qim,['/tmp/',num2str(fn),'_B.jpg']);
imwrite(cim11,['/tmp/',num2str(fn),'_HD.jpg']);
imwrite(cim12,['/tmp/',num2str(fn),'_H.jpg']);
imwrite(cim2,['/tmp/',num2str(fn),'_LS.jpg']);
imwrite(ch,['/tmp/',num2str(fn),'_ch.jpg']);
%}

rmpath('../homo_solver');

end

