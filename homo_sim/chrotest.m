function chrotest(chart,illu,Nt)
addpath('../homo_solver');

global deb;

switch chart
case 'SG'
    load illu20SG;
    load SG;
case 'DC'
    load illu20DC;
    load DC;
otherwise
    disp('chart not available');
end
chart_size = size(RGB);
chart_size = chart_size(1:2); % chart size

% apply illumination
if size(illu,2)<2
    %illu_sc = rand(3,1);
    illu_sc = [0.3,0.5,0.7];
    %illu_sc = ones(3,1);
    imi = img.*repmat(gradients(:,:,illu),[1,1,3]);
    for ch = 1:3, imi(:,:,ch) = imi(:,:,ch)*illu_sc(ch); end
else
    imi = img.*gradients(:,:,illu);
end

% apply noise
%imi = imnoise(imi,'poisson');

% downsampling
VER = version; % software version
if any(regexp(VER,'[a-z,A-Z]')) % it's MATLAB
    im = imresize(imi,chart_size,'box');
else % should be octave
    im = imresize(imi,chart_size,'linear');
end

% get their xyz
if deb
    figure('Name','two illuminations');
    subplot(1,2,1); imshow(img);
    subplot(1,2,2); imshow(imi);
end

RGB = RGB/max(RGB(:));
% convert to homogeous coodrinates
rgbhomog = bsxfun(@rdivide,RGB,RGB(:,:,3));
rgbhomo = bsxfun(@rdivide,im,im(:,:,3));
xyzhomog = bsxfun(@rdivide,XYZ,XYZ(:,:,3));

% flatten the data
rgbhomog = reshape(rgbhomog,[],3);
rgbhomo = reshape(rgbhomo,[],3);
xyzhomog = reshape(xyzhomog,[],3);

% original data
rgbg = reshape(RGB,[],3); % RGB GT
xyzg = reshape(XYZ,[],3); % XYZ GT

Ne = 2; % number of error types
Nm = 3; % number of test variables
Nc = 2^Nm; % number of test combinations
H_u = cell(Nc,Nt); % homography under uniform illumination
H_n = cell(Nc,Nt); % homography under non-uniform illumination
method = {@homo_als, @(p1,p2,p3,p4) ransacfithomo_als(p1,p2,0.005,p3,p4)};
% compute uniform homography
for t = 1:Nt
    for nor = 0:1 % normalisation
    for init = 0:1 % initilisation
    for m = 0:1 % methods
        cid = polyval([nor,init,m],2)+1;
        % compute uniform illumination
        H_u{cid,t} = method{m+1}(rgbhomog',xyzhomog',nor,init);
        % compute non-uniform illumination
        H_n{cid,t} = method{m+1}(rgbhomog',xyzhomog',nor,init);
    end
    end
    end
end

% evaluation
E = zeros(Nt,2,Nc,Ne); % homography under uniform illumination
for i = 1:Nc
    for t = 1:Nt
        % compute mapped xyz
        mxyz_u = rgbhomog*(H_u{i}');
        mxyz_n = rgbhomo*(H_n{i}');
        % compute mapping errors
        E(t,1,i,:) = uea_homo_error(mxyz_u,xyzhomog);
        E(t,2,i,:) = uea_homo_error(mxyz_n,xyzhomog);
    end
end

% display the results
ttitle = {'Uniforum Illum.','Non-Uniforum Illum.'};
tlabel = {'\Delta xy','\Delta uv'};
figure('Name','Mapping Error');
for ic = 1:2 % illumination condition
    subplot(1,2,ic);
    te = squeeze(mean(E(:,ic,:,:),1));
    bar(0:(Nc-1),te,'stacked');
    legend(tlabel);
    title(ttitle{ic});
    xlim([-1,Nc]);
end

rmpath('../homo_solver');
