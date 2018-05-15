function MP = gevers_index(dbname,gamma,Nbin)
%GEVERS_INDEX indexes images by color using m1m2m3 features. 
%
%   GEVERS_INDEX(DBNAME,GAMMA,NBIN) returns a matching percentile matirx for a
%   dataset.
%
%   Parameters:
%   DBNAME: Name of database
%   GAMMA: display gamma (default: 1)
%   NBIN: number of histogram bins (defualt: 32)

%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.

%   References:
%   Gevers, T. and Smeulders, A.W., 1998, January. Image indexing using 
%   composite color and shape invariant features. In Computer Vision, 1998. 
%   Sixth International Conference on (pp. 576-581). IEEE.


addpath('../utility');

if ~exist('gamma','var'), gamma = 1; end
if ~exist('Nbin','var'), Nbin = 32; end

% configuration
mlpath = 'model_cache_gevers/'; % model cache path
dcpath = 'db_cache_gevers/'; % dataset cache

[pm,rm,objsel] = dbparser(dbname); % load data path
objsel = 1:size(pm,1); % enforce all object test

No = size(objsel,2); % number of objects
Nc = size(pm,2); % number of conditions
MP = NaN(No,Nc); % MP for each query
[t1,t2] = gevers_hist(dbname,Nbin);
binrng = {t1,t2};

% build model caches
mpath = [mlpath,dbname];
if ~exist(mpath,'dir'), mkdir(mpath); end

for oi = objsel
    mlName = [mlpath,dbname,'/',num2str(oi),'.mat'];
    if ~exist(mlName,'file')
        rim = imread(rm{oi}.d); % load ref image
        rim = im2double(rim).^gamma; % de-gamma
        % m1m2
        pr = cm(rim);
        % build histogram
        mh = rgb2hist(pr,binrng);
        save(mlName,'mh');
    end
end

% builde cache dir
qpath = [dcpath,dbname,'/'];
if ~exist(qpath,'dir'), mkdir(qpath); end

% query tests

% each query image is tested agains all models
for oo = 1:No % for each object
    oi = objsel(oo);
    for ci = 1:Nc % for each condition
        % test if the test condition exists
        if isempty(pm{oi,ci}), continue; end

        oim = imread(pm{oi,ci}.d); % load query image
        oim = im2double(oim).^gamma;
        po = cm(oim);

        % if query cache does not exist
        qName = [qpath,num2str(oi),'_',num2str(ci),'.mat'];
        qh = rgb2hist(po,binrng);
        % m1m2
        save(qName,'qh');

        score_t = zeros(No,1); % cache for matching score

        if ~exist(qName,'file'), error('file not found'); end
        load(qName,'qh');

        % go through all models to compare
        parfor pp = 1:No % for each model to compare
            Pi = objsel(pp);
            % load referencing chromaticity feature profile
            mlName = [mlpath,dbname,'/',num2str(Pi),'.mat'];
            if ~exist(mlName,'file'), error('file not found'); end
            l = load(mlName,'mh');       

            % histogram intersection
            score_t(pp) = sum(min(l.mh,qh))/sum(l.mh);
        end

        % get match percentile for this query
        [~,ind] = sort(score_t,'descend');
        MP(oo,ci) = (No-find(ind==oo))/(No-1);
        fprintf('(%d,%d) = %.4f\n',oi,ci,MP(oo,ci));
    end
end

rmpath('../utility');

end

function ohist = rgb2hist(rgb,binrng)
    addpath('../utility');
    ohist = histcn(rgb(:,1:2),binrng{:}); % use RG tuple
    rmpath('../utility');
end

function crgb = cm(im)

    imsz = size(im);    

    % replace low intensity pixels
    oim = im;
    oim(oim<1/255) = 1/255;

    lim = log(oim); % convert to log domain
    % log difference image
    rim = zeros([imsz(1:2),2]);
    dim = zeros([imsz(1:2),2]);
    rim(:,:,1) = lim(:,:,1) - lim(:,:,2);
    rim(:,:,2) = lim(:,:,1) - lim(:,:,3);
    % edge detection
    for ch = 1:2
        dim(:,:,ch) = canny(rim(:,:,ch),1.0);
    end

    crgb = reshape(dim,[],2);

end

function BW = canny(I,sigma)
% canny gradient generator

    si = ceil(sigma*3)*2+1;
    % Convolution of image by Gaussian Coefficient
    B = fspecial('gaussian',si,sigma);
    A = conv2(I, B, 'same');

    % Convolution by image by horizontal and vertical filter
    x = -si:1:si;
    y = -si:1:si;
    gaussx = -(x/(sigma*sigma)).*exp(-(x.*x+y.*y)/(2*sigma*sigma));
    gaussy = gaussx';
    % Convolution by image by horizontal and vertical filter
    Filtered_X = conv2(A, gaussx, 'same');
    Filtered_Y = conv2(A, gaussy, 'same');

    % Calculate directions/orientations
    arah = atan2(Filtered_Y, Filtered_X);
    arah = arah*180/pi;

    pan = size(A,1);
    leb = size(A,2);

    % Adjustment for negative directions, making all directions positive
    for i=1:pan
        for j=1:leb
            if (arah(i,j)<0) 
                arah(i,j) = 360+arah(i,j);
            end;
        end;
    end;

    arah2 = zeros(pan, leb);

    % Adjusting directions to nearest 0, 45, 90, or 135 degree
    for i = 1  : pan
        for j = 1 : leb
            if ((arah(i, j) >= 0 ) && (arah(i, j) < 22.5) ||...
                (arah(i, j) >= 157.5) && (arah(i, j) < 202.5) ||...
                (arah(i, j) >= 337.5) && (arah(i, j) <= 360))
                arah2(i, j) = 0;
            elseif ((arah(i, j) >= 22.5) && (arah(i, j) < 67.5) ||...
                (arah(i, j) >= 202.5) && (arah(i, j) < 247.5))
                arah2(i, j) = 45;
            elseif ((arah(i, j) >= 67.5 && arah(i, j) < 112.5) ||...
                (arah(i, j) >= 247.5 && arah(i, j) < 292.5))
                arah2(i, j) = 90;
            elseif ((arah(i, j) >= 112.5 && arah(i, j) < 157.5) ||...
                (arah(i, j) >= 292.5 && arah(i, j) < 337.5))
                arah2(i, j) = 135;
            end;
        end;
    end;

    % Calculate magnitude
    magnitude = (Filtered_X.^2) + (Filtered_Y.^2);
    magnitude2 = sqrt(magnitude);

    BW = zeros (pan, leb);

    % Non-Maximum Supression
    for i = 2:pan-1
        for j = 2:leb-1
            if (arah2(i,j) == 0)
                BW(i,j) = (magnitude2(i,j) == ...
                max([magnitude2(i,j), magnitude2(i,j+1), magnitude2(i,j-1)]));
            elseif (arah2(i,j) == 45)
                BW(i,j) = (magnitude2(i,j) == ...
                max([magnitude2(i,j), magnitude2(i+1,j-1), magnitude2(i-1,j+1)]));
            elseif (arah2(i,j) == 90)
                BW(i,j) = (magnitude2(i,j) == ...
                max([magnitude2(i,j), magnitude2(i+1,j), magnitude2(i-1,j)]));
            elseif (arah2(i,j) == 135)
                BW(i,j) = (magnitude2(i,j) == ...
                max([magnitude2(i,j), magnitude2(i+1,j+1), magnitude2(i-1,j-1)]));
            end;
        end;
    end;

    BW = BW.*magnitude2;

end

