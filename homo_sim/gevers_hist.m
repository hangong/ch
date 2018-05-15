function [t1,t2] = gevers_hist(dbname,Nbin_o,gamma)
%GEVERS_HIST estimates the optimum rg histogram intervals for m1m2m3 
% feature counts.
%
%   GEVERS_HIST(DBNAME,GAMMA,NBIN) two histogram interval ticks for m1m2m3
%   indexing.
%
%   Parameters:
%   DBNAME: Name of database
%   GAMMA: display gamma (default: 1)
%   NBIN_O: number of histogram bins (defualt: 32)

%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.

%   References:
%   Gevers, T. and Smeulders, A.W., 1998, January. Image indexing using 
%   composite color and shape invariant features. In Computer Vision, 1998. 
%   Sixth International Conference on (pp. 576-581). IEEE.

addpath('../utility');

if ~exist('gamma','var'), gamma = 1; end
if ~exist('Nbin_o','var'), Nbin_o = 32; end

[pm,rm,objsel] = dbparser(dbname); % load data path
objsel = 1:size(pm,1); % enforce all object test

Nbin = 1000;
H1 = zeros(1,Nbin);
H2 = zeros(1,Nbin);

m_max = 0.6;
hr = linspace(0,m_max,Nbin+1);
ht = hr(1:end-1)+m_max/Nbin/2;

if ~exist('m1m2m3_hist_cache.mat','file')
    for oi = objsel
        rim = imread(rm{oi}.d); % load ref image
        rim = im2double(rim).^gamma; % de-gamma
        % m1m2
        pr = cm(rim);
        % accumulate values
        t = pr(:,1);
        t = t(isfinite(t)&t>0);
        H1 = H1 + histcounts(t,hr); 
        t = pr(:,2);
        t = t(isfinite(t)&t>0);
        H2 = H2 + histcounts(t,hr); 
    end

    save('m1m2m3_hist_cache.mat','H1','H2');
else
    load('m1m2m3_hist_cache.mat','H1','H2');
end

AHr1 = cumsum(H1);
AHr2 = cumsum(H2);

hH1 = linspace(AHr1(end)/Nbin_o,AHr1(end),Nbin_o);
hH2 = linspace(AHr2(end)/Nbin_o,AHr2(end),Nbin_o);

m1 = H1>0;
m2 = H2>0;

t1 = [eps,interp1(AHr1(m1),ht(m1),hH1,'linear','extrap')]; t1(end) = m_max;
t2 = [eps,interp1(AHr2(m2),ht(m2),hH2,'linear','extrap')]; t2(end) = m_max;

rmpath('../utility');

function crgb = cm(im)

    imsz = size(im);

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


end
