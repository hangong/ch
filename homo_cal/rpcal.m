function M = rpcal(rgb,xyz,Morder)
%RPCAL estimates a root-polynomial color correction matrix
%
%   Parameters:
%   RGB: Input RGB with shadings
%   XYZ: Input XYZ ground truth
%   Morder: order of root-polynomial

%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.

%   References:
%   Gong, H. and Finlayson, G.D., 2017. Root-Polynomial Color Homography
%   Color Correction, AIC Congress.

if nargin<3, Morder = 2; end

terms = build_terms(Morder); % build terms

Nterms = size(terms,1);
Npatch = size(rgb,1);

% extend the poly terms
pRGB = zeros(Npatch,Nterms);
powermul = @(P1,P2) prod(bsxfun(@(A,B) B.^A,P1,P2),2);
for ip = 1:Npatch
    pRGB(ip,:) = powermul(terms,rgb(ip,:));
end

I = eye(Nterms)*1e-7;
M.matrix = (pRGB'*pRGB + I'*I)\(pRGB'*xyz);
M.terms = terms;
M.cfun = @cfun;

function terms = build_terms(Mo)
% BUILD_TERMS computes the angular error
%
% INPUT
% Mo - maximum order of terms
%
% OUTPUT
% terms - coefficients of terms 

terms = zeros(0,3);
for nR = Mo:-1:0
    for nG = Mo-nR:-1:0
        nB = Mo-nR-nG;
        terms(end+1,:) = [nR,nG,nB];
    end
end

terms = terms/Mo;

end

function cXYZ = cfun(cRGB,cM,cterms)
    Nt = size(cterms,1);
    Np = size(cRGB,1);

    prgb = zeros(Np,Nt);
    pmul = @(P1,P2) prod(bsxfun(@(A,B) B.^A,P1,P2),2);
    for i = 1:Np
        prgb(i,:) = pmul(cterms,cRGB(i,:));
    end
    cXYZ = prgb*cM; % convert
end

end
