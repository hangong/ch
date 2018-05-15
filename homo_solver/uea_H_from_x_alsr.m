function [M,err,pD] = uea_H_from_x_alsr(P,Q,csz,max_iter,solver,Morder,ind)
%UEA_H_FROM_X_ALSR computes a Root-Polynomial partial homography matrix.
% 
%   M = UEA_H_FROM_X_ALSR(P,Q,csz) returns a RP colour correction matrix
%   which maps the intensity tripples from P to Q. P and Q are 3xN matrices
%   where N is the number of intentsities. csz is a 1x2 vector which indicates
%   the size of image.
%
%   M = UEA_H_FROM_X_ALSR(P,Q,csz,max_iter,solver,Morder,tol) allows users
%   to specify optional parameters as the follows:
%
%       Parameter name   Value
%       --------------   -----
%       'max_iter'       maximum number of iterations (default 30).
%       'solver'         names of shading solvers include 'DCT', 'DLT', 'edge'.
%                        (default: 'DLT')
%                        'DCT' uses 2D DCT to approximate shading field.
%                        'DLT' uses 2D DLT to approximate shading field.
%       'Morder'         order of roor-polynomial intentisty terms (default 2).
%
%   [M,err,pD] = UEA_H_FROM_X_ALSR(P,Q,csz) also returns color correction
%   error and diagonal shading matrix, in the final iteration.
%

%   Copyright 2015-2017 Han Gong <gong@fedoraproject.org>,
%   University of East Anglia  

%   References:
%      TBC

if (size(P) ~= size(Q))
 error ('Input point sets are different sizes!')
end

% size of valid intensities
[r,c] = size(P);

if nargin<4, max_iter = 30; end
if nargin<5, solver = 'DCT'; end
if nargin<6, Morder = 2; end
if nargin<7, ind = 1:prod(csz); end

P = P'; Q = Q';

M.terms = build_terms(Morder); % build RP terms exponential factors
M.cfun = @cfun; % RP homography mapping function

Nterms = size(M.terms,2);
Npatch = size(P,2);

% initialisation
D = speye(c);
N = P;
errs = Inf(max_iter+1,1); % error history

gradient = create_basis(csz,solver); % DCT gradient basis
GRADbasis = zeros(size(gradient));
num = size(gradient,3);

% extend the poly terms
pP = zeros(Nterms,Npatch);
for it = 1:Nterms
    ina = find(M.terms(:,it)==1);
    if isempty(ina)
        pP(it,:) = prod(P.^repmat(M.terms(:,it),[1,Npatch]));
    else
        pP(it,:) = P(ina,:);
    end
end

% solve the homography using ALS
n_it = 1;
while ( n_it-1<max_iter)
    n_it = n_it+1; % increase number of iteration

    switch solver
    case {'DCT','DLT'}
        imnew = reshape(N',[csz,3]);
        coef = gradbuilder(imnew,gradient,Q');

        for j = 1:num % build shading gradient
            GRADbasis(:,:,j) = gradient(:,:,j).*coef(j);
        end
        scs = sum(GRADbasis,3); % built shading
        D = diag(scs(:));
    end

    tpP = pP*D;
    H = Q(:,ind)/(tpP(:,ind)); % update H
    N = H*pP; % apply perspective transform

    errs(n_it) = meansqr(N*D-Q); % mean square error
end

if nargout>2, pD = D; end

%figure; plot(errs,'g');
%fprintf('ALS_RP %d: %f\n',n_it,errs(n_it));
%figure; imagesc(reshape(diag(D),4,6));

M.matrix = H;
err = errs(n_it);

function terms = build_terms(Mo)
% BUILD_TERMS computes the angular error
%
% INPUT
% Mo - maximum order of terms
%
% OUTPUT
% terms - exponential coefficients of terms 

    terms = zeros(0,3);
    for nR = Mo:-1:0
        for nG = Mo-nR:-1:0
            nB = Mo-nR-nG;
            terms(end+1,:) = [nR,nG,nB];
        end
    end

    terms = terms/Mo;
    terms = terms'; % transpose

    [~,tind] = sort(max(terms),2,'descend');
    terms = terms(:,tind);

end

function cXYZ = cfun(cRGB,cM,cterms)

    cRGB = cRGB';
    
    Nt = size(cterms,2);
    Np = size(cRGB,2);

    prgb = zeros(Nt,Np);
    for tit = 1:Nt
        tina = find(cterms(:,tit)==1);
        if isempty(tina)
            prgb(tit,:) = prod(cRGB.^repmat(cterms(:,tit),[1,Np]));
        else
            prgb(tit,:) = cRGB(tina,:);
        end
    end

    cXYZ = (cM*prgb)'; % convert

end

function basis = create_basis(sz,solver)
% function to cache DCT basis
%
% INPUT
% sz - size of the chart
% solver - name of solver
%
% OUTPUT
% basis_dct - DCT basis planes

switch solver
case 'DLT'
    f = @dltmtx;
case 'DCT'
    f = @dctmtx;
end

P1 = f(sz(1));
P2 = f(sz(2));

order = 2;
basis = zeros([sz,sum(1:order)]);
sd = find(sz==1);

basis(:,:,1) = P1(1,:)'*P2(1,:);
if isempty(sd) | sd~=1, basis(:,:,2) = P1(2,:)'*P2(1,:); end
if isempty(sd) | sd~=2, basis(:,:,3) = P1(1,:)'*P2(2,:); end
%basis(:,:,4) = P1(1,:)'*P2(3,:);
%basis(:,:,5) = P1(2,:)'*P2(2,:);
%basis(:,:,6) = P1(3,:)'*P2(1,:);
%basis(:,:,7) = P1(4,:)'*P2(1,:);
%basis(:,:,8) = P1(3,:)'*P2(2,:);
%basis(:,:,9) = P1(2,:)'*P2(3,:);
%basis(:,:,10) = P1(1,:)'*P2(4,:);
%basis(:,:,11) = P1(1,:)'*P2(5,:);

end

function im_shaded = fullimagegradienter(im,gradient)
% FULLIMAGEGRADIENTER applies a shading gradient to the image
%
% INPUT
% im - unshaded image
% gradient - DCT shading basis
% 
% im_shaded - shaded image

    im_shaded = zeros(size(im,1),size(im,2),3);
    for lc = 1:3
        im_shaded(:,:,lc) = im(:,:,lc).*gradient;      
    end
end

function coef = gradbuilder(im,gradient,XYZ)
% GRADBUILDER estimates the DCT approximated gradient
%
% INPUT
% im - unshaded image
% gradient - DCT shading basis
% XYZ - reference XYZs
%
% OUTPUT
% coef - coefficents of DCT basis

    Nb = size(gradient,3); % number of DCT basis
    A = zeros(numel(im),Nb);
    for lj = 1:Nb
        im_shaded = fullimagegradienter(im,gradient(:,:,lj));
        A(:,lj) = im_shaded(:);
    end

    A1 = A'*A;
    coef = (A1 + eye(size(A,2))*mean(diag(A1))*1e-3) \ (A'*XYZ(:));
end

end
