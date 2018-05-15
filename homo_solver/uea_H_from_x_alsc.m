function [H,pD] = uea_H_from_x_alsc(P,Q,csz,max_iter)
%UEA_H_FROM_X_ALSC computes a 3x3 homography matrix with smooth shading 
%constraints.
% 
%   M = UEA_H_FROM_X_ALSC(P,Q,csz) returns a 3x3 colour correction matrix
%   which maps the intensity tripples from P to Q. P and Q are 3xN matrices
%   where N is the number of intentsities. csz is a 1x2 vector which indicates
%   the size of image.
%
%   M = UEA_H_FROM_X_ALSC(P,Q,csz,max_iter,solver,tol) allows users
%   to specify optional parameters as the follows:
%
%       Parameter name   Value
%       --------------   -----
%       'max_iter'       maximum number of iterations (default 30).
%
%   [M,err,pD] = UEA_H_FROM_X_ALSC(P,Q,csz) also returns color correction
%   error and diagonal shading matrix, in the final iteration.
%

%   Copyright 2015 Han Gong <gong@fedoraproject.org>, University of East Anglia  

%   References:
%      TBC

if nargin<4, max_iter = 50; end

[Nch,Npx] = size(P);

ind1 = sum(P>0 & P<Inf,1)==Nch;
ind2 = sum(Q>0 & Q<Inf,1)==Nch;
vind = ind1 & ind2;

if (size(P) ~= size(Q))
 error ('Input point sets are different sizes!')
end

% definition for Graham
N = P;
D = speye(Npx);

% constract M
M1 = ShadingDiff(csz);
% constract RGB index
Dr = [1:Npx*Nch]'; % row index
Dc = repmat(1:Npx,[Nch,1]); Dc = Dc(:); % column index

% solve the homography using ALS
n_it = 1; d_err = Inf;
while (n_it-1<max_iter)
    n_it = n_it+1; % increase number of iteration

    D = SolveD2(N,Q);

    %figure;
    %ImD = full(reshape(diag(D),csz));
    %ImD = bFilter(ImD);
    %D = spdiags(ImD(:),0,Npx,Npx);
    %imagesc(min(ImD,2));

    P_d = P*D;
    M = Q(:,vind)/P_d(:,vind); % update M
    N = M*P;
end

H = M;

pD = D;

%plot(errs); hold on;
%fprintf('ALS %d: %f\n',n_it,errs(n_it));

    function D = SolveD2(lP,lQ)
        p = lP; q = lQ;
        p(:,~vind) = 1; q(:,~vind) = 1;

        A = sparse(Dr,Dc,p(:),Npx*Nch,Npx);
        B = q(:);
        A1 = A'*A;

        % compute D
        lambda = mean(diag(A1))./mean(diag(M1));
        lambda = 1e3*lambda;
        D = (A1+lambda*M1)\(A'*B);
        %D = A1\(A'*B);

        D = spdiags(D,0,Npx,Npx);
    end

    function lM = ShadingDiff(lsz)
    % minimise an edge image

        nel = prod(lsz);
        snel = prod(lsz-1);

        ind = zeros(lsz); ind(:) = 1:nel;
        cdx = ind(2:lsz(1),2:lsz(2)); % centre
        tdx = ind(1:lsz(1)-1,2:lsz(2)); % top
        ldx = ind(2:lsz(1),1:lsz(2)-1); % left

        % flatten index
        cdx = cdx(:); tdx = tdx(:); ldx = ldx(:);

        sMv = sparse(cdx,tdx,-ones(1,snel),nel,nel) + ... % y-dir 
              sparse(cdx,cdx,ones(1,snel),nel,nel);
        sMh = sparse(cdx,ldx,-ones(1,snel),nel,nel) + ... % x-dir 
              sparse(cdx,cdx,ones(1,snel),nel,nel);

        lM = sMv'*sMv + sMh'*sMh;

    end

end
