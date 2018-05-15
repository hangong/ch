function [H,err,pD] = uea_H_from_x_als(p1,p2,max_iter,tol,k)

% [H,rms,pa] = uea_H_from_x_als(H0,p1,p2,max_iter,tol,k)
%
% Compute H using alternating least square

% An initial estimate of
% H is required, which would usually be obtained using
% vgg_H_from_x_linear. It is not necessary to precondition the
% supplied points.
%
% The format of the xs is
% [x1 x2 x3 ... xn ; 
%  y1 y2 y3 ... yn ;
%  w1 w2 w3 ... wn]

if nargin<3, max_iter = 50; end
if nargin<4, tol = 1e-20; end
if nargin<5, k = 'lin'; end

[Nch,Npx] = size(p1);

ind1 = sum(p1>0 & p1<Inf,1)==Nch;
ind2 = sum(p2>0 & p2<Inf,1)==Nch;
vind = ind1 & ind2;
kind = find(vind);

if (size(p1) ~= size(p2))
 error ('Input point sets are different sizes!')
end

% definition for Graham
P = p1;
Q = p2;
N = P;
D = speye(Npx);

errs = Inf(max_iter+1,1); % error history

% solve the homography using ALS
n_it = 1; d_err = Inf;
while ( n_it-1<max_iter && d_err>tol )
    n_it = n_it+1; % increase number of iteration

    D = SolveD1(N,Q);

    P_d = P*D;
    if strcmp(k,'lin')
        M = Q(:,vind)/P_d(:,vind); % update M
    else
        K = mean(diag(P_d*P_d'))./1e3;
        M = ((P_d*P_d'+K*eye(Nch))\P_d*(Q'))';
    end
    N = M*P;

    NDiff = (N*D-Q).^2; % difference
    errs(n_it) = mean(mean(NDiff(:,vind))); % mean square error
    d_err = errs(n_it-1) - errs(n_it); % 1 order error
end

H = M;
err = errs(n_it);

pD = D;

%plot(errs); hold on;
%fprintf('ALS %d: %f\n',n_it,errs(n_it));

%figure; imagesc(reshape(diag(D),4,6));

    function D = SolveD1(pp,qq)

        [nCh,nPx] = size(pp);

        p = pp; q = qq;

        % constract RGB index
        Dr = [1:nPx*nCh]'; % row index
        Dc = repmat(1:nPx,[nCh,1]); Dc = Dc(:); % column index

        A = sparse(Dr,Dc,p(:),nPx*nCh,nPx);
        B = q(:);
        A1 = A'*A;

        % compute D
        D = A1\(A'*B);
        D = spdiags(D,0,nPx,nPx);
    end


end
