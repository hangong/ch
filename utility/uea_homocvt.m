function q = uea_homocvt(p,H,C)
% HOMOCVT applies homography transform

if nargin<3, C = eye(size(H)); end

pr = reshape(p,[],size(H,1))';
prC = C*pr; % chromaticity array
qC = H*prC; % apply inverse homography
q = (C\qC)'; % convert back to rgb
