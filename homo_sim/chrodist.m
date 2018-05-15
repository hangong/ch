function chist = chrodist(rgb,Nbin)
% CHIST computes normalised rg chromaticity 2D distribution

% Copyright 2015 Han Gong, University of East Anglia

% exclude saturated pixels
ex_msk = rgb<0.01;
ex_msk = ~logical(sum(ex_msk,1));
ex_msk = ex_msk & mean(rgb,1)>0.05;

C = [1,0,0;0,1,0;1,1,1]; % base converse matrix

pC = C*rgb(:,ex_msk); % chromaticity array
hpC = bsxfun(@rdivide,pC,pC(3,:)); % homogenous ones

% quantitise the chromaticities
hr = linspace(0,1,Nbin+1);
chist = histcn(hpC(1:2,:)',hr,hr)';
% make the probalities sum up to one
chist = chist/max(chist(:));

end