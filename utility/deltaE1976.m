function dE = deltaE1976(lab_ref,lab_est)
% DELTAE1976 computes the 1976 colour difference
%
% INPUT
% lab_ref - reference LABs
% lab_est - estimated LABs
% 
% OUTPUT
% dE - colour difference
%
% Copyright 2015 Han Gong, University of East Anglia

dE = sqrt(sum((lab_ref - lab_est).^2,2));
