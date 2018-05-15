% ISCOLINEAR_N - are 3 points colinear (ND version)
%
% Usage:  r = iscolinear_n(P)
%
% Arguments:
%        P - Points in ND (Nx3).

%
% Returns:
%        r = 1 if points are co-linear, 0 otherwise

% Copyright (c) 2015 Han Gong 
% University of East Anglia
% 

function r = iscolinear_n(P)

    if ~(size(P,1)>=3)                              
        error('points must have the same dimension of at least 3');
    end
    
	r =  norm(cross(P(:,2)-P(:,1), P(:,3)-P(:,1))) < eps;
