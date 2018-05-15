% UEA_RANSACFITHOMOGRAPHY - fits MD homography using RANSAC
%
% Usage:   [H,inlier] = uea_alsransac_luv(x1, x2, t, mit)
%
% Arguments:
%          x1  - 2xN or MxN set of homogeneous points.  If the data is
%                2xN it is assumed the homogeneous scale factor is 1.
%          x2  - 2xN or MxN set of homogeneous points such that x1<->x2.
%          t   - The distance threshold between data point and the model
%                used to decide whether a point is an inlier or not. 
%                Note that point coordinates are normalised to that their
%                mean distance from the origin is sqrt(2).  The value of
%                t should be set relative to this, say in the range 
%                0.001 - 0.01
%          mit - Max number of iteration
%
% Note that it is assumed that the matching of x1 and x2 are putative and it
% is expected that a percentage of matches will be wrong.
%
% Returns:
%          inliers - An array of indices of the elements of x1, x2 that were
%                    the inliers for the best model.
%          H       - The MxM homography such that x2 = H*x1.
%
% See Also: ransac, uea_H_from_x_als

% Han Gong <gong@fedoraproject.org>
% School of Computing Sciences
% University of East Anglia

function [H,inliers] = uea_alsransac_luv(x1, x2, white, t, mit)

    if nargin<5, mit = 500; end

    if ~all(size(x1)==size(x2))
        error('Data sets x1 and x2 must have the same dimension');
    end
    
    [nd,npts] = size(x1);
    if nd < 3
        error('x1 and x2 must have at least 3 rows');
    end

    s = nd+1; % minmum number of points to solve homography
    if npts < s
        error('Must have at least %d points to fit homography',s);
    end

    % generate points combinations
    xcomb = combnk(1:s,3);
    ncomb = size(xcomb,1);

    fittingfn = @wrap_als;
    distfn    = @homogdist;
    degenfn   = @isdegenerate;

    % x1 and x2 are 'stacked' to create a 6xN array for ransac
    [~,inliers] = ransac([x1;x2], fittingfn, distfn,...
                         degenfn, s, t, 0, 100, mit);
    % setxor(1:npts',inliers)
    % bad_patch = [2,5,7,10,11,15,22];
    % inliers = setxor(1:npts',bad_patch);

    % Now do a final least squares fit on the data points considered to
    % be inliers.
    if numel(inliers)>=4
        H = uea_H_from_x_als(x1(:,inliers),x2(:,inliers));
    else
        H = uea_H_from_x_als(x1,x2);
    end

%----------------------------------------------------------------------
% Function to evaluate the forward transfer error of a homography with
% respect to a set of matched points as needed by RANSAC.

function [inliers, H] = homogdist(H, x, t)

    lx1 = x(1:nd,:);   % Extract x1 and x2 from x
    lx2 = x((nd+1):end,:);

    % Calculate, in both directions, the transfered points    
    Hx1 = H*lx1;
    
    % Calculate lab distance
    luv_ref = HGxyz2luv(lx2',white)'; % reference LUV
    luv_est = HGxyz2luv(Hx1',white)'; % reference LUV

    uv_ref = bsxfun(@rdivide,luv_ref(2:3,:),max(luv_ref(1,:),eps));
    uv_est = bsxfun(@rdivide,luv_est(2:3,:),max(luv_est(1,:),eps));

    d = sqrt(sum((uv_ref-uv_est).^2,1));
    inliers = find(d<t);
end
    
%----------------------------------------------------------------------
% Function to determine if a set of 4 pairs of matched points give rise
% to a degeneracy in the calculation of a homography as needed by RANSAC.
% This involves testing whether any 3 of the 4 points in each set is
% colinear. 
     
function r = isdegenerate(x)
    lx1 = x(1:nd,:);   % Extract x1 and x2 from x
    lx2 = x((nd+1):end,:);

    ir1 = arrayfun(@(i) iscolinear_n(lx1(:,xcomb(i,:))), 1:ncomb);
    ir2 = arrayfun(@(i) iscolinear_n(lx2(:,xcomb(i,:))), 1:ncomb);

    r = any([ir1,ir2]);
end
    
function H = wrap_als(x)
    H = uea_H_from_x_als(x(1:nd,:),x((nd+1):end,:));
end

end
