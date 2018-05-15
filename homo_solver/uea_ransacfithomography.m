% UEA_RANSACFITHOMOGRAPHY - fits MD homography using RANSAC
%
% Usage:   [inliers,H] = uea_ransacfithomography(x1, x2, t, mit)
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
% See Also: ransac, homography2d, homography1d

% Han Gong <gong@fedoraproject.org>
% School of Computing Sciences
% University of East Anglia

function [inliers,H] = uea_ransacfithomography(x1, x2, t, mit)

    if nargin<4, mit = 5000; end

    if ~all(size(x1)==size(x2))
        error('Data sets x1 and x2 must have the same dimension');
    end
    
    [nd,npts] = size(x1);
    if nd < 3 
        error('x1 and x2 must have at least 3 rows');
    end

    s = nd+1; % minmum number of points to solve homography
    if npts < s
        em = sprintf('Must have at least %d points to fit homography',s);
        error(em);
    end
    
    % generate points combinations
    xcomb = combnk(1:s,3);
    ncomb = size(xcomb,1);

    fittingfn = @wrap_vgg_homographynd;
    distfn    = @homogdistnd;
    degenfn   = @isdegenerate;
    
    % x1 and x2 are 'stacked' to create a 6xN array for ransac
    [~,inliers] = ransac([x1;x2], fittingfn, distfn,...
                         degenfn, s, t, 0, 100, mit);

    % Now do a final least squares fit on the data points considered to
    % be inliers.
    if nargout > 1
        if numel(inliers)>=4
            H = uea_H_from_x_als(x1(:,inliers),x2(:,inliers));
        else
            H = uea_H_from_x_als(x1,x2);
        end
    else
        H = [];
    end

%----------------------------------------------------------------------
% Function to evaluate the symmetric transfer error of a homography with
% respect to a set of matched points as needed by RANSAC.

function [inliers, H] = homogdistnd(H, x, t)

    x1 = x(1:nd,:);   % Extract x1 and x2 from x
    x2 = x((nd+1):end,:);    
    
    % Calculate, in both directions, the transfered points    
    Hx1    = H*x1;
    invHx2 = H\x2;
    
    % Normalise so that the homogeneous scale parameter for all coordinates
    % is 1.
    
    x1     = hnormalise(x1);
    x2     = hnormalise(x2);
    Hx1    = hnormalise(Hx1);
    invHx2 = hnormalise(invHx2); 
    
    d2 = sum((x1-invHx2).^2)  + sum((x2-Hx1).^2);
    inliers = find(abs(d2) < t);
end
    
%----------------------------------------------------------------------
% Function to determine if a set of 4 pairs of matched points give rise
% to a degeneracy in the calculation of a homography as needed by RANSAC.
% This involves testing whether any 3 of the 4 points in each set is
% colinear. 
     
function r = isdegenerate(x)
    x1 = x(1:nd,:);   % Extract x1 and x2 from x
    x2 = x((nd+1):end,:);

    ir1 = arrayfun(@(i) iscolinear_n(x1(:,xcomb(i,:))), 1:ncomb);
    ir2 = arrayfun(@(i) iscolinear_n(x2(:,xcomb(i,:))), 1:ncomb);

    r = any([ir1,ir2]);
end
    
function H = wrap_vgg_homographynd(x)
    H = uea_H_from_x_als(x(1:nd,:),x((nd+1):end,:));
end

end
