% RGB2LUV performs RGB to Luv color conversion.
%   LUV = RGB2LUV( RGB ) converts the RGB color vectors into LUV. RGB must
%   range in (0.0,1.0).
%
% See also LUV2RGB
%

% 6.891 Problem Set 4: Problem 1
% C. Mario Christoudias

function luv = rgb2luv( rgb )

% convert to XYZ

XYZ = [0.4125,  0.3576,  0.1804; ...
       0.2125,  0.7154,  0.0721; ...
       0.0193,  0.1192,  0.9502];

xyz = XYZ * rgb;

% convert to Luv

luv = xyz;

Yn = 1;
Lt = 0.008856;
Un_prime = 0.19784977571475;
Vn_prime = 0.46834507665248;
L0 = xyz(2,:) / Yn;

warning off MATLAB:divideByZero;
constant = xyz(1,:) + 15 * xyz(2,:) + 3 * xyz(3,:);
u_prime = (constant ~= 0) .* ((4 * xyz(1,:)) ./ constant) + (constant == 0) * 4.0;
v_prime = (constant ~= 0) .* ((9 * xyz(2,:)) ./ constant) + (constant == 0) * 9.0/15.0;

luv(1,:) = (L0 > Lt) .* (116.0 * (L0 .^ (1/3)) - 16.0) + (L0 <= Lt) .* (903.3 * L0);
luv(2,:) = 13 * luv(1,:) .* (u_prime - Un_prime);
luv(3,:) = 13 * luv(1,:) .* (v_prime - Vn_prime);

% be rid of NaNs
luv(find(isnan(luv))) = 0;
