function M = lscal(rgb,xyz)
% LSCAL computes the colour correction matrix by using
% the linear least squares method

M = rgb\xyz;

