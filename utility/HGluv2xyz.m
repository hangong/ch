function [xyz] = HGluv2xyz(luv,white)
% white - white xyz

% http://www.brucelindbloom.com/index.html?Eqn_Luv_to_XYZ.html

if (size(luv,2)~=3)
   disp('luv must be n by 3'); return;   
end
xyz = zeros(size(luv));

% compute u' v' for white
upw = 4*white(1)/(white(1) + 15*white(2) + 3*white(3));
vpw = 9*white(2)/(white(1) + 15*white(2) + 3*white(3));

% compute Y
index = (luv(:,1) > 8);
xyz(:,2) = xyz(:,2) + index.*(((luv(:,1)+16)/116).^3);
xyz(:,2) = xyz(:,2) + (1-index).*(luv(:,1)./903.3);

% compute a, b, c, d
a = (1/3)*(52*luv(:,1)./(luv(:,2)+13*luv(:,1)*upw) - 1);
b = -5*xyz(:,2);
c = -1/3;
d = xyz(:,2).*(39*luv(:,1)./(luv(:,3)+13*luv(:,1)*vpw)-5);

% compute X,Z
xyz(:,1) = (d-b)./(a-c);
xyz(:,3) = xyz(:,1).*a+b;

end

