function tab = gentab(de,Method,de_str)
%GENTAB generates an error result statistics table.
%
%   Copyright 2018 Han Gong <gong@fedoraproject.org>
%   University of East Anglia

if nargin<3, disp_flag = false; else, disp_flag = true; end

if length(Method) ~= size(de,3)
    error('Dimension of methods dismatches the error matrix.');
end

summary = @(f,x) squeeze(mean(f(x,1),2));

de_mean = summary(@mean,de);
de_median = summary(@median,de);
de_95 = summary(@(x,d) quantile(x,.95,d),de);
de_max = summary(@(x,d) max(x,[],d),de);

Variable = {'Mean','Median','pct95','Max'};
tab = table(de_mean,de_median,de_95,de_max,...
    'RowNames',Method,'VariableNames',Variable);

% display table
if disp_flag
    disp(de_str);
    disp(tab);
end

end
