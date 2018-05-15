function [p1,p2] = distmatchbyfeat(Io,Ir)

res = size(Io,1);

% SIFT
%{
% detect and match features
[p1,d1] = vl_sift(single(Io));
[p2,d2] = vl_sift(single(Ir));

% match features
[matches, scores] = vl_ubcmatch(d1,d2,2.5) ;
[drop, perm] = sort(scores, 'descend');
matches = matches(:, perm);
scores  = scores(perm);
p1 = p1(1:2,matches(1,:));
p2 = p2(1:2,matches(2,:));
%}

% ASIFT
tpath = 'asiftjunk_dm/';
if exist(tpath,'dir')
    rmdir(tpath,'s');
end
mkdir(tpath);

% write images
imwrite(Io,[tpath,'Io.png']);
imwrite(Ir,[tpath,'Ir.png']);

cmd = ['./detectASIFTfeature ',[tpath,'Io.png '],[tpath,'key1.txt '],'0']; % command
[status,~] = system(cmd);

cmd = ['./detectASIFTfeature ',[tpath,'Ir.png '],[tpath,'key2.txt '],'0']; % command
[status,~] = system(cmd);

cmd = ['./matchASIFTfeature ', [tpath,'key1.txt '],[tpath,'key2.txt '], [tpath,'match.txt ']]; % command
[status,~] = system(cmd);

fid = fopen([tpath,'match.txt'],'r');
npts = fscanf(fid,'%d',1);

if npts>0
    P = fscanf(fid,'%f');
    P = reshape(P,[],npts);

    p1 = P([1,2],:);
    p2 = P([3,4],:);

    % make compatible points
    p1 = cat(1,p1/res,ones(1,npts));
    p2 = cat(1,p2/res,ones(1,npts));
else
    p1 = []; p2 = [];
end
 
fclose(fid);
