function MP = homo_index(dbname,gamma,check_comp)
%HOMO_INDEX indexes images by color using homographical ASIFT features 
%   matching.
%
%   HOMO_INDEX(DBNAME,GAMMA) returns a matching percentile (MP) matirx for a
%   dataset.
%
%   Parameters:
%   DBNAME: Name of database
%   GAMMA: display gamma (default: 1)
%   CHECK_COMP: check homography compatibility (default: true). The MP 
%   values for imcompatible cases will be NaN.

%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.

%   References:
%   Finlayson, G.D., Gong, H. and Fisher, R., 2017. Color Homography: theory
%   and applications. IEEE Transactions on Pattern Analysis and Machine
%   Intelligence.


addpath('../homo_solver');
addpath('../utility');

if ~exist('gamma','var'), gamma = 1; end
if ~exist('check_comp','var'), check_comp = true; end

%% configuration
mlpath = 'model_cache2/'; % model cache path
dcpath = 'db_cache2/'; % dataset cache
C = [1,0,0;0,1,0;1,1,1]; % base conversion matrix
res = 320;
Nl = 4; % fixed gamma

% ASIFT junk path
c = clock;
c = arrayfun(@num2str,c,'UniformOutput',false);
tpath = [cell2mat(c),'/'];
mkdir(tpath);

[pm,rm] = dbparser(dbname); % load data path
objsel = 1:size(pm,1);

No = size(objsel,2); % number of objects
Nc = size(pm,2); % number of conditions
MP = NaN(No,Nc); % MP for each query
if exist('MP.mat','file'), load('MP.mat','MP'); end % recovery

% build model caches
mpath = [mlpath,dbname];
if ~exist(mpath,'dir'), mkdir(mpath); end

parfor oo = 1:No
    oi = objsel(oo);
    mlName = [mlpath,dbname,'/',num2str(oi),'.txt'];
    if ~exist(mlName,'file')
        rim = imread(rm{oi}.d); % load ref image
        rim = im2double(rim).^gamma;
        pr = reshape(rim,[],3)';
        cr = chrodist(pr,res).^(1/Nl); % build hist
        % extract ASIFT feature and build cache
        cr_path = [tpath,'cr',num2str(oi),'.png'];
        imwrite(cr,cr_path); % write image
        cmd = ['./detectASIFTfeature ', cr_path,' ',mlName,' 0'];
        [status,~] = system(cmd);
        delete(cr_path);
    end
end

% query tests

% builde cache dir
qpath = [dcpath,dbname,'/'];
if ~exist(qpath,'dir'), mkdir(qpath); end

% each query image is tested agains all models
try
for oo = 1:No % for each object
    oi = objsel(oo);
    for ci = 1:Nc % for each condition
        % test if the test condition exists
        if isempty(pm{oi,ci}), continue; end
        % test if MP exists
        if ~isnan(MP(oi,ci)), continue; end
        
        % load query image
        oim = imread(pm{oi,ci}.d);
        oim = im2double(oim).^gamma;
        po = reshape(oim,[],3)';
        % test if the input image is suitable for homography match
        if (~check_comp) || CheckSpd(po)
            co = chrodist(po,res).^(1/Nl); % build hist

            % if query cache does not exist
            % extract ASIFT feature of current query
            qName = [qpath,num2str(oi),'_',num2str(ci),'.txt'];
            if ~exist(qName,'file')
                imwrite(co,[tpath,'co.png']); % write image
                cmd = ['./detectASIFTfeature ',[tpath,'co.png '],qName,' 0'];
                [~,~] = system(cmd);
            end
            
            score_t = zeros(No,1); % cache for matching score
            % go through all models to compare
            parfor pp = 1:No % for each model to compare (parfor)
                pi = objsel(pp);

                % load referencing chromaticity feature profile
                mlName = [mlpath,dbname,'/',num2str(pi),'.txt'];
                if ~exist(mlName,'file'), error('file not found'); end
                fid = fopen(mlName,'r');
                npts_r = fscanf(fid,'%d',1); % number of features
                fclose(fid);

                if ~exist(qName,'file'), error('file not found'); end
                fid = fopen(qName,'r');
                npts_o = fscanf(fid,'%d',1); % number of features
                fclose(fid);

                % try match the features
                mat_cache = [tpath, num2str(oi),'_',num2str(pi),'match.txt'];
                cmd = ['./matchASIFTfeature ',qName,' ',mlName,' ', mat_cache];
                [~,~] = system(cmd);
                fid = fopen(mat_cache,'r');
                npts = fscanf(fid,'%d',1); % number of matches
                if npts>4
                    P = fscanf(fid,'%f');
                    P = reshape(P,[],npts);

                    p1 = P([1,2],:);
                    p2 = P([3,4],:);

                    %make compatible points
                    p1 = cat(1,p1/res,ones(1,npts));
                    p2 = cat(1,p2/res,ones(1,npts));

                    %compute homograhpy between two images
                    rng(0);
                    [inliers] = uea_ransacfithomography(p1,p2,0.01,2000);

                    npts_matched = size(inliers,2);
                else
                    npts_matched = 0;
                end
                fclose(fid);
                delete(mat_cache);

                score_t(pp) = (npts_matched.^2)/(npts_o*npts_r); % matching score
            end

            % get match percentile for this query
            if sum(score_t)>0
                [~,ind] = sort(score_t,'descend');
                MP(oo,ci) = (No-find(ind==oo))/(No-1);
            end
        end

        % if no feature (very unlikely) is found, then the probability is a
        % random selection i.e. MP=0.5.
        if (~check_comp) || isnan(MP(oo,ci)), MP(oo,ci) = 0.5; end
        fprintf('(%d,%d) = %.4f\n',oi,ci,MP(oo,ci));
    end
end
catch e
    save('MP.mat','MP');
    rethrow(e);
end

rmdir(tpath,'s');

rmpath('../utility');
rmpath('../homo_solver');

end

function ret = CheckSpd(rgb)
% returns color structure complexity given RGBs

    % exclude saturated pixels
    ex_msk = rgb<0.01;
    ex_msk = ~logical(sum(ex_msk,1));
    ex_msk = ex_msk & mean(rgb,1)>0.1;

    C = [1,0,0;0,1,0;1,1,1]; % base converse matrix
    pC = C*rgb(:,ex_msk); % chromaticity array
    hpC = bsxfun(@rdivide,pC(1:2,:),pC(3,:)); % homogenous ones

    ret = true;
    [~,d] = eig(cov(hpC')); % check eigen value
    d = diag(d);
    %sqrt(d(1)*d(2))
    if sqrt(d(1)*d(2))<1e-3, ret = false; end
end
