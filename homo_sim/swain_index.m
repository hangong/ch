function MP = swain_index(dbname,gamma,ND,BinLen)
%SWAIN_INDEX indexes images by color using rg chromaticity and histogram
%   intersection.
%
%   SWAIN_INDEX(DBNAME,GAMMA) returns a matching percentile (MP) matirx 
%   for a dataset.
%
%   Parameters:
%   DBNAME: Name of database
%   GAMMA: display gamma (default: 1)
%   ND: number of dimention (default:2)
%   BinLen: hitogram bin size (default: 16)

%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.

%   References:
%   Swain, M.J. and Ballard, D.H., 1991. Color indexing. International 
%   journal of computer vision, 7(1), pp.11-32.

addpath('../utility');

if nargin<2, gamma = 1; end
if nargin<3, ND = 2; end
if nargin<4, BinLen = 16; end

%% configuration
mlpath = 'model_cache_swain/'; % model cache path
dcpath = 'db_cache_swain/'; % dataset cache

[pm,rm,objsel] = dbparser(dbname); % load data path
%objsel = 1:size(pm,1);

No = size(objsel,2); % number of objects
Nc = size(pm,2); % number of conditions
MP = NaN(No,Nc); % MP for each query

binrng = cell(1,ND);
for c = 1:ND, binrng{c} = linspace(0,1,BinLen+1); end

% build model caches
mpath = [mlpath,dbname];
if ~exist(mpath,'dir'), mkdir(mpath); end

for oi = objsel
    mlName = [mlpath,dbname,'/',num2str(oi),'.mat'];
    if ~exist(mlName,'file')
        rim = imread(rm{oi}.d); % load ref image
        rim = im2double(rim).^gamma; % de-gamma
        pr = reshape(rim,[],3)';
        % build histogram
        mh = rgb2hist(pr,binrng);
        save(mlName,'mh');
    end
end

% builde cache dir
qpath = [dcpath,dbname,'/'];
%if ~exist(qpath,'dir') mkdir(qpath); end

% query tests

% each query image is tested agains all models
for oo = 1:No % for each object
    oi = objsel(oo);
    for ci = 1:Nc % for each condition
        % test if the test condition exists
        if isempty(pm{oi,ci}), continue; end

        oim = imread(pm{oi,ci}.d); % load query image
        oim = im2double(oim).^gamma;
        po = reshape(oim,[],3)';

        % if query cache does not exist
        qName = [qpath,num2str(oi),'_',num2str(ci),'.mat'];
        qh = rgb2hist(po,binrng);
        %save(qName,'qh');

        score_t = zeros(No,1); % cache for matching score

        %if ~exist(qName,'file'), error('file not found'); end
        %load(qName,'qh');

        % go through all models to compare
        parfor pp = 1:No % for each model to compare
            
            pi = objsel(pp);
            % load referencing chromaticity feature profile
            mlName = [mlpath,dbname,'/',num2str(pi),'.mat'];
            if ~exist(mlName,'file'), error('file not found'); end
            l = load(mlName,'mh');       

            % histogram intersection
            score_t(pp) = sum(min(l.mh(:),qh(:)))/sum(l.mh(:));
        end

        % get match percentile for this query
        [~,ind] = sort(score_t,'descend');
        MP(oo,ci) = (No-find(ind==oo))/(No-1);
        fprintf('(%d,%d) = %.4f\n',oi,ci,MP(oo,ci));
    end
end

rmpath('../utility');

function ohist = rgb2hist(rgb,binrng)

    % exclude saturated index
    ind = sum(rgb==0 | rgb==1,1)==0;
    nrgb = rgb(:,ind);

    sumrgb = sum(nrgb);

    opc = bsxfun(@rdivide, nrgb(1:ND,:), sumrgb);
    ohist = histcn(opc',binrng{:});
end

end
