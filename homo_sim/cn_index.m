function MP = cn_index(dbname,gamma,ND,BinLen)
%CN_INDEX indexes images by color using comprehensive normalization
%   matching.
%
%   CN_INDEX(DBNAME,GAMMA) returns a matching percentile (MP) matirx for a
%   dataset.
%
%   Parameters:
%   DBNAME: Name of database
%   GAMMA: display gamma (default: 1)
%   ND: number of dimention (default:2)
%   BinLen: hitogram bin size (default: 16)

%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.

%   References:
%   Finlayson, G.D., Schiele, B. and Crowley, J.L., 1998, June. 
%   Comprehensive colour image normalization. In European conference 
%   on computer vision (pp. 475-490).


addpath('../utility');

if nargin<2, gamma = 1; end
if nargin<3, ND = 2; end
if nargin<4, BinLen = 16; end

%% configuration
mlpath = 'model_cache_cn/'; % model cache path
dcpath = 'db_cache_cn/'; % dataset cache

[pm,rm,objsel] = dbparser(dbname); % load data path
objsel = 1:size(pm,1);

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
        pr = reshape(rim,[],3);
        % comprehensive normalisation
        npr = cn(pr);
        % build histogram
        mh = rgb2hist(npr,binrng);
        save(mlName,'mh');
    end
end

% builde cache dir
qpath = [dcpath,dbname,'/'];
if ~exist(qpath,'dir') mkdir(qpath); end

% query tests

% each query image is tested agains all models
for oo = 1:No % for each object
    oi = objsel(oo);
    for ci = 1:Nc % for each condition
        % test if the test condition exists
        if isempty(pm{oi,ci}), continue; end

        qName = [qpath,num2str(oi),'_',num2str(ci),'.mat'];

        if ~exist(qName,'file')
            oim = imread(pm{oi,ci}.d); % load query image
            oim = im2double(oim).^gamma;
            po = reshape(oim,[],3);
            npo = cn(po); % normalise

            % if query cache does not exist
            qh = rgb2hist(npo,binrng);
            %save(qName,'qh');
        end

        %if ~exist(qName,'file'), error('file not found'); end
        %load(qName,'qh');

        score_t = zeros(No,1); % cache for matching score
        
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
    ohist = histcn(rgb(:,1:ND),binrng{:}); % use RG tuple
end

end

function crgb = cn(rgb)

    % exclude saturated index
    ind = sum(rgb==0 | rgb==1,2)==0;

    nrgb = rgb(ind,:);
    [p q] = size(nrgb);
    for i = 1:10
        % normalisation
        nrgb = nrgb./(ones(p,1)*mean(nrgb));
        nrgb = nrgb./(sum(nrgb,2)*ones(1,3));
    end

    crgb = nrgb;
end

