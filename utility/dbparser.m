function [pm,rm,ObjSel] = dbparser(dbname)
% parser illumination dataset and build a path
% matrix for each testcase (object by condition)

% Copyright 2015 Han Gong
% University of East Anglia
% http://www2.cmp.uea.ac.uk/~ybb15eau

dbpath = ['../data/',dbname,'/']; % dataset path

dbname = strrep(dbname,'/','_');
if exist([dbname,'.mat'])
    load([dbname,'.mat'],'pm','rm','ObjSel');
    return
end

switch dbname
case 'uea_uncalibrated'
    dCond = getAll(dbpath,'d'); % dir of conditions
    nCond = numel(dCond); % number of conditions
    nObj = 28; % we know that there are 28 objects
    pm = cell(nObj,nCond);
    rm = cell(nObj,1); % cell array for reference image paths
    for ci = 1:nCond
        for oi = 1:nObj
            pm{oi,ci}.d = [dbpath,dCond{ci},'/',num2str(oi),'.tif'];
            rm{oi}.d = [dbpath,dCond{2},'/',num2str(oi),'.tif'];
        end
    end

    ObjSel = 1:nObj;
    ObjSel = ObjSel(~ismember(ObjSel,[6,7,8]));
case 'SFU_objects_test'
    dObj = getAll(dbpath,'d'); % dir of objects
    dObj = sort_nat(dObj);
    nObj = numel(dObj); % number of conditions

    ObjSel = 1:nObj;
    ObjSel = ObjSel(~ismember(ObjSel,[1,2,31]));
    nObj = numel(dObj); % number of conditions

    % define pattern of conditions
    Conds = {'ph-ulm','solux-3500','solux3500+3202','solux-4100','solux-4700',...
             'solux-4700+3202','syl-50MR16Q','syl-50MR16Q+3202','syl-cwf','syl-wwf'};
    nCond = numel(Conds);

    pm = cell(nObj,nCond); % cell array for image paths
    rm = cell(nObj,1); % cell array for reference image paths
    refPat = Conds{7}; % refernece pattern

    for oi = 1:nObj 
        for ci = 1:nCond
            file = [dbpath,dObj{oi},'/',dObj{oi},'_',Conds{ci}];
            if exist([file,'.tif'],'file')
                % build image data paths
                pm{oi,ci}.d = [file,'.tif']; % data path
            end
            reffile = [dbpath,dObj{oi},'/',dObj{oi},'_',refPat];
            rm{oi}.d = [reffile,'.tif']; % data path
        end
    end

case {'aloi_col','aloi_ill','aloi_view'}
    dObj = getAll(dbpath,'d'); % dir of objects
    dObj = sort_nat(dObj);
    ObjSel = [2,13,15,18,23,25,36,46,75,76,77,78,85,95,110,111,124,126,...
              136,137,144,154,157,169,177,185,195,196,197,198,199,200,204,205,206,...
              207,208,210,212,218,222,248,259,263,264,274,276,277,284,...
              288,290,296,297,298,301,305,312,313,320,322,325,327,340,349,...
              353,355,362,369,371,381,385,390,391,393,396,420,426,437,447,...
              452,454,464,478,485,489,509,514,526,529,541,547,550,561,...
              562,564,584,588,592,599,602,604,606,608,610,611,615,616,618,...
              619,626,627,630,637,638,640,648,674,676,677,678,681,684,...
              687,688,696,697,698,699,761,725,726,729,730,736,741,743,748,...
              749,750,755,765,766,769,771,772,773,776,777,778,782,789,784,795,796,...
              798,804,805,806,815,822,825,828,831,832,837,841,842,846,853,...
              854,857,858,861,864,868,872,874,876,890,906,907,908,909,910,927,...
              930,945,946,958,959,960,961,962,965,967,969,970,986,989,990,992,994];
    nObj = numel(dObj); % number of conditions

    % assume that nCond is the same for all objects
    nCond = numel(getAll([dbpath,dObj{1}],'f'));

    pm = cell(nObj,nCond); % cell array for image paths
    rm = cell(nObj,1);
    mskpath = '../data/aloi/mask/'; % mask path
    refpath = '../data/aloi/view/';

    for oi = 1:nObj
        % get all cond files of an object
        fCond = getAll([dbpath,dObj{oi}],'f');
        fCond = sort_nat(fCond);
        for ci = 1:numel(fCond)
            [~,fCondName] = fileparts(fCond{ci});
            % regexp for right mask
            [si,ei] = regexp(fCondName,'[c,i,r]{1}\d+');
            cat = fCondName(si); % category
            CondName = fCondName(si:ei); % object name
            if strcmp(cat,'i') || ...
              (strcmp(cat,'r') && strcmp(CondName,'0')) % frontal view for illumination test
                cat = 'c'; CondName = '1';
            end
            % build image data paths
            pm{oi,ci}.d = [dbpath,dObj{oi},'/',fCond{ci}]; % data path
            pm{oi,ci}.m = [mskpath,dObj{oi},'/',dObj{oi},'_',cat,...
                           CondName,'.png']; % mask path
        end

        % build reference image paths (frontal view)
        rm{oi}.d = [refpath,dObj{oi},'/',dObj{oi},'_r0.png'];
        rm{oi}.m = [mskpath,dObj{oi},'/',dObj{oi},'_c1.png'];
    end
otherwise
    error('unkown option');
end

save([dbname,'.mat'],'pm','rm','ObjSel');
