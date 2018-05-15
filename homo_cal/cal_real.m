%   This is the main script for color correction evaluation. The test is 
%   based on a 24-patch color checker.
%
%   Copyright 2018 Han Gong <gong@fedoraproject.org>, University of East
%   Anglia.

%   References:
%   Finlayson, G.D., Gong, H. and Fisher, R., 2017. Color Homography: theory
%   and applications. IEEE Transactions on Pattern Analysis and Machine
%   Intelligence.

addpath('../utility');

% configration
dbpath = '../data/HG_ColourChecker/'; % path of rawdata

fmethod = {@alshomocal;@ransachomocal_luv;@lscal;...
           @(p1,p2) alsRPcal(p1,p2,[4,6])};
Method = {'ALS';'ALS_RANSAC';'LS';'ALS_RP'};

% discover a list of images for conversion
fl = getAll([dbpath,'patch_real'],'f'); % get all files
fn = sort_nat(fl);

Npic = numel(fn);
Npatch = 24;
Nmethod = numel(fmethod);

% non-uniform shading errors
de00_n = zeros(Npatch,Npic,Nmethod);
de76_n = zeros(Npatch,Npic,Nmethod);
deluv_n = zeros(Npatch,Npic,Nmethod);
dergb_n = zeros(Npatch,Npic,Nmethod);

% uniform shading errors
de76_u = zeros(Npatch,Npic,Nmethod);
deluv_u = zeros(Npatch,Npic,Nmethod);
dergb_u = zeros(Npatch,Npic,Nmethod);

% relative difference of correction matrix
md = zeros(1,Npic,Nmethod);

for i = 1:Npic
    
    % ref cat
    cat = regexp(fn{i},'^[^_]+','match');
    % load data
    load([dbpath,'patch_real/',fn{i}]);
    % load reference data
    load([dbpath,'ref_real-',cat{1},'.mat']);

    xyz_std = ref.XYZ./ref.XYZ(4,2); % refernece XYZ
    lab_ref = HGxyz2lab(xyz_std,xyz_std(4,:)); % reference LAB
    luv_ref = HGxyz2luv(xyz_std,xyz_std(4,:)); % reference LUV
    rgb_ref = xyz2rgb(xyz_std,'WhitePoint',xyz_std(4,:)); % reference LUV

    % compute colour correction matrix
    fsv = reshape(cap.sv,[],3);
    fsv_uniform = reshape(cap.sv_uniform,[],3);
    for m = 1:Nmethod
        % compute the color correction transform
        switch Method{m}
        case {'ALS_RANSAC'}
            M_n = fmethod{m}(fsv,xyz_std,xyz_std(4,:),fsv_uniform);
            M_u = fmethod{m}(fsv_uniform,xyz_std,xyz_std(4,:),fsv_uniform);
        otherwise
            M_n = fmethod{m}(fsv,xyz_std);
            M_u = fmethod{m}(fsv_uniform,xyz_std);
        end

        % compute xyz using the ground truth RGBs
        switch Method{m}
        case {'ALS','ALS_RANSAC'}
            xyz_est_n = uea_homocvt(fsv_uniform,M_n);
            xyz_est_u = uea_homocvt(fsv_uniform,M_u);
        case {'ALS_RP'}
            xyz_est_n = M_n.cfun(fsv_uniform',M_n.matrix,M_n.terms)';
            xyz_est_u = M_u.cfun(fsv_uniform',M_u.matrix,M_u.terms)';
        otherwise
            xyz_est_n = fsv_uniform*M_n;
            xyz_est_u = fsv_uniform*M_u;
        end

        % normalize by a white patch's green intensity
        %XYZ_est_n = xyz_est_n./xyz_est_u(4,2);
        %XYZ_est_u = xyz_est_u./xyz_est_u(4,2);
        XYZ_est_n = xyz_est_n./xyz_est_n(4,2);
        XYZ_est_u = xyz_est_u./xyz_est_u(4,2);
        
        % Evaluation
        % non-uniform test DE LAB
        lab_est_n = HGxyz2lab(XYZ_est_n,xyz_std(4,:));
        de76_n(:,i,m) = deltaE1976(lab_ref,lab_est_n);
        % uniform test DE LAB
        lab_est_u = HGxyz2lab(XYZ_est_u,xyz_std(4,:));
        de76_u(:,i,m) = deltaE1976(lab_ref,lab_est_u);

        % LUV error
        luv_est_n = HGxyz2luv(XYZ_est_n,xyz_std(4,:));
        luv_est_u = HGxyz2luv(XYZ_est_u,xyz_std(4,:));
        deluv_n(:,i,m) = deltaE1976(luv_ref,luv_est_n);
        deluv_u(:,i,m) = deltaE1976(luv_ref,luv_est_u);
        
        % RGB error
        rgb_est_n = xyz2rgb(XYZ_est_n,'WhitePoint',xyz_std(4,:));
        rgb_est_u = xyz2rgb(XYZ_est_u,'WhitePoint',xyz_std(4,:));
        dergb_n(:,i,m) = deltaE1976(rgb_ref,rgb_est_n);
        dergb_u(:,i,m) = deltaE1976(rgb_ref,rgb_est_u);
    end
end

% print evaluation results (non-uniform shading)
trgb_n = gentab(dergb_n,Method,'RGB (Non-Uniform)');
t76_n = gentab(de76_n,Method,'DeltaE LAB 1976 (Non-Uniform)');
tluv_n = gentab(deluv_n,Method,'DeltaE LUV (Non-Uniform)');

rmpath('../utility');
