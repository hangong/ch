function H = H_from_feat(oim,rim,res,plotflag,Nl)

if nargin<4, plotflag = false; end
if nargin<3, res = 500; end

% dimension check
if ndims(oim)==2
    co = oim;
else
    co = chrodist(reshape(oim,[],3)',res);
end
if ndims(oim)==2
    cr = rim;
else
    cr = chrodist(reshape(rim,[],3)',res);
end

% compute chromaticity distributions for both images
if nargin<5
    Nl = 3; Nl_max = 8;
else
    Nl_max = Nl+1;
end;

npts = 0;
while (npts<4 && Nl<Nl_max)
    Nl = Nl+1;
    Io = co.^(1/Nl);
    Ir = cr.^(1/Nl);
    % match features
    [p1,p2] = distmatchbyfeat(Io,Ir);
    npts = size(p1,2);
end

% compute homography according to matched features
if npts==0
    H = [];
else
    [inliers, H] = uea_ransacfithomography(p1,p2,0.01);
end

% debug information
if plotflag
    C = [1,0,0;0,1,0;1,1,1]; % base conversion matrix
    disp('Estimated H');
    disp(H);

    %% display matches
    figure;
    subplot(1,2,1);
    imshowpair(Io,Ir);
    title('Inital Gamuts');
    fprintf('Matches: %d, Nl: %d\n',npts,Nl);

    %% display gamut match
    if ~isempty(H)
        poim = reshape(oim,[],3)'; % rendered intensities
        poC = C*poim; % chromaticity array
        peC = H*poC; % apply homography
        pe = C\peC; % convert back to rgb
        Ie = chrodist(pe,res).^(1/Nl);

        prim = reshape(rim,[],3)';
        prC = C*prim;

        eim = reshape(pe',size(oim));

        % align exposure
        sc = mean(oim(:))/mean(eim(:));
        eim = eim*sc;
        %eim = eim/max(eim(:));

        % prepare chromaticity image display
        oimC = reshape(poC',size(oim));
        rimC = reshape(prC',size(rim));
        eimC = reshape(peC',size(oim));
        oimC = bsxfun(@rdivide, oimC, oimC(:,:,3));
        rimC = bsxfun(@rdivide, rimC, rimC(:,:,3));
        eimC = bsxfun(@rdivide, eimC, eimC(:,:,3));

        cr = chrodist(reshape(rim,[],3)',res).^(1/Nl);
        subplot(1,2,2);
        imshowpair(Ie,cr);
        title('Aligned Gamuts');
    end

    figure;
    imshow(1-[Io,Ir]);
    if ~isempty(H)
        xx = round([p1(1,:); p2(1,:) + size(Io,2)/res]*res);
        yy = round([p1(2,:); p2(2,:)]*res);
        line(xx(:,inliers), yy(:,inliers),'LineWidth',2);

        gamma = 1;
        figure;
        subplot(2,3,1);
        imshow(oim.^gamma);
        title('original');
        subplot(2,3,2);
        imshow(rim.^gamma);
        title('reference');
        subplot(2,3,3);
        imshow(max(eim,0).^gamma);
        title('estimated');
        subplot(2,3,4); % second row for chromaticity plot
        imshow(oimC);
        title('original');
        subplot(2,3,5);
        imshow(rimC);
        title('reference');
        subplot(2,3,6);
        imshow(eimC);
        title('estimated');
    end
end

end

