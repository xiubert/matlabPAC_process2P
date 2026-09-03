function outFile = plotROIoverlay(dataPath,animal,opts)
% plotROIoverlay  Render an ROI set over its mean image, for review.
%
%   outFile = plotROIoverlay(dataPath,animal)
%   outFile = plotROIoverlay(dataPath,animal,'cond','all','outDir',D)
%
%   A headless run has no moment where a human sees the ROIs, so it has to
%   leave one behind. This draws every ROI's outline over the image it was
%   segmented from, labelled by ID, and — when the set came from
%   consensusROIsets — coloured by how many tifs detected each cell, which is
%   the fastest way to see whether the marginal ROIs are real.
%
%   Works retrospectively: everything it needs is inside the saved ROI bundle
%   and the session mean tif, so it can be run over a finished cohort without
%   re-segmenting.
%
%   Inputs
%     dataPath  animal data folder holding <animal>_moCorrROI_<cond>.mat.
%     animal    animal ID.
%
%   Name-value
%     'cond'     condition to draw. Default: every condition found.
%     'meanImg'  H x W image to draw on. Default: NoRMCorred/<animal>_
%                cellposeMean.tif when present, else the condition's mean
%                projection via condMeanImg, else a blank field.
%     'outDir'   where the PNG goes. Default: the artifact dir.
%     'artifactDir'  where the ROI files and session mean live. Default '' =
%                dataPath (the flat legacy layout). See animalPaths.
%     'showIDs'  label each ROI with its ID. Default true.
%     'resolution' export DPI. Default 150.
%     'visible'  show the figure. Default false (it is closed after export).
%
%   Output
%     outFile   cellstr of the PNGs written, one per condition.
%
%   See also cellposeROIset, consensusROIsets, labelImg2moCorROI, condMeanImg

arguments
    dataPath (1,:) char
    animal   (1,:) char
    opts.cond       (1,:) char = ''
    opts.meanImg          double = []
    opts.outDir     (1,:) char = ''
    opts.showIDs    (1,1) logical = true
    opts.resolution (1,1) double = 150
    opts.visible    (1,1) logical = false
    opts.artifactDir(1,:) char = ''
end

artDir = opts.artifactDir;
if isempty(artDir); artDir = dataPath; end
outDir = opts.outDir;
if isempty(outDir); outDir = artDir; end
if ~isfolder(outDir); mkdir(outDir); end

if isempty(opts.cond)
    d = dir(fullfile(artDir,[animal '_moCorrROI_*.mat']));
    conds = regexprep({d.name},['^' animal '_moCorrROI_(.*)\.mat$'],'$1');
else
    conds = {opts.cond};
end
if isempty(conds)
    error('plotROIoverlay:noROIfiles','No %s_moCorrROI_*.mat in %s',animal,dataPath);
end

outFile = {};
for c = 1:numel(conds)
    cond = conds{c};
    f = fullfile(artDir,[animal '_moCorrROI_' cond '.mat']);
    S = load(f);
    if ~isfield(S,'moCorROI') || isempty(S.moCorROI)
        warning('plotROIoverlay:noROI','%s holds no ROIs; skipping.',f);
        continue
    end
    roi = S.moCorROI;
    roi = roi(~[roi.deleted]);
    sz  = size(roi(1).mask);

    img = resolveImage(opts.meanImg,dataPath,artDir,animal,cond,sz);

    fh = figure('Visible',onOff(opts.visible),'Color','w',...
        'Position',[100 100 900 780]);
    %drawn as explicit RGB, not imagesc: the vote colorbar below needs its own
    %clim, and setting clim on a scalar-mapped image would rescale the anatomy
    %to the vote range and wash it out
    lo = prctile(img(:),1); hi = prctile(img(:),99.5);
    g  = min(max((img-lo)./max(hi-lo,eps),0),1);
    image(repmat(g,1,1,3)); axis image off; hold on

    hasVotes = isfield(roi,'votes') && ~isempty(roi(1).votes);
    if hasVotes
        v  = [roi.votes];
        nS = max(v);
        if isfield(S,'consensusInfo') && isfield(S.consensusInfo,'nSets') ...
                && ~isempty(S.consensusInfo.nSets)
            nS = S.consensusInfo.nSets;
        end
        cmap = parula(256);
        %colour by detection count: a cell seen in one or two tifs is the one
        %worth a second look, and it should not look like a cell seen in all
        cIdx = @(k) max(1,min(256,round(255*(v(k)-1)/max(1,nS-1))+1));
    end

    for k = 1:numel(roi)
        xy = roi(k).ROIcurveOrderedXY;
        if hasVotes; col = cmap(cIdx(k),:); else; col = [1 1 0]; end
        plot(xy(1,:),xy(2,:),'-','Color',col,'LineWidth',1.1);
        if opts.showIDs
            text(roi(k).pos(1),roi(k).pos(2)-3,roi(k).ID,'Color',col,...
                'FontSize',6,'FontWeight','bold','Clipping','on');
        end
    end

    ttl = sprintf('%s  |  %s  |  %d ROIs',animal,cond,numel(roi));
    if hasVotes
        ttl = sprintf('%s  |  votes %d-%d of %d (median %g)',...
            ttl,min(v),max(v),nS,median(v));
        colormap(gca,cmap);
        clim([1 max(2,nS)]);          % safe: the image is truecolor
        cb = colorbar('eastoutside');
        cb.Label.String = 'tifs detecting this cell';
    end
    if isfield(S,'roiParams') && isfield(S.roiParams,'dilatePx')
        ttl = sprintf('%s  |  dilatePx %d',ttl,S.roiParams.dilatePx);
    end
    title(ttl,'Interpreter','none','FontSize',10);

    outFile{end+1} = fullfile(outDir,sprintf('%s_ROIoverlay_%s.png',animal,cond)); %#ok<AGROW>
    exportgraphics(fh,outFile{end},'Resolution',opts.resolution);
    close(fh);
end
end

% ------------------------------------------------------------------------
function img = resolveImage(given,dataPath,artDir,animal,cond,sz)
if ~isempty(given)
    img = given;
    return
end
%the segmented session mean, written by cellposeROIset
%written next to the ROI files by cellposeROIset; older runs left it under
%NoRMCorred/, so look there too
meanTif = fullfile(artDir,[animal '_cellposeMean.tif']);
if ~isfile(meanTif)
    meanTif = fullfile(dataPath,'NoRMCorred',[animal '_cellposeMean.tif']);
end
if isfile(meanTif)
    try
        img = double(imread(meanTif));
        if isequal(size(img),sz); return, end
    catch
    end
end
%otherwise project this condition off the motion-corrected tifs
try
    S = load(fullfile(dataPath,[animal '_tifCondSplitLegend.mat']),'tifList');
    if isfield(S.tifList,cond)
        img = condMeanImg(dataPath,S.tifList.(cond));
        if isequal(size(img),sz); return, end
    end
catch
end
warning('plotROIoverlay:noMeanImage',...
    ['No mean image for %s/%s; drawing the outlines on a blank field. The '...
     'shapes are still reviewable, the anatomy is not.'],animal,cond);
img = zeros(sz);
end
% ------------------------------------------------------------------------
function s = onOff(tf)
if tf; s = 'on'; else; s = 'off'; end
end
