function out = cellposeROIset(dataPath,animal,tifList,tifFiles,opts)
% cellposeROIset  Segment with Cellpose and write one _moCorrROI_<cond>.mat per condition.
%
%   out = cellposeROIset(dataPath,animal,tifList,tifFiles)
%   out = cellposeROIset(...,'mode','consensus','minVotes',2,'dilatePx',2)
%
%   The headless replacement for processAnimal2P §4-5 (TIFcatROIgui + save).
%   Two things have to be true of the result: the masks must be cells, and an
%   ID must mean the same cell in every treatment condition. This gets the
%   first from Cellpose and the second by segmenting ONE shared set in a
%   common reference frame and placing it into each condition, so
%   intersectROIfiles (§6) becomes an invariant check rather than the thing
%   doing the matching.
%
%   Two modes, both ending in the same shared-set placement:
%     'consensus' (default) -- segment EVERY tif's mean projection separately
%         and keep the cells enough tifs agree on (consensusROIsets). Cellpose
%         detection on this data is unstable tif to tif, and averaging the
%         tifs first destroys the texture it depends on, so voting across
%         per-tif segmentations recovers more cells than any single pass:
%         on TO0003, 18/18 hand-drawn cells at minVotes 2 versus 16/18 for one
%         pass over the session mean. Every cell also gets a detection count.
%     'sessionMean' -- one segmentation of the frame-weighted mean of all
%         tifs. Cheaper (one GPU call), no vote information.
%
%   Inputs
%     dataPath   animal data folder; holds NoRMCorred/ and receives the files.
%     animal     animal ID, for the filenames.
%     tifList    condition -> tif struct map from §2. Pass only conditions
%                that share a frame size: a 256x128 crop condition is derived
%                from a 256x256 one afterwards via remapROIfile, not here.
%     tifFiles   the full tif list, for tifIDXinAllTifList.
%
%   Name-value -- what to segment
%     'mode'          'consensus' (default) | 'sessionMean'
%     'minVotes'      consensus only. Below 1 it is a FRACTION of the tif
%                     count; 1 or more an absolute count. Default 2: a cell
%                     must appear in at least two tifs. A strict inner join
%                     over every tif yields nothing, because the worst tif
%                     governs.
%     'conds'         restrict to these conditions. Default: all in tifList.
%     'tifMeans'      struct of precomputed per-tif means (one cell array per
%                     condition), to skip re-reading NoRMCorred/. Optional.
%
%   Name-value -- ROI construction (see labelImg2moCorROI)
%     'dilatePx'      Default 2, calibrated on TO0003 so segmented masks come
%                     to ~0.86 of the hand-drawn area and dF/F amplitudes stay
%                     comparable to the existing corpus.
%     'minAreaPx'     Default 20.  'maxAreaPx' Default Inf.
%     'edgeMarginPx'  Default 3, so FISSA's neuropil annulus is not clipped.
%   Name-value -- registration
%     'maxShiftPx'    Default 15 (NoRMCorre's max_shift).
%     'refCond'       Default: the condition with the most tifs.
%   Name-value -- other
%     'cellposeArgs'  name-value pairs forwarded to cellposeSegment. The
%                     defaults here set diameter 15 and cellprob -1, which is
%                     what the TO0003 calibration settled on; cpsam's
%                     automatic sizing fails outright on a smoothed mean.
%     'saveMeanTif'   write <animal>_cellposeMean.tif for QC. Default true.
%     'verbose'       Default true.
%
%   Output (struct)
%     .mode .conds .refCond .shifts .sessionMean .moCorROI (reference frame)
%     .roiInfo .consensusInfo .cellposeParams .files .perCond .meanTif
%
%   Each bundle carries moCorROI, moCorSeqN, nTifs, tifIDXinAllTifList and
%   moCorTifNames (the FISSA driver errors without the last), plus
%   cellposeParams and roiParams for provenance.
%
%   See also cellposeSegment, consensusROIsets, labelImg2moCorROI,
%   registerConditionMeans, remapROIfile, intersectROIfiles

arguments
    dataPath   (1,:) char
    animal     (1,:) char
    tifList    (1,1) struct
    tifFiles         struct
    opts.mode         (1,:) char {mustBeMember(opts.mode,{'consensus','sessionMean'})} = 'consensus'
    opts.minVotes     (1,1) double {mustBePositive} = 2
    opts.conds              cell = {}
    opts.tifMeans           struct = struct()
    opts.dilatePx     (1,1) double {mustBeNonnegative,mustBeInteger} = 2
    opts.minAreaPx    (1,1) double {mustBeNonnegative} = 20
    opts.maxAreaPx    (1,1) double {mustBePositive}    = Inf
    opts.edgeMarginPx (1,1) double {mustBeNonnegative,mustBeInteger} = 3
    opts.maxShiftPx   (1,1) double {mustBeNonnegative} = 15
    opts.refCond      (1,:) char = ''
    opts.cellposeArgs       cell = {'diameter',15,'cellprobThreshold',-1}
    opts.saveMeanTif  (1,1) logical = true
    opts.verbose      (1,1) logical = true
end

if ~isfolder(dataPath)
    error('cellposeROIset:noDataPath','Not a folder: %s',dataPath);
end

allCond = fieldnames(tifList);
conds = opts.conds;
if isempty(conds)
    conds = allCond;
else
    conds = conds(:);
    missing = conds(~ismember(conds,allCond));
    if ~isempty(missing)
        error('cellposeROIset:condNotInTifList',...
            'Condition(s) %s are not in tifList.',strjoin(missing',', '));
    end
end
if isempty(conds)
    error('cellposeROIset:noConditions','No conditions to segment.');
end

%--- per-tif mean projections --------------------------------------------
tifMeans = cell(numel(conds),1);
tifFrames = cell(numel(conds),1);
for c = 1:numel(conds)
    cond = conds{c};
    if isfield(opts.tifMeans,cond)
        tifMeans{c} = opts.tifMeans.(cond)(:)';
        tifFrames{c} = ones(1,numel(tifMeans{c}));
    else
        if opts.verbose
            fprintf('cellposeROIset: projecting %d tifs for ''%s''...\n',...
                numel(tifList.(cond)),cond);
        end
        [tifMeans{c},tifFrames{c}] = projectCondition(dataPath,tifList.(cond));
    end
end

%--- reference condition + registration ----------------------------------
if isempty(opts.refCond)
    [~,refIdx] = max(cellfun(@(c) numel(tifList.(c)),conds));
else
    refIdx = find(strcmp(conds,opts.refCond),1);
    if isempty(refIdx)
        error('cellposeROIset:badRefCond',...
            'refCond ''%s'' is not one of: %s',opts.refCond,strjoin(conds',', '));
    end
end

condMeans = cell(numel(conds),1);
for c = 1:numel(conds)
    condMeans{c} = weightedMean(tifMeans{c},tifFrames{c});
end

if isscalar(condMeans)
    shifts      = [0 0];
    sessionMean = condMeans{1};
    regInfo     = struct('refIdx',1,'shiftMagPx',0,'method','single condition',...
        'names',{conds'},'maxShiftPx',opts.maxShiftPx);
    if opts.verbose
        fprintf('cellposeROIset: one condition (''%s''), no registration needed.\n',conds{1});
    end
else
    [shifts,sessionMean,regInfo] = registerConditionMeans(condMeans,...
        'maxShiftPx',opts.maxShiftPx,'refIdx',refIdx,'names',conds',...
        'verbose',opts.verbose);
end

%--- segment -------------------------------------------------------------
consensusInfo = struct();
switch opts.mode
    case 'sessionMean'
        segImg = fillNaN(sessionMean);
        cpArgs = buildCPargs(animal,'cellposeMean',tifList,conds{refIdx},opts);
        [labelImg,cellposeParams] = cellposeSegment(segImg,cpArgs{:});
        if opts.verbose
            fprintf('cellposeROIset: %d labels from the session mean (%.1f s)\n',...
                cellposeParams.nLabels,cellposeParams.runSeconds);
        end
        [roiRef,roiInfo] = labelImg2moCorROI(labelImg,...
            'dilatePx',opts.dilatePx,'minAreaPx',opts.minAreaPx,...
            'maxAreaPx',opts.maxAreaPx,'edgeMarginPx',opts.edgeMarginPx);

    case 'consensus'
        %every tif, shifted into the reference frame first so one vote map
        %serves all conditions
        roiSets = {};
        nSeg = 0;
        tSeg = tic;
        cellposeParams = struct();
        for c = 1:numel(conds)
            for k = 1:numel(tifMeans{c})
                img = tifMeans{c}{k};
                if any(shifts(c,:) ~= 0)
                    img = imtranslate(img,shifts(c,:),'FillValues',NaN);
                end
                img = fillNaN(img);
                cpArgs = buildCPargs(animal,sprintf('%s_t%03d',conds{c},k),...
                    tifList,conds{c},opts);
                [L,P] = cellposeSegment(img,cpArgs{:});
                nSeg = nSeg + 1;
                if nSeg == 1, cellposeParams = P; end
                if P.nLabels > 0
                    roiSets{end+1} = labelImg2moCorROI(L,...
                        'dilatePx',opts.dilatePx,'minAreaPx',opts.minAreaPx,...
                        'maxAreaPx',opts.maxAreaPx,...
                        'edgeMarginPx',opts.edgeMarginPx); %#ok<AGROW>
                else
                    roiSets{end+1} = labelImg2moCorROI(zeros(size(img))); %#ok<AGROW>
                end
            end
        end
        cellposeParams.nSegmentations = nSeg;
        cellposeParams.totalSeconds   = toc(tSeg);
        if opts.verbose
            fprintf('cellposeROIset: segmented %d tifs in %.0f s (%d..%d labels each)\n',...
                nSeg,cellposeParams.totalSeconds,...
                min(cellfun(@numel,roiSets)),max(cellfun(@numel,roiSets)));
        end

        [roiRef,consensusInfo] = consensusROIsets(roiSets,...
            'minVotes',opts.minVotes,'minAreaPx',opts.minAreaPx,...
            'verbose',opts.verbose);
        roiInfo = struct('nLabels',consensusInfo.nClusters,'nKept',numel(roiRef),...
            'dropped',struct('empty',0,'minArea',consensusInfo.droppedArea,...
                             'maxArea',0,'edge',0),...
            'areaPx',arrayfun(@(r) nnz(r.mask),roiRef),'droppedIDs',[]);
        if ~isempty(roiInfo.areaPx)
            roiInfo.equivDiamPx = 2*sqrt(roiInfo.areaPx/pi);
        else
            roiInfo.equivDiamPx = [];
        end
        labelImg = roiSetToLabel(roiRef,size(sessionMean));
end

if isempty(roiRef)
    error('cellposeROIset:noROI',...
        ['Segmentation produced no usable ROIs. Check the projections and the '...
         'area/edge/vote settings (mode=%s).'],opts.mode);
end
if opts.verbose
    fprintf('cellposeROIset: %d ROIs in the shared set (median %.1f px diameter)\n',...
        numel(roiRef),median(roiInfo.equivDiamPx));
end

roiParams = struct('mode',opts.mode,'minVotes',opts.minVotes,...
    'dilatePx',opts.dilatePx,'minAreaPx',opts.minAreaPx,...
    'maxAreaPx',opts.maxAreaPx,'edgeMarginPx',opts.edgeMarginPx,...
    'registration',regInfo,'shifts',shifts,'conds',{conds'},...
    'refCond',conds{refIdx},'source','cellposeROIset');

%--- QC artefact ----------------------------------------------------------
meanTifPath = '';
if opts.saveMeanTif
    srcTif = refTif(tifList,conds{refIdx},dataPath);
    meanTifPath = fullfile(dataPath,[animal '_cellposeMean.tif']);
    try
        if isfile(srcTif)
            writeTifWithHeader(fillNaN(sessionMean),meanTifPath,srcTif,...
                'class','int16','copyright','cellposeROIset_sessionMean');
        end
    catch ME
        warning('cellposeROIset:meanTifFailed',...
            'Could not write %s (%s); continuing.',meanTifPath,ME.message);
        meanTifPath = '';
    end
end

%--- place the shared set into every condition ---------------------------
files   = cell(numel(conds),1);
perCond = struct('name',{},'nROI',{},'nDropped',{});
for c = 1:numel(conds)
    cond = conds{c};
    roiK = shiftROIset(roiRef,-shifts(c,:),size(sessionMean),opts.edgeMarginPx);
    if isempty(roiK)
        error('cellposeROIset:allShiftedOut',...
            'Every ROI fell outside the frame for condition ''%s''.',cond);
    end

    S = struct();
    S.moCorROI           = roiK;
    S.moCorSeqN          = find(strcmp(allCond,cond),1);
    S.nTifs              = numel(tifList.(cond));
    S.tifIDXinAllTifList = ismember({tifFiles.name}',{tifList.(cond).name}');
    S.moCorTifNames      = strrep({tifList.(cond).name}','.tif','_NoRMCorre.tif');
    S.cellposeParams     = cellposeParams;
    S.roiParams          = roiParams;
    S.consensusInfo      = consensusInfo;
    S.cellposeLabelImg   = uint16(labelImg);

    files{c} = fullfile(dataPath,[animal '_moCorrROI_' cond '.mat']);
    save(files{c},'-struct','S');

    perCond(c).name     = cond;
    perCond(c).nROI     = numel(roiK);
    perCond(c).nDropped = numel(roiRef) - numel(roiK);
    if opts.verbose
        fprintf('  %-24s %3d ROIs -> %s\n',cond,numel(roiK),...
            [animal '_moCorrROI_' cond '.mat']);
    end
end

out = struct('mode',opts.mode,'conds',{conds'},'refCond',conds{refIdx},...
    'shifts',shifts,'sessionMean',sessionMean,'labelImg',labelImg,...
    'moCorROI',roiRef,'roiInfo',roiInfo,'roiParams',roiParams,...
    'consensusInfo',consensusInfo,'cellposeParams',cellposeParams,...
    'files',{files},'perCond',perCond,'meanTif',meanTifPath);
end

% ------------------------------------------------------------------------
function [means,frames] = projectCondition(dataPath,condTifs)
means  = cell(1,numel(condTifs));
frames = zeros(1,numel(condTifs));
for k = 1:numel(condTifs)
    f = fullfile(dataPath,'NoRMCorred',...
        strrep(condTifs(k).name,'.tif','_NoRMCorre.tif'));
    if ~isfile(f)
        error('cellposeROIset:missingTif','Not found: %s',f);
    end
    m = flattenTif(f,'write',false);
    if isstruct(m)
        error('cellposeROIset:multiChannel',...
            '%s still holds %d channels; split them first.',f,numel(m));
    end
    means{k}  = double(m);
    frames(k) = numel(imfinfo(f));
end
end
% ------------------------------------------------------------------------
function m = weightedMean(means,frames)
acc = zeros(size(means{1}));
for k = 1:numel(means)
    acc = acc + means{k}.*frames(k);
end
m = acc./sum(frames);
end
% ------------------------------------------------------------------------
function a = fillNaN(a)
%Cellpose cannot take NaN; the only NaNs here are the thin border a shifted
%condition could not contribute to, and the background median keeps that
%border from reading as an object
bad = ~isfinite(a);
if any(bad(:))
    a(bad) = median(a(~bad));
end
end
% ------------------------------------------------------------------------
function t = refTif(tifList,cond,dataPath)
t = fullfile(dataPath,'NoRMCorred',...
    strrep(tifList.(cond)(1).name,'.tif','_NoRMCorre.tif'));
end
% ------------------------------------------------------------------------
function args = buildCPargs(animal,tag,tifList,cond,opts)
args = [{'name',sprintf('%s_%s',animal,tag)},opts.cellposeArgs];
if ~any(strcmp(args(1:2:end),'srcTif'))
    src = fullfile(tifList.(cond)(1).folder,tifList.(cond)(1).name);
    if isfile(src)
        args = [args,{'srcTif',src}];
    end
end
end
% ------------------------------------------------------------------------
function L = roiSetToLabel(roi,sz)
L = zeros(sz);
for k = 1:numel(roi)
    L(roi(k).mask) = k;
end
end
% ------------------------------------------------------------------------
function roiOut = shiftROIset(roiIn,dxdy,imgSize,edgeMarginPx)
%translate a whole ROI set by an integer [dx dy], dropping any ROI that would
%no longer sit fully inside the frame. Mirrors remapROItoAcq's regeneration of
%the derived geometry fields so FISSA gets a well-formed polygon either way.
dx = round(dxdy(1));
dy = round(dxdy(2));
H  = imgSize(1);
W  = imgSize(2);

if dx == 0 && dy == 0
    roiOut = roiIn;
    return
end

keep = false(1,numel(roiIn));
out  = roiIn;
for k = 1:numel(roiIn)
    m = imtranslate(roiIn(k).mask,[dx dy],'FillValues',0);
    if ~any(m(:)) || nnz(m) ~= nnz(roiIn(k).mask)
        continue    %a mask clipped by the translation is not the same cell
    end
    rws = find(any(m,2));
    cls = find(any(m,1));
    if rws(1) <= edgeMarginPx || rws(end) > H-edgeMarginPx || ...
       cls(1) <= edgeMarginPx || cls(end) > W-edgeMarginPx
        continue
    end
    keep(k) = true;

    r = out(k);
    r.mask              = m;
    r.XYvertices        = [r.XYvertices(:,1)+dx, r.XYvertices(:,2)+dy];
    r.pos(1)            = r.pos(1) + dx;
    r.pos(2)            = r.pos(2) + dy;
    r.label.Position(1) = r.label.Position(1) + dx;
    r.label.Position(2) = r.label.Position(2) + dy;
    r.ROIxyCoord        = mask2polyCoord(m);
    r.ROIcurveOrderedXY = [r.ROIcurveOrderedXY(1,:)+dx; r.ROIcurveOrderedXY(2,:)+dy];
    out(k) = r;
end
roiOut = out(keep);
end
