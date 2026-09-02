function out = cellposeROIset(dataPath,animal,condMeans,tifList,tifFiles,opts)
% cellposeROIset  Segment once, write one _moCorrROI_<cond>.mat per condition.
%
%   out = cellposeROIset(dataPath,animal,condMeans,tifList,tifFiles)
%   out = cellposeROIset(...,'dilatePx',1,'edgeMarginPx',3)
%
%   The headless replacement for processAnimal2P §4-5 (TIFcatROIgui + save).
%   Rather than segmenting each condition separately -- which would give every
%   condition its own arbitrary label numbering and destroy pre/post cell
%   tracking -- this registers the condition mean images to each other,
%   segments the aligned average ONCE, and places that single ROI set into
%   each condition's frame. IDs then refer to the same cell in every
%   condition by construction, and intersectROIfiles (§6) becomes an
%   invariant check rather than the thing doing the matching.
%
%   Inputs
%     dataPath   animal data folder; the ROI files are written here.
%     animal     animal ID, for the filenames.
%     condMeans  struct, one field per condition (names matching tifList),
%                each an H x W mean projection of that condition's
%                motion-corrected stack. All must be the same size: a 256x128
%                crop condition does NOT belong here, it is derived from a
%                256x256 condition afterwards via remapROIfile.
%     tifList    the condition -> tif struct map from §2.
%     tifFiles   the full tif list, for tifIDXinAllTifList.
%
%   Name-value -- ROI construction (forwarded to labelImg2moCorROI)
%     'dilatePx'      grow masks; CALIBRATE, do not guess. Default 0.
%     'minAreaPx'     Default 20.  'maxAreaPx'  Default Inf.
%     'edgeMarginPx'  keep ROIs this far inside the frame so FISSA's neuropil
%                     annulus is not clipped. Default 3.
%   Name-value -- registration
%     'maxShiftPx'    Default 15 (NoRMCorre's max_shift). See
%                     registerConditionMeans.
%     'refCond'       condition to register to. Default: the one with the
%                     most tifs (the best-averaged image).
%   Name-value -- other
%     'cellposeArgs'  cell array of name-value pairs forwarded verbatim to
%                     cellposeSegment (backend, model, diameter, ...).
%     'saveMeanTif'   write the segmented session mean next to the ROI files
%                     as <animal>_cellposeMean.tif, for QC. Default true.
%     'verbose'       Default true.
%
%   Output (struct)
%     .conds .refCond .shifts .sessionMean .labelImg .moCorROI (reference
%     frame) .roiInfo .cellposeParams .files (per-condition paths)
%     .perCond (name, nROI, nDropped)
%
%   Each written bundle carries the variables the rest of the pipeline reads
%   -- moCorROI, moCorSeqN, nTifs, tifIDXinAllTifList, moCorTifNames -- plus
%   cellposeParams and roiParams for provenance. moCorTifNames is not
%   optional: the FISSA driver errors without it.
%
%   See also cellposeSegment, labelImg2moCorROI, registerConditionMeans,
%   remapROIfile, intersectROIfiles, processAnimal2Pheadless

arguments
    dataPath   (1,:) char
    animal     (1,:) char
    condMeans  (1,1) struct
    tifList    (1,1) struct
    tifFiles         struct
    opts.dilatePx     (1,1) double {mustBeNonnegative,mustBeInteger} = 0
    opts.minAreaPx    (1,1) double {mustBeNonnegative} = 20
    opts.maxAreaPx    (1,1) double {mustBePositive}    = Inf
    opts.edgeMarginPx (1,1) double {mustBeNonnegative,mustBeInteger} = 3
    opts.maxShiftPx   (1,1) double {mustBeNonnegative} = 15
    opts.refCond      (1,:) char = ''
    opts.cellposeArgs       cell = {}
    opts.saveMeanTif  (1,1) logical = true
    opts.verbose      (1,1) logical = true
end

if ~isfolder(dataPath)
    error('cellposeROIset:noDataPath','Not a folder: %s',dataPath);
end

conds = fieldnames(condMeans);
if isempty(conds)
    error('cellposeROIset:noConditions','condMeans has no conditions.');
end
missing = conds(~isfield(tifList,conds));
if ~isempty(missing)
    error('cellposeROIset:condNotInTifList',...
        'Condition(s) %s are in condMeans but not in tifList.',strjoin(missing',', '));
end

%--- reference condition: the best-averaged image -------------------------
allCond = fieldnames(tifList);
if isempty(opts.refCond)
    nT = cellfun(@(c) numel(tifList.(c)),conds);
    [~,refIdx] = max(nT);
else
    refIdx = find(strcmp(conds,opts.refCond),1);
    if isempty(refIdx)
        error('cellposeROIset:badRefCond',...
            'refCond ''%s'' is not one of: %s',opts.refCond,strjoin(conds',', '));
    end
end

%--- register + average ---------------------------------------------------
meanCell = cellfun(@(c) double(condMeans.(c)),conds,'uni',0);
if isscalar(meanCell)
    shifts      = [0 0];
    sessionMean = meanCell{1};
    regInfo     = struct('refIdx',1,'shiftMagPx',0,'method','single condition',...
        'names',{conds'},'maxShiftPx',opts.maxShiftPx);
    if opts.verbose
        fprintf('cellposeROIset: one condition (''%s''), no registration needed.\n',conds{1});
    end
else
    [shifts,sessionMean,regInfo] = registerConditionMeans(meanCell,...
        'maxShiftPx',opts.maxShiftPx,'refIdx',refIdx,'names',conds',...
        'verbose',opts.verbose);
end

%--- segment the session mean --------------------------------------------
%NaN appears only in the thin border a condition could not contribute to
%after shifting; Cellpose cannot take NaN, and filling with the background
%median keeps that border from reading as an object.
segImg = sessionMean;
nanPix = ~isfinite(segImg);
if any(nanPix(:))
    segImg(nanPix) = median(segImg(~nanPix));
end

srcTif = fullfile(tifList.(conds{refIdx})(1).folder,tifList.(conds{refIdx})(1).name);
cpArgs = [{'name',sprintf('%s_cellposeMean',animal)},opts.cellposeArgs];
if ~any(strcmp(cpArgs(1:2:end),'srcTif')) && isfile(srcTif)
    cpArgs = [cpArgs,{'srcTif',srcTif}];
end
[labelImg,cellposeParams] = cellposeSegment(segImg,cpArgs{:});

if opts.verbose
    fprintf('cellposeROIset: Cellpose returned %d labels in %.1f s\n',...
        cellposeParams.nLabels,cellposeParams.runSeconds);
end

%--- reference-frame ROI set ---------------------------------------------
[roiRef,roiInfo] = labelImg2moCorROI(labelImg,...
    'dilatePx',opts.dilatePx,'minAreaPx',opts.minAreaPx,...
    'maxAreaPx',opts.maxAreaPx,'edgeMarginPx',opts.edgeMarginPx);
if isempty(roiRef)
    error('cellposeROIset:noROI',...
        ['Segmentation produced no usable ROIs (%d labels, all dropped by QC). '...
         'Check the session mean image and the area/edge settings.'],roiInfo.nLabels);
end
if opts.verbose
    fprintf('cellposeROIset: %d ROIs kept (dropped %d small, %d large, %d edge)\n',...
        roiInfo.nKept,roiInfo.dropped.minArea,roiInfo.dropped.maxArea,roiInfo.dropped.edge);
end

roiParams = struct('dilatePx',opts.dilatePx,'minAreaPx',opts.minAreaPx,...
    'maxAreaPx',opts.maxAreaPx,'edgeMarginPx',opts.edgeMarginPx,...
    'registration',regInfo,'shifts',shifts,'conds',{conds'},...
    'refCond',conds{refIdx},'source','cellposeROIset');

%--- optional QC artefact -------------------------------------------------
meanTifPath = '';
if opts.saveMeanTif
    meanTifPath = fullfile(dataPath,[animal '_cellposeMean.tif']);
    try
        if isfile(srcTif)
            writeTifWithHeader(segImg,meanTifPath,srcTif,...
                'class','int16','copyright','cellposeROIset_sessionMean');
        end
    catch ME
        warning('cellposeROIset:meanTifFailed',...
            'Could not write %s (%s); continuing.',meanTifPath,ME.message);
        meanTifPath = '';
    end
end

%--- place the set into every condition and save -------------------------
files   = cell(numel(conds),1);
perCond = struct('name',{},'nROI',{},'nDropped',{});
for k = 1:numel(conds)
    cond = conds{k};
    %shifts(k,:) maps condition k INTO the reference frame, so a
    %reference-frame ROI goes the other way
    roiK = shiftROIset(roiRef,-shifts(k,:),size(labelImg),opts.edgeMarginPx);
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
    S.cellposeLabelImg   = uint16(labelImg);

    files{k} = fullfile(dataPath,[animal '_moCorrROI_' cond '.mat']);
    save(files{k},'-struct','S');

    perCond(k).name     = cond;
    perCond(k).nROI     = numel(roiK);
    perCond(k).nDropped = numel(roiRef) - numel(roiK);
    if opts.verbose
        fprintf('  %-24s %3d ROIs -> %s\n',cond,numel(roiK),...
            [animal '_moCorrROI_' cond '.mat']);
    end
end

out = struct('conds',{conds'},'refCond',conds{refIdx},'shifts',shifts,...
    'sessionMean',sessionMean,'labelImg',labelImg,'moCorROI',roiRef,...
    'roiInfo',roiInfo,'roiParams',roiParams,'cellposeParams',cellposeParams,...
    'files',{files},'perCond',perCond,'meanTif',meanTifPath);
end

% ------------------------------------------------------------------------
function roiOut = shiftROIset(roiIn,dxdy,imgSize,edgeMarginPx)
%translate a whole ROI set by an integer [dx dy], dropping any ROI that would
%no longer sit fully inside the frame (with the same margin the QC used).
%Mirrors remapROItoAcq's regeneration of the derived geometry fields so FISSA
%receives a well-formed polygon either way.
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
    if ~any(m(:))
        continue
    end
    rws = find(any(m,2));
    cls = find(any(m,1));
    if rws(1) <= edgeMarginPx || rws(end) > H-edgeMarginPx || ...
       cls(1) <= edgeMarginPx || cls(end) > W-edgeMarginPx
        continue
    end
    %a mask clipped by the translation is not the same cell any more
    if nnz(m) ~= nnz(roiIn(k).mask)
        continue
    end
    keep(k) = true;

    r = out(k);
    r.mask                = m;
    r.XYvertices          = [r.XYvertices(:,1)+dx, r.XYvertices(:,2)+dy];
    r.pos(1)              = r.pos(1) + dx;
    r.pos(2)              = r.pos(2) + dy;
    r.label.Position(1)   = r.label.Position(1) + dx;
    r.label.Position(2)   = r.label.Position(2) + dy;
    r.ROIxyCoord          = mask2polyCoord(m);
    r.ROIcurveOrderedXY   = [r.ROIcurveOrderedXY(1,:)+dx; r.ROIcurveOrderedXY(2,:)+dy];
    out(k) = r;
end
roiOut = out(keep);
end
