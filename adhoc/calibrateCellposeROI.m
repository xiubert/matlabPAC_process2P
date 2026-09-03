function R = calibrateCellposeROI(dataPath,opts)
% calibrateCellposeROI  Calibrate Cellpose ROI settings against hand-drawn ROIs.
%
%   R = calibrateCellposeROI(dataPath)
%   R = calibrateCellposeROI(dataPath,'dilateSweep',0:3,'iouThreshold',0.3)
%
%   Segments an animal's motion-corrected session mean with Cellpose, sweeps
%   the mask dilation, and compares each result against that animal's
%   hand-drawn moCorROI set. Answers one question: what dilatePx makes
%   segmented masks the same SIZE as the hand-drawn convention, so dF/F
%   amplitudes stay comparable to the existing corpus.
%
%   What is scored and what is not
%     Scored     : median IoU and area ratio over MATCHED cells, and the
%                  correlation of the raw F traces those masks produce.
%     Not scored : cell COUNT. Hand labelling under-counts, so Cellpose ROIs
%                  with no hand-drawn partner are candidate missed cells, not
%                  false positives. They are reported, and their area
%                  distribution is compared against the matched ones -- if
%                  they look like the matched masks they are probably cells;
%                  if they are much smaller they are probably debris.
%
%   Inputs
%     dataPath  animal folder holding <animal>_moCorrROI_<cond>.mat,
%               <animal>_tifCondSplitLegend.mat and NoRMCorred/.
%
%   Name-value
%     'cond'          condition to calibrate on. Default: the one whose
%                     hand-drawn ROI file exists (errors if ambiguous).
%     'dilateSweep'   dilations to try. Default 0:3.
%     'iouThreshold'  Default 0.3.
%     'edgeMarginPx'  Default 3.  'minAreaPx' Default 20.
%     'traceTif'      index of the condition's tif used for the trace
%                     comparison. Default 1. Set 0 to skip.
%     'cellposeArgs'  forwarded to cellposeSegment.
%     'plot'          show mean image + both ROI sets. Default true.
%
%   Output (struct)
%     .animal .cond .sessionMean .labelImg .cellposeParams
%     .human (the hand-drawn set) .sweep (table: dilatePx, nCellpose,
%     nMatched, medianIoU, areaRatio, medEquivDiam) .best (chosen dilatePx)
%     .trace (r per matched cell, at the best dilation)
%
%   See also cellposeROIset, compareROIsets, labelImg2moCorROI, condMeanImg

arguments
    dataPath (1,:) char
    opts.cond          (1,:) char = ''
    opts.dilateSweep         double = 0:3
    opts.iouThreshold  (1,1) double = 0.3
    opts.edgeMarginPx  (1,1) double = 3
    opts.minAreaPx     (1,1) double = 20
    opts.traceTif      (1,1) double = 1
    opts.cellposeArgs        cell = {}
    opts.plot          (1,1) logical = true
end

if ~isfolder(dataPath)
    error('calibrateCellposeROI:noDataPath','Not a folder: %s',dataPath);
end
animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');
if isempty(animal)
    error('calibrateCellposeROI:noAnimal',...
        'Could not derive an animal ID from %s',dataPath);
end

%--- the hand-drawn set ---------------------------------------------------
roiFiles = dir(fullfile(dataPath,[animal '_moCorrROI_*.mat']));
if isempty(roiFiles)
    error('calibrateCellposeROI:noHumanROI',...
        'No %s_moCorrROI_*.mat in %s -- nothing to calibrate against.',animal,dataPath);
end
if isempty(opts.cond)
    if ~isscalar(roiFiles)
        error('calibrateCellposeROI:ambiguousCond',...
            'Several ROI files (%s); pass ''cond''.',...
            strjoin({roiFiles.name},', '));
    end
    cond = regexp(roiFiles(1).name,['(?<=' animal '_moCorrROI_).*(?=\.mat)'],'match','once');
else
    cond = opts.cond;
end
humanFile = fullfile(dataPath,[animal '_moCorrROI_' cond '.mat']);
Sh = load(humanFile,'moCorROI');
human = Sh.moCorROI;
human = human(~[human.deleted]);
fprintf('Hand-drawn set: %d ROIs from %s\n',numel(human),[animal '_moCorrROI_' cond '.mat']);

%--- the session mean -----------------------------------------------------
Sc = load(fullfile(dataPath,[animal '_tifCondSplitLegend.mat']),'tifList');
if ~isfield(Sc.tifList,cond)
    error('calibrateCellposeROI:condNotInList',...
        'Condition ''%s'' is not in the split legend (%s).',...
        cond,strjoin(fieldnames(Sc.tifList)',', '));
end
condTifs = Sc.tifList.(cond);
fprintf('Projecting %d motion-corrected tifs...\n',numel(condTifs));
sessionMean = condMeanImg(dataPath,condTifs);

%--- segment once ---------------------------------------------------------
srcTif = fullfile(dataPath,'NoRMCorred',...
    strrep(condTifs(1).name,'.tif','_NoRMCorre.tif'));
cpArgs = [{'name',sprintf('%s_calib',animal)},opts.cellposeArgs];
if ~any(strcmp(cpArgs(1:2:end),'srcTif')) && isfile(srcTif)
    cpArgs = [cpArgs,{'srcTif',srcTif}];
end
[labelImg,cpParams] = cellposeSegment(sessionMean,cpArgs{:});
fprintf('Cellpose: %d labels in %.1f s\n',cpParams.nLabels,cpParams.runSeconds);

%--- sweep ----------------------------------------------------------------
nSw = numel(opts.dilateSweep);
[nCP,nMatch,medIoU,aRatio,medDiam,medAreaOnlyB,medAreaMatched] = deal(nan(nSw,1));
roiBySweep = cell(nSw,1);
cmpBySweep = cell(nSw,1);

for s = 1:nSw
    d = opts.dilateSweep(s);
    [roi,info] = labelImg2moCorROI(labelImg,'dilatePx',d,...
        'minAreaPx',opts.minAreaPx,'edgeMarginPx',opts.edgeMarginPx);
    cmp = compareROIsets(human,roi,'iouThreshold',opts.iouThreshold);

    roiBySweep{s} = roi;
    cmpBySweep{s} = cmp;
    nCP(s)     = numel(roi);
    nMatch(s)  = cmp.nMatched;
    medIoU(s)  = cmp.medianIoU;
    aRatio(s)  = cmp.areaRatio;
    medDiam(s) = median(info.equivDiamPx);
    if ~isempty(cmp.areaOnlyB),  medAreaOnlyB(s)   = median(cmp.areaOnlyB); end
    if ~isempty(cmp.pairs),      medAreaMatched(s) = median(cmp.areaB(cmp.pairs(:,2))); end
end

sweep = table(opts.dilateSweep(:),nCP,nMatch,medIoU,aRatio,medDiam,...
    medAreaMatched,medAreaOnlyB,...
    'VariableNames',{'dilatePx','nCellpose','nMatched','medianIoU',...
    'areaRatio','medEquivDiamPx','medAreaMatched','medAreaOnlyB'});

fprintf('\n--- dilation sweep (%d hand-drawn ROIs, IoU>=%.2f) ---\n',...
    numel(human),opts.iouThreshold);
disp(sweep)

[~,bestIdx] = max(medIoU);
best = opts.dilateSweep(bestIdx);
fprintf(['Best median IoU at dilatePx = %d (IoU %.3f, %d/%d hand-drawn ROIs '...
    'matched, area ratio %.2f)\n'],best,medIoU(bestIdx),nMatch(bestIdx),...
    numel(human),aRatio(bestIdx));

cmpBest = cmpBySweep{bestIdx};
roiBest = roiBySweep{bestIdx};
if ~isempty(cmpBest.onlyB)
    fprintf(['Cellpose-only ROIs: %d. Median area %.0f px vs %.0f px for '...
        'matched -- %s\n'],numel(cmpBest.onlyB),median(cmpBest.areaOnlyB),...
        median(cmpBest.areaB(cmpBest.pairs(:,2))),...
        sizeVerdict(median(cmpBest.areaOnlyB),median(cmpBest.areaB(cmpBest.pairs(:,2)))));
end

%--- do the masks measure the same cells? --------------------------------
trace = struct('r',[],'relDiff',[],'tif','');
if opts.traceTif > 0 && ~isempty(cmpBest.pairs)
    tifIdx = min(opts.traceTif,numel(condTifs));
    tf = fullfile(dataPath,'NoRMCorred',...
        strrep(condTifs(tifIdx).name,'.tif','_NoRMCorre.tif'));
    fprintf('\nTrace comparison on %s...\n',condTifs(tifIdx).name);
    stack = double(justLoadTif(tf));
    nF = size(stack,3);
    P  = size(cmpBest.pairs,1);
    [rr,dd] = deal(nan(P,1));
    for p = 1:P
        mH = human(cmpBest.pairs(p,1)).mask;
        mC = roiBest(cmpBest.pairs(p,2)).mask;
        [fH,fC] = deal(zeros(1,nF));
        for k = 1:nF
            fr = stack(:,:,k);
            fH(k) = mean(fr(mH));
            fC(k) = mean(fr(mC));
        end
        rr(p) = corr(fH(:),fC(:));
        dd(p) = median(abs(fC-fH))/median(fH);
    end
    trace.r       = rr;
    trace.relDiff = dd;
    trace.tif     = condTifs(tifIdx).name;
    fprintf('Matched-cell raw F: median r = %.3f (min %.3f), median |dF|/F = %.3f\n',...
        median(rr),min(rr),median(dd));
end

%--- picture --------------------------------------------------------------
if opts.plot
    figure('Name',sprintf('%s Cellpose vs hand-drawn (dilatePx=%d)',animal,best));
    imagesc(sessionMean); axis image; colormap(gray); hold on
    for k = 1:numel(human)
        c = human(k).ROIcurveOrderedXY;
        plot(c(1,:),c(2,:),'c-','LineWidth',1.2);
    end
    for k = 1:numel(roiBest)
        c = roiBest(k).ROIcurveOrderedXY;
        isMatched = ismember(k,cmpBest.pairs(:,2));
        if isMatched
            plot(c(1,:),c(2,:),'y-','LineWidth',1);
        else
            plot(c(1,:),c(2,:),'r--','LineWidth',1);
        end
    end
    title(sprintf(['%s %s | cyan = hand-drawn (%d), yellow = matched (%d), '...
        'red = Cellpose-only (%d)'],animal,cond,numel(human),...
        cmpBest.nMatched,numel(cmpBest.onlyB)));
end

R = struct('animal',animal,'cond',cond,'sessionMean',sessionMean,...
    'labelImg',labelImg,'cellposeParams',cpParams,'human',human,...
    'sweep',sweep,'best',best,'roiBest',roiBest,'cmpBest',cmpBest,...
    'trace',trace,'iouThreshold',opts.iouThreshold);
end

% ------------------------------------------------------------------------
function s = sizeVerdict(aOnly,aMatched)
r = aOnly/aMatched;
if r > 0.7
    s = 'comparable, so probably cells the hand-drawing missed';
elseif r > 0.4
    s = 'somewhat smaller; inspect them before trusting the extra ROIs';
else
    s = 'much smaller, so more likely debris -- consider raising minAreaPx';
end
end
