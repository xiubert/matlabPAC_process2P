function [shifts,sessionMean,info] = registerConditionMeans(condMeans,opts)
% registerConditionMeans  Align per-condition mean images so one ROI set fits all.
%
%   [shifts,sessionMean,info] = registerConditionMeans(condMeans)
%   [...] = registerConditionMeans(condMeans,'maxShiftPx',15,'refIdx',1)
%
%   Each motion-correction condition is registered to its OWN template by
%   NoRMCorre, so pre- and post-treatment stacks of the same field can sit a
%   few pixels apart. A human redraws ROIs per condition and keeps the same
%   cell under the same ID by eye; a segmenter cannot. This measures the
%   residual translation between condition mean images so the headless path
%   can segment ONCE and place that single ROI set into every condition --
%   which is what makes IDs mean the same cell across conditions.
%
%   Translation only, and rounded to whole pixels. The acquisitions share a
%   zoom and a field of view (processAnimal2P §2 already split conditions by
%   geometry), so the residual really is a small shift; and ROI masks are
%   logical, so a sub-pixel estimate would have to be interpolated back into
%   a binary mask, which costs more than the <1 px it would buy.
%
%   Inputs
%     condMeans  cell array of H x W mean images, one per condition, all the
%                same size (they must be -- conditions of differing geometry
%                are separate groups and 256x128 crops go through
%                remapROItoAcq instead).
%
%   Name-value
%     'maxShiftPx'  refuse a shift larger than this, in pixels. Default 15 --
%                   NoRMCorre's own max_shift in this pipeline. A larger
%                   offset means the two acquisitions are not looking at the
%                   same patch of tissue, and forcing one shared ROI set onto
%                   them would sample the wrong cells. Erroring is the point.
%     'refIdx'      condition all others are registered to. Default 1.
%     'names'       condition names, for readable messages. Default {}.
%     'verbose'     print the estimated shifts. Default true.
%
%   Outputs
%     shifts       nCond x 2 integer [dx dy] mapping each condition's frame
%                  INTO the reference frame. Row refIdx is [0 0].
%                  To place a reference-frame ROI into condition k, translate
%                  it by -shifts(k,:).
%     sessionMean  H x W mean of every condition mean after alignment (NaN
%                  where a condition contributed no pixel). This is the image
%                  to segment: a cell active in only one condition still
%                  shows up, instead of being diluted by an unregistered
%                  average.
%     info         .refIdx .shiftMagPx .method .names
%
%   See also cellposeROIset, remapROItoAcq, intersectROIfiles

arguments
    condMeans          cell
    opts.maxShiftPx (1,1) double {mustBeNonnegative} = 15
    opts.refIdx     (1,1) double {mustBePositive,mustBeInteger} = 1
    opts.names            cell = {}
    opts.verbose    (1,1) logical = true
end

nCond = numel(condMeans);
if nCond == 0
    error('registerConditionMeans:noConditions','condMeans is empty.');
end
if opts.refIdx > nCond
    error('registerConditionMeans:badRef',...
        'refIdx %d exceeds the %d conditions supplied.',opts.refIdx,nCond);
end

names = opts.names;
if isempty(names)
    names = arrayfun(@(k) sprintf('cond%d',k),1:nCond,'uni',0);
elseif numel(names) ~= nCond
    error('registerConditionMeans:badNames',...
        '%d names for %d conditions.',numel(names),nCond);
end

sz = size(condMeans{1});
for k = 2:nCond
    if ~isequal(size(condMeans{k}),sz)
        error('registerConditionMeans:sizeMismatch',...
            ['Condition ''%s'' is %s but ''%s'' is %s. Conditions of different '...
             'geometry cannot share one ROI set; they are separate groups, and '...
             'a 256x128 crop is handled by remapROItoAcq.'],...
            names{k},mat2str(size(condMeans{k})),names{1},mat2str(sz));
    end
end

fixed  = prepImage(condMeans{opts.refIdx});
shifts = zeros(nCond,2);
method = 'imregcorr(translation)';

for k = 1:nCond
    if k == opts.refIdx, continue, end
    moving = prepImage(condMeans{k});
    t = estimateShift(moving,fixed);
    shifts(k,:) = round(t);
end

shiftMag = hypot(shifts(:,1),shifts(:,2));
bad = find(shiftMag > opts.maxShiftPx);
if ~isempty(bad)
    msg = strjoin(arrayfun(@(k) sprintf('%s: [%d %d] (%.1f px)',...
        names{k},shifts(k,1),shifts(k,2),shiftMag(k)),bad(:)','uni',0),'; ');
    error('registerConditionMeans:shiftTooLarge',...
        ['Condition mean images are %s from the reference (''%s''), above '...
         'maxShiftPx = %g.\nA shift this large means the acquisitions are not '...
         'the same field of view, so one shared ROI set would sample different '...
         'cells in each condition. Draw ROIs per condition, or raise '...
         'maxShiftPx only if you have checked the images overlap.'],...
        msg,names{opts.refIdx},opts.maxShiftPx);
end

%--- aligned average ------------------------------------------------------
acc = zeros(sz);
cnt = zeros(sz);
for k = 1:nCond
    a = imtranslate(double(condMeans{k}),shifts(k,:),'FillValues',NaN);
    valid = ~isnan(a);
    acc(valid) = acc(valid) + a(valid);
    cnt = cnt + valid;
end
sessionMean = acc./cnt;
sessionMean(cnt==0) = NaN;

if opts.verbose
    fprintf('registerConditionMeans: reference ''%s''\n',names{opts.refIdx});
    for k = 1:nCond
        fprintf('  %-24s shift [%+d %+d] px (%.1f)\n',...
            names{k},shifts(k,1),shifts(k,2),shiftMag(k));
    end
end

info = struct('refIdx',opts.refIdx,'shiftMagPx',shiftMag(:)',...
    'method',method,'names',{names},'maxShiftPx',opts.maxShiftPx);
end

% ------------------------------------------------------------------------
function a = prepImage(a)
%registration wants a finite, comparably-scaled image
a = double(a);
a(~isfinite(a)) = 0;
r = max(a(:)) - min(a(:));
if r > 0
    a = (a - min(a(:)))./r;
end
end
% ------------------------------------------------------------------------
function t = estimateShift(moving,fixed)
%[tx ty] taking moving into fixed. imregcorr returns a geometric transform
%whose class changed across releases (affine2d -> transltform2d), so read the
%translation out of whichever we get.
tform = imregcorr(moving,fixed,'translation');
if isprop(tform,'Translation')            % transltform2d (R2022b+)
    t = double(tform.Translation);
elseif isprop(tform,'T')                  % affine2d
    t = double(tform.T(3,1:2));
else
    error('registerConditionMeans:unknownTform',...
        'Cannot read a translation out of a %s returned by imregcorr.',class(tform));
end
t = t(:)';
end
