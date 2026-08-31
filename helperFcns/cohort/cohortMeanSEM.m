function out = cohortMeanSEM(M)
% cohortMeanSEM  Between-cell mean and SEM of a cells x frames matrix.
%
%   out = cohortMeanSEM(M)
%
%   The single place cohort-level mean/SEM is computed, so every plotter
%   inherits the same small-n policy.
%
%   Input:
%     M - nCells x nFrames numeric matrix. One row per CELL (trial-averaged
%         first). Must be 2-D; a cells x frames matrix with nCells == 1 is
%         valid and handled explicitly.
%
%   Output (struct):
%     .mean       1 x nFrames, mean across cells (NaN-omitting)
%     .sem        1 x nFrames, SEM across cells; ALL NaN when n < 2
%     .n          number of cells (rows of M)
%     .nPerFrame  1 x nFrames count of non-NaN cells contributing per frame
%     .showBand   true when n >= 2, i.e. when the SEM band is meaningful
%     .empty      true when n == 0
%
%   WHY THIS EXISTS (two small-n traps this closes):
%
%   1. SEMcalc(x) called WITHOUT an explicit dimension transposes a row
%      vector -- `if isvector(data) && size(data,2)>1, data = data(:); end` --
%      so a 1 x nFrames single-cell trace returns a 1 x 1 SCALAR "SEM across
%      time" instead of a per-frame SEM. Silent and wrong. This function
%      always reduces along dim 1 and never reshapes.
%
%   2. SEM across a single cell is 0, not NaN, so a naive plot draws a
%      zero-width band that reads as "no variability" rather than "n = 1".
%      Here n < 2 yields NaN and showBand = false, so callers omit the band
%      and label the n instead.
%
%   See also groupN, gatherCellTraces, SEMcalc.

if ~isnumeric(M) || ~ismatrix(M)
    error('cohortMeanSEM:badInput', ...
        'M must be a 2-D numeric cells x frames matrix (got %s, %s).', ...
        class(M), mat2str(size(M)));
end

n = size(M,1);
out.n = n;
out.empty = (n == 0);

if n == 0
    out.mean = zeros(1,size(M,2));
    out.sem  = nan(1,size(M,2));
    out.nPerFrame = zeros(1,size(M,2));
    out.showBand = false;
    return
end

out.mean      = mean(M,1,'omitnan');
out.nPerFrame = sum(~isnan(M),1);

if n < 2
    out.sem = nan(1,size(M,2));
    out.showBand = false;
else
    out.sem = std(M,0,1,'omitnan') ./ sqrt(max(out.nPerFrame,1));
    out.sem(out.nPerFrame < 2) = NaN;   % frames with <2 contributing cells
    out.showBand = true;
end
end
