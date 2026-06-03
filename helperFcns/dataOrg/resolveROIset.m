function mc = resolveROIset(groupStim, roiSets, roiTifs, roiCounts)
% resolveROIset  Pick the moCorROI set that owns a stim group's tifs.
%
%   mc = resolveROIset(groupStim, roiSets, roiTifs, roiCounts)
%
%   A stim group's traces were extracted in a specific condition's ROI order,
%   so the moCorROI used to label them MUST be that same condition's set. The
%   primary key is tif MEMBERSHIP: the condition whose tif list contains all of
%   this group's tifs. Resolving by ROI count alone is unsafe -- two conditions
%   can share a count (e.g. every matched 256x256 cell happens to lie inside
%   the 256x128 crop), which would silently mislabel the group's traces with
%   the wrong condition's ROI ID order.
%
%   Inputs
%     groupStim  struct array of the stim group's tifFileList.stim entries
%                (each has .name and .moCorRawFroi).
%     roiSets    cell array, one moCorROI struct array per condition.
%     roiTifs    cell array, original tif names each ROI set covers ({} if a
%                legacy file lacked moCorTifNames).
%     roiCounts  numeric, numel(roiSets{f}) per condition.
%
%   Resolution
%     1. MEMBERSHIP: exactly one condition whose tifs contain the whole group.
%     2. FALLBACK (group spans conditions, or legacy files without tif lists):
%        match by trace row-count -- safe because conditions sharing a count
%        via intersectROIfiles also share ID order.
%
%   See also stimParam2ROI, anmlROIbyStimTable

gnames = {groupStim.name};

% 1. membership
own = false(numel(roiSets),1);
for f = 1:numel(roiSets)
    own(f) = ~isempty(roiTifs{f}) && all(ismember(gnames, roiTifs{f}));
end
if nnz(own) == 1
    mc = roiSets{find(own,1)};
    return
end

% 2. fallback: trace row-count
n = size(groupStim(1).moCorRawFroi,1);
hit = find(roiCounts==n,1,'first');
if isempty(hit)
    error('resolveROIset:roiCountMismatch',...
        'No moCorrROI set with %d ROIs to match this stim group''s traces.',n);
end
mc = roiSets{hit};
end
