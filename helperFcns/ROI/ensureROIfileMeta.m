function fixed = ensureROIfileMeta(roiDir,animal,tifList,tifFiles)
% ensureROIfileMeta  Backfill the bookkeeping FISSA needs into ROI files.
%
%   fixed = ensureROIfileMeta(roiDir,animal,tifList,tifFiles)
%
%   A <animal>_moCorrROI_<cond>.mat must carry three variables besides the ROIs
%   themselves: nTifs, tifIDXinAllTifList and moCorTifNames. The last is what
%   FISSAviaMatlab_prePostTreatment.py uses to build each group's ordered image
%   list, and it hard-errors without it -- ROI files drawn before that driver
%   existed have only the ROIs, so an otherwise valid animal fails at the
%   Python step with a message about re-saving section 5.
%
%   Everything missing is recomputed from the condition-split legend, which is
%   the same source section 5 used, so this is a backfill rather than a guess.
%   Files that already carry the variables are left untouched.
%
%   Inputs
%     roiDir    folder holding <animal>_moCorrROI_<cond>.mat -- the animal
%               folder in the flat layout, analysis/<run>/ otherwise. Not
%               necessarily the animal folder, hence not called dataPath.
%     animal    animal ID.
%     tifList   condition -> tif struct map (the §2 legend).
%     tifFiles  the full tif list, for tifIDXinAllTifList.
%
%   Output
%     fixed     cellstr of the conditions whose files were rewritten.
%
%   See also intersectROIfiles, remapROIfile, TIFcatROIgui

fixed = {};
conds = fieldnames(tifList);
for c = 1:numel(conds)
    cond = conds{c};
    f = fullfile(roiDir,[animal '_moCorrROI_' cond '.mat']);
    if ~isfile(f); continue, end

    S = load(f);
    need = {};
    if ~isfield(S,'nTifs') || isempty(S.nTifs)
        S.nTifs = numel(tifList.(cond));
        need{end+1} = 'nTifs'; %#ok<AGROW>
    end
    if ~isfield(S,'tifIDXinAllTifList') || isempty(S.tifIDXinAllTifList)
        S.tifIDXinAllTifList = ismember({tifFiles.name}',{tifList.(cond).name}');
        need{end+1} = 'tifIDXinAllTifList'; %#ok<AGROW>
    end
    if ~isfield(S,'moCorTifNames') || isempty(S.moCorTifNames)
        S.moCorTifNames = strrep({tifList.(cond).name}','.tif','_NoRMCorre.tif');
        need{end+1} = 'moCorTifNames'; %#ok<AGROW>
    end
    if isempty(need); continue, end

    if numel(S.moCorTifNames) ~= numel(tifList.(cond))
        warning('ensureROIfileMeta:countMismatch',...
            ['%s lists %d moCorTifNames but condition ''%s'' has %d tifs; '...
             'leaving the file alone.'],f,numel(S.moCorTifNames),cond,...
            numel(tifList.(cond)));
        continue
    end

    save(f,'-struct','S');
    fixed{end+1} = cond; %#ok<AGROW>
    fprintf('  %s: added %s\n',[animal '_moCorrROI_' cond '.mat'],strjoin(need,', '));
end
end
