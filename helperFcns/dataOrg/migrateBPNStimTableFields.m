function migrateBPNStimTableFields(matFile)
% migrateBPNStimTableFields  Rename legacy BPN stim-table variables in a
%                            saved *_anmlROI_BPNstimTable.mat file.
%
%   migrateBPNStimTableFields(matFile)
%
%   Loads matFile, renames any legacy variable names on the anmlROIbyStim
%   and stimTable tables to the new <StimFamily><unit><Quantity> scheme:
%       sonset        -> BPNsOnset
%       dBampl        -> BPNdBAmpl
%       mspulseLen    -> BPNmsPulseLen
%       msStimLen     -> BPNmsStimLen
%       BPNfreqcutoff -> BPNkHzFreqCutoff
%
%   Safe to re-run: variables already on the new name are skipped. The
%   .mat file is only resaved if at least one rename occurred.
%
%   Inputs:
%     matFile - absolute path to a *_anmlROI_BPNstimTable.mat file.

oldNames = {'sonset','dBampl','mspulseLen','msStimLen','BPNfreqcutoff'};
newNames = {'BPNsOnset','BPNdBAmpl','BPNmsPulseLen','BPNmsStimLen','BPNkHzFreqCutoff'};

S = load(matFile);
varList = fieldnames(S);

changed = false;
for v = 1:numel(varList)
    val = S.(varList{v});
    if ~istable(val)
        continue
    end
    present = ismember(oldNames, val.Properties.VariableNames);
    if ~any(present)
        continue
    end
    val = renamevars(val, oldNames(present), newNames(present));
    S.(varList{v}) = val;
    changed = true;
    fprintf('migrateBPNStimTableFields: renamed %d variable(s) on %s in %s\n', ...
        sum(present), varList{v}, matFile);
end

if changed
    save(matFile, '-struct', 'S', '-v7.3');
    fprintf('migrateBPNStimTableFields: resaved %s\n', matFile);
else
    fprintf('migrateBPNStimTableFields: no legacy names found in %s; nothing to do.\n', matFile);
end
end
