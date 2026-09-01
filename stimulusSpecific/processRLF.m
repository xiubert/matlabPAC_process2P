% processRLF  Rate-level function and population dF/F for one BPN group.
%
% Entry point for BPN cohort plots. Select a BPN_Group<g>.mat (built by
% aggregateStimGroup) and this produces the RLF, the population dF/F traces
% per sound level, and peak dF/F per level -- each panel stating how many
% cells and mice it summarises.
%
% Superseded the previous version's hardcoded list of per-animal
% *_anmlROI_BPNstimTable.mat paths: assembling animals into a group is now
% aggregateStimGroup's job (with a manifest and a provenance stamp), so this
% script only plots. To build or rebuild a group file:
%
%   manifest = struct('group','D','family','BPN', ...
%                     'animals',["TO0006","TO0007"], ...
%                     'cohortRoot','/media/DATA/Ophys/Jinbo/TOMT', ...
%                     'outDir','/media/DATA/Ophys/Jinbo/TOMT/aggregate data');
%   aggregateStimGroup(manifest);
%
% Works for any group size, including a single animal.
%
% See also plotBPNgroup, aggregateStimGroup, loadStimGroup, tableRLF, plotRLF.

addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'helperFcns')));

%% params (EDIT IF NEEDED)
nConsec    = 2;        % min consecutive sig dB levels for RLF inclusion
levels     = [];       % [] -> every level present in the group
traceCells = 'all';    % 'all' | 'included' (RLF-included cells only)

%% pick the group file
if ~exist('groupFile','var') || isempty(groupFile)
    [f,pth] = uigetfile({'BPN_Group*.mat','BPN group files'}, ...
        'Select a BPN group file');
    if isequal(f,0); error('processRLF:noFile','No group file selected.'); end
    groupFile = fullfile(pth,f);
end

%% plot
out = plotBPNgroup(groupFile, ...
    'nConsec',nConsec, 'levels',levels, 'traceCells',traceCells);

%% report
fprintf('\n%s\n', groupFile);
fprintf('  %s\n', out.N.label);
disp(out.N.perAnimal);
fprintf('  RLF: %d of %d cells included (>= %d consecutive significant levels)\n', ...
    out.rlf.nIncluded, out.rlf.nTotal, out.rlf.nConsec);
if out.N.singleAnimal
    fprintf('  NOTE: single-animal group -- no across-animal inference.\n');
end
