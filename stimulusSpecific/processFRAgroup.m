% processFRAgroup  Cohort tuning plots for one FRA condition group.
%
% Entry point for FRA cohort figures, mirroring processRLF for BPN. Select an
% FRA_Group<g>.mat (built by aggregateStimGroup) and this produces the four
% population panels: best-frequency distribution, threshold distribution,
% bandwidth, and frequency response area colormaps.
%
% To build or rebuild a group file, resolving each animal's table wherever it
% lives (TO0004 and TO0005 keep theirs under Region1/):
%
%   root = '/media/DATA/Ophys/Jinbo/TOMT/animals';
%   files = arrayfun(@(a) string(fullfile(root,a,a+"_anmlROI_FRAtable.mat")), ...
%                    ["TO0009","TO0010","TO0011"]);
%   manifest = struct('group','B','family','FRA', ...
%                     'animals',["TO0009","TO0010","TO0011"], ...
%                     'cohortRoot',root, ...
%                     'outDir','/media/DATA/Ophys/Jinbo/TOMT/aggregate data', ...
%                     'files',{cellstr(files)});
%   aggregateStimGroup(manifest);
%
% The per-animal table is written by processFRA alongside <animal>_FRAmap.mat
% and holds one row per (ROI, frequency, level).
%
% TWO THINGS TO KNOW BEFORE READING THE FIGURES:
%
%   1. Bandwidth is BW20, not BW10. Levels are 20 dB apart (30/50/70), so
%      threshold + 10 dB is never a level that was presented.
%   2. Every run prints a noise floor: the significance rate in the real
%      response window over the rate the identical test gives on a silent
%      late window. A ratio near 1 means the mask underlying the threshold
%      and bandwidth panels carries little stimulus information.
%
% See also plotFRAgroup, tableFRAmetrics, aggregateStimGroup, processRLF.

addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'helperFcns')));

%% params (EDIT IF NEEDED)
nConsec  = 1;      % consecutive significant levels defining threshold
minBand  = 2;      % contiguous significant frequencies for a level to count
nExample = 8;      % example single-cell FRAs in the colormap figure
sigOnly  = true;   % colormaps from significant peaks only

%% pick the group file
if ~exist('groupFile','var') || isempty(groupFile)
    [f,pth] = uigetfile({'FRA_Group*.mat','FRA group files'}, ...
        'Select an FRA group file');
    if isequal(f,0); error('processFRAgroup:noFile','No group file selected.'); end
    groupFile = fullfile(pth,f);
end

%% plot
out = plotFRAgroup(groupFile, ...
    'nConsec',nConsec,'minBand',minBand,'nExample',nExample,'sigOnly',sigOnly);
