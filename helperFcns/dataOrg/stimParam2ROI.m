function [pulseLegend2P,stimGroupIDX,outputTables] = stimParam2ROI(dataPath)
% stimParam2ROI  Organize one animal's 2P traces by ROI x stim parameters.
%
%   [pulseLegend2P,stimGroupIDX,outputTables] = stimParam2ROI(dataPath)
%
%   For an animal data folder that already contains moCorROI traces,
%   tifFileList, and *_Pulses.mat files, this function:
%     1. Builds pulseLegend2P from the *_Pulses.mat files via
%        tifPulseLegend2P.
%     2. Classifies each tif into a stim group based on the pulseSet
%        name, producing stimGroupIDX.
%     3. For each stim group with more than one tif, builds a
%        tifStimParamTable, an anmlROIbyStim table, and a deduplicated
%        stimTable via stimParams2TifTable + anmlROIbyStimTable, and
%        saves them to a stim-specific .mat file in dataPath.
%
%   Layout of resulting per-stim tables: one row per (ROI, unique-stim)
%   pair, so all repetitions of a given stim on a given ROI share a row.
%
%   Inputs:
%     dataPath - absolute path to an animal data folder. The animal ID
%                (pattern '[A-Z]{2}\d{4}', e.g. AA0067) is parsed from the
%                path and used to locate the required input files:
%                  <animal>_moCorrROI_*.mat (loads moCorROI)
%                  <animal>_tifFileList.mat (loads tifFileList)
%                  *_Pulses.mat (consumed by tifPulseLegend2P)
%
%   Outputs:
%     pulseLegend2P - struct array from tifPulseLegend2P (one entry per
%                     *_Pulses.mat file). Also saved to
%                     <animal>_pulseLegend2P.mat.
%     stimGroupIDX  - struct with one field per stim group, each holding:
%                       .pulseLegend2P - logical mask into pulseLegend2P
%                       .tifFileList   - logical mask into tifFileList.stim
%                     Detected groups (matched against pulse.pulseSet):
%                       ptStimIDX         - 'PTinContrast'
%                       BPNStimIDX        - 'BPN'
%                       spontStimIDX      - 'spont'
%                       contrastChangeIDX - 'contrastChange'
%                       pupilReflexIDX    - 'LED trigger'
%                     Also saved to <animal>_stimGroupIDX.mat.
%     outputTables  - flat cell array of interleaved name/value pairs
%                     {'tifStimParamTable', T1, 'anmlROIbyStim', T2, ...}
%                     spanning every stim group that triggered processing.
%                     Empty if no group has >1 tif.
%
%   Side-effect .mat files written (only for groups with >1 tif):
%     <animal>_pulseLegend2P.mat
%     <animal>_stimGroupIDX.mat
%     <animal>_anmlROI_CGCstimTable.mat    (PTinContrast)
%     <animal>_anmlROI_BPNstimTable.mat    (BPN)
%     <animal>_anmlROI_SpontstimTable.mat  (spont)
%     <animal>_anmlROI_dContrastTable.mat  (contrastChange)
%
%   Notes:
%     - Pupil-reflex tifs are detected and indexed but not currently
%       processed into a per-stim table.
%     - The Spont branch has previously surfaced a tabular/unique error
%       when a non-char cell column made it into the dedup step; the
%       captured error trace is preserved in-source above the spont
%       block as a regression hint.

outputTables = cell(0);
animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');

tmp0 = cellstr(ls(fullfile(dataPath,[animal '_moCorrROI_*.mat'])));
tmp = load(fullfile(dataPath,tmp0{1}),...
    'moCorROI');
moCorROI = tmp.moCorROI;
clear tmp*

load(fullfile(dataPath,[animal '_tifFileList.mat']),'tifFileList')

% get pulse name for each tif
pulseLegend2P = tifPulseLegend2P(dataPath);

save(fullfile(dataPath,[animal '_pulseLegend2P.mat']),'pulseLegend2P','-v7.3')

% sort by stimulus group
% for PT in contrast
stimGroupIDX.ptStimIDX.pulseLegend2P = contains({pulseLegend2P.pulseSet},'PTinContrast')';
stimGroupIDX.ptStimIDX.tifFileList = ismember({tifFileList.stim.name},...
    {pulseLegend2P(stimGroupIDX.ptStimIDX.pulseLegend2P).tif})';
% for BPN
stimGroupIDX.BPNStimIDX.pulseLegend2P = contains({pulseLegend2P.pulseSet},'BPN')';
stimGroupIDX.BPNStimIDX.tifFileList = ismember({tifFileList.stim.name},...
    {pulseLegend2P(stimGroupIDX.BPNStimIDX.pulseLegend2P).tif})';
% for spont
stimGroupIDX.spontStimIDX.pulseLegend2P = contains({pulseLegend2P.pulseSet},'spont')';
stimGroupIDX.spontStimIDX.tifFileList = ismember({tifFileList.stim.name},...
    {pulseLegend2P(stimGroupIDX.spontStimIDX.pulseLegend2P).tif})';
% for contrast change
stimGroupIDX.contrastChangeIDX.pulseLegend2P = contains({pulseLegend2P.pulseSet},'contrastChange')';
stimGroupIDX.contrastChangeIDX.tifFileList = ismember({tifFileList.stim.name},...
    {pulseLegend2P(stimGroupIDX.contrastChangeIDX.pulseLegend2P).tif})';
% for pupil reflex
stimGroupIDX.pupilReflexIDX.pulseLegend2P = contains({pulseLegend2P.pulseSet},'LED trigger')';
stimGroupIDX.pupilReflexIDX.tifFileList = ismember({tifFileList.stim.name},...
    {pulseLegend2P(stimGroupIDX.pupilReflexIDX.pulseLegend2P).tif})';

save(fullfile(dataPath,[animal '_stimGroupIDX.mat']),'stimGroupIDX','-v7.3')

%sort pure tone in contrast stims
if sum(stimGroupIDX.ptStimIDX.tifFileList)>1

    % Build a per-tif stim-parameter table for the PT-in-contrast group:
    % one row per tif, scalar columns for per-tif params and cell columns
    % for per-pulse params. Feeds anmlROIbyStimTable's pivot/dedup step.
    tifStimParamTable = stimParams2TifTable(...
        tifFileList.stim(stimGroupIDX.ptStimIDX.tifFileList),dataPath);
    
    % sort ROIs by stim params
    [anmlROIbyStim,stimTable] = anmlROIbyStimTable(animal,...
        tifFileList.stim(stimGroupIDX.ptStimIDX.tifFileList),...
        moCorROI,...
        tifStimParamTable);
    outputTables{end+1} = 'tifStimParamTable';
    outputTables{end+1} = tifStimParamTable;
    outputTables{end+1} = 'anmlROIbyStim';
    outputTables{end+1} = anmlROIbyStim;
    outputTables{end+1} = 'stimTable';
    outputTables{end+1} = stimTable;
    
    save(fullfile(dataPath,[animal '_anmlROI_CGCstimTable.mat']),...
        'anmlROIbyStim','stimTable','tifStimParamTable',...
        'dataPath','-v7.3')
end

%BPN
if sum(stimGroupIDX.BPNStimIDX.tifFileList)>1

    % Build a per-tif stim-parameter table for the BPN group: one row per
    % tif, scalar columns for per-tif params (trigDelay, ISI, totalPulses)
    % and cell columns for per-pulse params (BPNsOnset, BPNdBAmpl, ...).
    % Feeds anmlROIbyStimTable's multi-pulse expansion + dedup.
    tifStimParamTable = stimParams2TifTable(...
        tifFileList.stim(stimGroupIDX.BPNStimIDX.tifFileList),dataPath);
    
    % sort ROIs by stim params
    [anmlROIbyStim,stimTable] = anmlROIbyStimTable(animal,...
        tifFileList.stim(stimGroupIDX.BPNStimIDX.tifFileList),...
        moCorROI,...
        tifStimParamTable);
    outputTables{end+1} = 'tifStimParamTable';
    outputTables{end+1} = tifStimParamTable;
    outputTables{end+1} = 'anmlROIbyStim';
    outputTables{end+1} = anmlROIbyStim;
    outputTables{end+1} = 'stimTable';
    outputTables{end+1} = stimTable;
    
    save(fullfile(dataPath,[animal '_anmlROI_BPNstimTable.mat']),...
        'anmlROIbyStim','stimTable','tifStimParamTable',...
        'dataPath','-v7.3')
end

% SPONT
%%%
% Error using tabular/unique (line 39)
% Unable to group rows using unique values of the table variable 'SpontMsStimLen'.  UNIQUE returned an error.
% 
% Error in anmlROIbyStimTable (line 86)
% [stimTable,~,ic] = unique(tmpT);
% 
% Error in stimParam2ROI (line 101)
%     [anmlROIbyStim,stimTable] = anmlROIbyStimTable(animal,...
% 
% Error in processAnimal2P (line 347)
% [pulseLegend2P,stimGroupIDX,ROIoutputTables] = stimParam2ROI(dataPath);
% 
% Caused by:
%     Error using cell/unique (line 85)
%     Cell array input must be a cell array of character vectors.
%%%
if sum(stimGroupIDX.spontStimIDX.tifFileList)>1

    % Build a per-tif stim-parameter table for the Spont group: one row
    % per tif, scalar columns for per-tif params (trigDelay, ISI,
    % totalPulses) and a per-pulse SpontMsStimLen cell column. Feeds
    % anmlROIbyStimTable's multi-pulse expansion + dedup.
    tifStimParamTable = stimParams2TifTable(...
        tifFileList.stim(stimGroupIDX.spontStimIDX.tifFileList),dataPath);
    
    % sort ROIs by stim params
    [anmlROIbyStim,stimTable] = anmlROIbyStimTable(animal,...
        tifFileList.stim(stimGroupIDX.spontStimIDX.tifFileList),...
        moCorROI,...
        tifStimParamTable);
    outputTables{end+1} = 'tifStimParamTable';
    outputTables{end+1} = tifStimParamTable;
    outputTables{end+1} = 'anmlROIbyStim';
    outputTables{end+1} = anmlROIbyStim;
    outputTables{end+1} = 'stimTable';
    outputTables{end+1} = stimTable;
    
    save(fullfile(dataPath,[animal '_anmlROI_SpontstimTable.mat']),...
        'anmlROIbyStim','stimTable','tifStimParamTable',...
        'dataPath','-v7.3')
end

%sort contrast change stimulus
if sum(stimGroupIDX.contrastChangeIDX.tifFileList)>1

    % Build a per-tif stim-parameter table for the contrast-change group:
    % one row per tif, scalar columns for per-tif params and cell columns
    % for per-pulse params. Feeds anmlROIbyStimTable's pivot/dedup step.
    dContrastTifParamTable = stimParams2TifTable(...
        tifFileList.stim(stimGroupIDX.contrastChangeIDX.tifFileList),dataPath);
    
    [anmlROIdContrast,dContrastTable] = anmlROIbyStimTable(animal,...
        tifFileList.stim(stimGroupIDX.contrastChangeIDX.tifFileList),...
        moCorROI,...
        dContrastTifParamTable);
    
    outputTables{end+1} = 'dContrastTifParamTable';
    outputTables{end+1} = dContrastTifParamTable;
    outputTables{end+1} = 'anmlROIdContrast';
    outputTables{end+1} = anmlROIdContrast;
    outputTables{end+1} = 'dContrastTable';
    outputTables{end+1} = dContrastTable;
    
    save(fullfile(dataPath,[animal '_anmlROI_dContrastTable.mat']),'anmlROIdContrast','dContrastTable',...
        'dContrastTifParamTable','dataPath','tifFileList','fissaScaleFactor',...
        'stimGroupIDX','pulseLegend2P','-v7.3')
end