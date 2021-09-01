%% stimParam2ROI
%IF HAVE STIM PARAMS in _Pulses.mat files, organize stim output my animal, ROI, stim params
function [pulseLegend2P,stimGroupIDX,outputTables] = stimParam2ROI(dataPath)
outputTables = cell(0);
%eg for a given ROI all traces from a given stim param will be in one row
animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');

tmp0 = cellstr(ls(fullfile(dataPath,[animal '_moCorrROI_*.mat'])));
tmp = load(fullfile(dataPath,tmp0{1}),...
    'moCorROI');
moCorROI = tmp.moCorROI;
clear tmp*

load(fullfile(dataPath,[animal '_tifFileList.mat']),...
    'tifFileList','FISSAoutput','fissaScaleFactor')

% get pulse name for each tif
pulseLegend2P = tifPulseLegend2P(dataPath);

% sort by stimulus group
stimGroupIDX.ptStimIDX.pulseLegend2P = contains({pulseLegend2P.pulseSet},'PTinContrast')';
stimGroupIDX.ptStimIDX.tifFileList = ismember({tifFileList.stim.name},...
    {pulseLegend2P(stimGroupIDX.ptStimIDX.pulseLegend2P).tif})';
stimGroupIDX.contrastChangeIDX.pulseLegend2P = contains({pulseLegend2P.pulseSet},'contrastChange')';
stimGroupIDX.contrastChangeIDX.tifFileList = ismember({tifFileList.stim.name},...
    {pulseLegend2P(stimGroupIDX.contrastChangeIDX.pulseLegend2P).tif})';
stimGroupIDX.pupilReflexIDX.pulseLegend2P = contains({pulseLegend2P.pulseSet},'LED trigger')';
stimGroupIDX.pupilReflexIDX.tifFileList = ismember({tifFileList.stim.name},...
    {pulseLegend2P(stimGroupIDX.pupilReflexIDX.pulseLegend2P).tif})';

%sort pure tone in contrast stims
if sum(stimGroupIDX.ptStimIDX.tifFileList)>1
    
    % get stim params for each tif
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
    
    save(fullfile(dataPath,[animal '_anmlROI_stimTable.mat']),'anmlROIbyStim','stimTable','tifStimParamTable',...
        'dataPath','FISSAoutput','tifFileList','fissaScaleFactor',...
        'stimGroupIDX','pulseLegend2P','-v7.3')
else
    save(fullfile(dataPath,[animal '_tifFileList.mat']),...
        'dataPath','FISSAoutput','tifFileList','fissaScaleFactor',...
        'stimGroupIDX','pulseLegend2P','-v7.3')
end

%sort contrast change stimulus
if sum(stimGroupIDX.contrastChangeIDX.tifFileList)>1
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