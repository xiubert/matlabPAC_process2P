function [tifStimParamTable,fileStruct] = stimParams2TifTable(fileStruct,varargin)
% stimParams2TifTable  Build a per-tif stim parameter table for one stim
%                      group by loading and parsing each tif's _Pulses.mat.
%
%   [tifStimParamTable, fileStruct] = stimParams2TifTable(fileStruct)
%   [tifStimParamTable, fileStruct] = stimParams2TifTable(fileStruct, dataDir)
%
%   For each tif in fileStruct, loads the matching <tif>_Pulses.mat,
%   dispatches to the appropriate stim-family extractor based on the
%   pulse(1).pulseset string, and assembles one row per tif into
%   tifStimParamTable. The output is the input table to
%   anmlROIbyStimTable downstream.
%
%   Inputs:
%     fileStruct - struct array (typically tifFileList.stim filtered to a
%                  single stim group). Required fields:
%                    name      - either the .tif filename or the
%                                _Pulses.mat filename (auto-converted)
%                    folder    - directory containing the _Pulses.mat
%                    treatment - treatment label (e.g. 'pre' / 'post')
%                  On output, .name holds the _Pulses.mat filename and a
%                  new .tif field holds the original .tif filename.
%     dataDir    - (optional) fallback directory used if fileStruct(i).folder
%                  does not resolve (handles
%                  MATLAB:load:couldNotReadFile).
%
%   Output:
%     tifStimParamTable - nTif x M table. Two trailing columns are
%                         appended regardless of stim family:
%                           tif       - original tif filename
%                           treatment - treatment label
%                         Other columns come from the stim-family
%                         extractor (see Dispatch below).
%
%   Column-type convention (preserved across all stim families):
%     - Scalar columns  : parameters that apply to the whole tif
%                         (e.g. trigDelay, ISI, totalPulses).
%     - Cell columns    : per-pulse parameters for multi-pulse pulsesets
%                         (N-by-1 cell values, one entry per pulse).
%     anmlROIbyStimTable relies on this convention to expand multi-pulse
%     tifs into one row per (tif, pulse).
%
%   Dispatch (matched against pulse(1).pulseset substring):
%     'PTinContrast'   -> extractStimParams           (PT / DRC)
%     'contrastChange' -> extractContrastChangeParams
%     'BPN'            -> extractBPNStimParams
%     'spont'          -> extractSpontParams
%
%   Notes:
%     - The first pulseset substring match wins; pulseset strings should
%       be designed so the four tokens are mutually exclusive.
%     - All tifs in fileStruct are assumed to share the same stim family.
%       The current loop body re-resolves stimParamExtractFcn from the
%       FIRST file's pulseset rather than the current file's, so passing
%       a mixed-family fileStruct will silently apply the wrong extractor
%       to later tifs.
%     - Empty per-tif fields are coerced to NaN via emptyCells2NaN before
%       table assembly so heterogeneous rows can still be stacked.

if ~isempty(varargin)
    dataDir = varargin{1};
end

%get pulse files
if all(contains({fileStruct.name}','.tif'))
    tmp0 = {fileStruct.name};
    [fileStruct.tif] = tmp0{:};
    tmp = cellfun(@(c) strrep(c,'.tif','_Pulses.mat'),{fileStruct.name},'UniformOutput',0);
    [fileStruct.name] = tmp{:};
    clear tmp*
end

%load first pulse file
try
    tempStim = load(fullfile(fileStruct(1).folder,fileStruct(1).name));
catch ME %find via MException.last
    if (strcmp(ME.identifier,'MATLAB:load:couldNotReadFile'))
        tempStim = load(fullfile(dataDir,fileStruct(1).name));
    else
        rethrow(ME)
    end
end

% determine pulse type for respective stimParam extract function
if contains(tempStim.pulse(1).pulseset,'PTinContrast')
    stimParamExtractFcn = @extractStimParams;
elseif contains(tempStim.pulse(1).pulseset,'contrastChange')
    stimParamExtractFcn = @extractContrastChangeParams;
elseif contains(tempStim.pulse(1).pulseset,'BPN')
    stimParamExtractFcn = @extractBPNStimParams;
elseif contains(tempStim.pulse(1).pulseset,'spont')
    stimParamExtractFcn = @extractSpontParams;
end

% initialize table with first pulse
tempParams = stimParamExtractFcn(tempStim.params, tempStim.pulse);
tifStimParamTable = cell2table(emptyCells2NaN(struct2cell(tempParams)'),'VariableNames',fieldnames(tempParams));

% load subsequent pulses
for nTif = 2:length(fileStruct)
    try
        stim = load(fullfile(fileStruct(nTif).folder,fileStruct(nTif).name));
        if contains(tempStim.pulse(1).pulseset,'PTinContrast')
            stimParamExtractFcn = @extractStimParams;
        elseif contains(tempStim.pulse(1).pulseset,'contrastChange')
            stimParamExtractFcn = @extractContrastChangeParams;
        elseif contains(tempStim.pulse(1).pulseset,'BPN')
            stimParamExtractFcn = @extractBPNStimParams;
        elseif contains(tempStim.pulse(1).pulseset,'spont')
            stimParamExtractFcn = @extractSpontParams;
        end
        
    catch ME %find via MException.last
        if (strcmp(ME.identifier,'MATLAB:load:couldNotReadFile'))
            stim = load(fullfile(dataDir,fileStruct(nTif).name));
        else
            rethrow(ME)
        end
    end
    
    try
        tifStimParamTable(nTif,:) = cell2table(emptyCells2NaN(struct2cell(stimParamExtractFcn(stim.params, stim.pulse))'));
    catch
        error('can''t parse stim params')
    end
    clear stim
end

tifStimParamTable.tif = convertCharsToStrings({fileStruct.tif}');
tifStimParamTable.treatment = convertCharsToStrings({fileStruct.treatment}');

end