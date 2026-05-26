function [tifStimParamTable,fileStruct] = stimParams2TifTable(fileStruct,varargin)
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
% Important
% parameters specific to each TIF file are stored as scalar variables,
% while parameters extracted from the pulse(s) are stored as cell variables

end