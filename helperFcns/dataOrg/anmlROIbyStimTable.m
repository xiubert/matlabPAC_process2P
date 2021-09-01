function [anmlROIbyStim,stimTable] = anmlROIbyStimTable(animal,tifFileListStim,moCorROI,tifStimParamTable)
%sort responses from all ROI for a given stim
%into a table of animal ROIs by stim parameters
tifStimParamTable = tifStimParamTable(:,~contains(tifStimParamTable.Properties.VariableNames,'tif'));

if isfield(tifFileListStim,'fissaFroi')
        FISSA = true;
else
        FISSA = false;
end

if any(contains(tifStimParamTable.Properties.VariableNames,'rawPulse'))
    [stimTable,~,ic] = unique(removevars(tifStimParamTable,'rawPulse'));
else
    [stimTable,~,ic] = unique(tifStimParamTable);
end

%sort stim table by treatment
[sT,sIDX] = sortrows(stimTable,'treatment');
sT = flipud(sT); %bc post alphabetically comes before pre
sIDX = flipud(sIDX);
stimIDX = zeros(size(ic));
for k = 1:length(sIDX)
    %stimIDX is the corresponding stim number for each tif file
    stimIDX(ic==sIDX(k)) = k;
    if ~isequaln(stimTable(sIDX(k),:),sT(k,:))
        error('something went wrong here')
    end    
end
clear stimTable ic sIDX
stimTable = sT;
clear sT

%begin to create table of ROIs for each stim
TanmlROI = table(repmat(animal,[length(moCorROI) 1]),char({moCorROI.ID}),...
    'VariableNames',{'animal','roiID'});
TanmlROI = repmat(TanmlROI,[size(stimTable,1) 1]);
stimID = ones(length(moCorROI),1);
roiTstim = repmat(stimTable(1,:),[length(moCorROI) 1]);
for nStim = 2:size(stimTable,1)
    stimID = [stimID; ones(length(moCorROI),1)*nStim];
    roiTstim = [roiTstim; repmat(stimTable(nStim,:),[length(moCorROI) 1])];
end
TanmlROI = [addvars(TanmlROI,stimID) roiTstim];

%get roiF data by tif file
tifRawF = {tifFileListStim.rawFroi}';
tifMoCorRawF = {tifFileListStim.moCorRawFroi}';
if FISSA
    tifFissaFroi = {tifFileListStim.fissaFroi}';
    tifSCALEDfissaFroi = {tifFileListStim.SCALEDfissaFroi}';
end

%organize pulseID info
if any(contains(tifStimParamTable.Properties.VariableNames,'rawPulse'))
    tifPulseID = str2double(cellfun(@(c) c(end), tifStimParamTable.rawPulse,'uni',0));
    pulseID = cell(size(stimTable,1),1);
end

%init vars
frameRate = zeros(size(stimTable,1),1);
[rawFroi,moCorRawFroi] = deal(cell(size(TanmlROI,1),1));
if FISSA
    [fissaFroi,SCALEDfissaFroi] = deal(cell(size(TanmlROI,1),1));
end

%fill data from each unique stimulus condition into ROIs at respective
%table locations
for nStim = 1:size(stimTable,1) %faster if loop over stim instead of ROI
    
    %get tifs associated with corresponding stim condition
    idxTif = find(nStim==stimIDX);
     
    frameRate(nStim) = tifFileListStim(idxTif(1)).frameRate;
    if any(contains(tifStimParamTable.Properties.VariableNames,'rawPulse'))
        pulseID{nStim} = tifPulseID(idxTif);
    end
    
    %for each ROI fill roiF matrix with respective data from stim
    for nROI = 1:length(moCorROI)
        %table may not include all ROI, ensure ROI corresponds to roiID
        roiIDstr = moCorROI(nROI).ID;
        %find table for for respective ROI and stim
        anmlROIidx = find(and(str2double(string(TanmlROI{:,'roiID'}))==str2double(roiIDstr),...
            TanmlROI{:,'stimID'}==nStim));
        
        %fill roiF data
        rawFroi{anmlROIidx} = cell2mat(cellfun(@(c) c(nROI,:),tifRawF(idxTif),'uni',0));
        moCorRawFroi{anmlROIidx} = cell2mat(cellfun(@(c) c(nROI,:),tifMoCorRawF(idxTif),'uni',0));
        if FISSA
            fissaFroi{anmlROIidx} = cell2mat(cellfun(@(c) c(nROI,:),tifFissaFroi(idxTif),'uni',0));
            SCALEDfissaFroi{anmlROIidx} = cell2mat(cellfun(@(c) c(nROI,:),tifSCALEDfissaFroi(idxTif),'uni',0));
        end
        clear anmlROIidx roiIDstr  
    end
    clear idxTif
end

%add framerate to anmlROIbyStim table
FR = zeros(size(TanmlROI,1),1);
for nStim = 1:size(stimTable,1)
    FR(TanmlROI.stimID==nStim) = frameRate(nStim);
end
TanmlROI = addvars(TanmlROI,FR,'After','trigDelay','NewVariableNames',{'frameRate'});

%add FISSA output
if FISSA
    anmlROIbyStim = addvars(TanmlROI,rawFroi,moCorRawFroi,fissaFroi,SCALEDfissaFroi);
else
    anmlROIbyStim = addvars(TanmlROI,rawFroi,moCorRawFroi);
end

%add pulseID info
if any(contains(tifStimParamTable.Properties.VariableNames,'rawPulse'))
    stimTable = addvars(stimTable,pulseID,frameRate);
else
    stimTable = addvars(stimTable,frameRate);
end