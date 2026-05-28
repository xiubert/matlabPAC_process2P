function [anmlROIbyStim,stimTable] = anmlROIbyStimTable(animal,tifFileListStim,moCorROI,tifStimParamTable)
% anmlROIbyStimTable  Build a per-(ROI, unique-stim) response table for
%                     one animal, plus a deduplicated stim parameter table.
%
%   [anmlROIbyStim, stimTable] = anmlROIbyStimTable(animal, ...
%                                tifFileListStim, moCorROI, tifStimParamTable)
%
%   Pivots per-tif fluorescence traces into a long-form table keyed by
%   (animal, roiID, stimID), so all repetitions of a unique stimulus
%   condition for a given ROI end up vertically concatenated in a single
%   row's cell entries. Also returns the deduplicated stimulus table,
%   sorted so 'post' rows precede 'pre' rows on the 'treatment' column.
%
%   Inputs:
%     animal           - animal ID string (e.g. 'AA0067'); replicated into
%                        every output row.
%     tifFileListStim  - struct array (the .stim substruct of tifFileList)
%                        subset to the stim group being processed. Each
%                        element must carry:
%                          frameRate           - scalar (Hz)
%                          rawFroi             - nROI x nFrame
%                          moCorRawFroi        - nROI x nFrame
%                          fissaFroi           - (optional) nROI x nFrame
%                          SCALEDfissaFroi     - (optional) nROI x nFrame
%                        If fissaFroi is present on the struct, FISSA
%                        columns are included in the output.
%     moCorROI         - struct array of motion-corrected ROIs (.ID field
%                        used as the roiID column).
%     tifStimParamTable- per-tif stim parameter table from
%                        stimParams2TifTable. May be single-pulse
%                        (scalar columns) or multi-pulse (cell columns
%                        carrying N-by-1 per-pulse values plus a
%                        totalPulses column with totalPulses>1).
%
%   Outputs:
%     anmlROIbyStim - (nROI * nUniqueStim) x M table:
%                       animal, roiID, stimID, <stim param columns>,
%                       frameRate, rawFroi, moCorRawFroi,
%                       fissaFroi, SCALEDfissaFroi (last two only if
%                       FISSA inputs were present).
%                     The *Froi columns are cells holding nRep x nFrame
%                     matrices (rows = repetitions of the stim on that
%                     ROI, columns = time).
%     stimTable     - nUniqueStim x M table: one row per unique stim
%                     condition, with stim params + frameRate (+ pulseID
%                     if rawPulse was present on the input table).
%                     stimID values in anmlROIbyStim index into this
%                     table after the post-before-pre sort.
%
%   Processing steps:
%     1. Equalize trace lengths: for each of {rawFroi, moCorRawFroi,
%        fissaFroi, SCALEDfissaFroi} present on tifFileListStim, trim all
%        per-tif vectors to the shortest length across the set.
%     2. Strip any 'tif*' columns from tifStimParamTable (provenance
%        columns not used downstream).
%     3. Build tmpT: one row per (tif, pulse) pair. Multi-pulse tifs
%        (totalPulses>1) are expanded by replicating scalar columns and
%        unpacking N-by-1 cell columns into N rows. Trial-level
%        identifiers ('trigDelay','ISI','totalPulses') are not replicated.
%     4. unique(tmpT) -> stimTable + ic mapping. Sort by 'treatment' and
%        flipud so 'post' precedes 'pre' alphabetically. stimIDX is the
%        resulting 1..nUniqueStim index per tif.
%     5. For each stim, gather the cell-array of trial traces per ROI.
%        Multi-pulse tifs are sliced into per-pulse windows using
%        frameRate, trigDelay, and ISI before concatenation; single-pulse
%        tifs are taken whole.
%     6. Stitch the per-ROI cells back onto TanmlROI (which carries
%        animal+roiID+stimID + the replicated stim params).
%
%   Notes:
%     - 'rawPulse' is removed from tmpT before dedup so pulse-index
%       suffixes do not split otherwise-identical stim rows; the
%       per-trial pulse indices are preserved as a pulseID cell column
%       on stimTable.
%     - frameRate is assumed identical across all tifs in tifFileListStim
%       (the function uses tifFileListStim(1).frameRate for every stim).
%     - This function errors if the sort-vs-unique reconciliation in the
%       post/pre flip disagrees ('something went wrong here'); that
%       indicates a unique/sortrows mismatch on the 'treatment' column.

% list of fields inside stim to equalize
fields = {'rawFroi','moCorRawFroi','fissaFroi','SCALEDfissaFroi'};

for fi = 1:numel(fields)
    f = fields{fi};
    % collect lengths for entries that have a non-empty numeric vector in tifFileList(k).stim.(f)
    L = zeros(1,numel(tifFileListStim));
    for k = 1:numel(tifFileListStim)
        L(k) = numel(tifFileListStim(k).(f));   
    end
    m = min(L);  % target length (smallest)
    % Trim each valid vector to length m
    for k = 1:numel(tifFileListStim)
        v = tifFileListStim(k).(f);
        if numel(v) > m
            tifFileListStim(k).(f) = v(1:m);
        end
    end
end

tifStimParamTable = tifStimParamTable(:,~contains(tifStimParamTable.Properties.VariableNames,'tif'));

if isfield(tifFileListStim,'fissaFroi')
        FISSA = true;
else
        FISSA = false;
end

%%
allVars = tifStimParamTable.Properties.VariableNames;
isCell = varfun(@iscell, tifStimParamTable, 'OutputFormat', 'uniform');
cellVars = allVars(isCell);
scalarVars = allVars(~isCell);
tmpT = table();
if ismember('totalPulses', tifStimParamTable.Properties.VariableNames) && any(tifStimParamTable.totalPulses > 1)
    % multiple pulses for each tif; expand the table
    tmpS = struct();
    for v = allVars
        tmpS.(v{1}) = {};
    end
    for r = 1:height(tifStimParamTable)
        % determine how many elements to expand for this row
        % assume all cellVars have same length for this row; use first cellVar
        n = numel(tifStimParamTable.(cellVars{1}){r});
        % append scalar variables repeated n times except 'trigDelay', 'ISI' and 'totalPulses'
        repVars = setdiff(scalarVars, {'trigDelay','ISI','totalPulses'}, 'stable');
        for s = 1:numel(repVars)
            name = repVars{s};
            val = tifStimParamTable.(name)(r);
            replicated = repmat({val}, n, 1);
            tmpS.(name) = [tmpS.(name); replicated];
        end
        % append cell variables by taking contents of the cell (expected n�1)
        for c = 1:numel(cellVars)
            name = cellVars{c};
            cellContents = tifStimParamTable.(name){r}; % should be an n�1 cell
            tmpS.(name) = [tmpS.(name); cellContents];
        end
    end
    % Convert accumulated cell columns back to table
    tableVars = setdiff(allVars, {'trigDelay','ISI','totalPulses'}, 'stable');
    for v = tableVars
        name = v{1};
        col = tmpS.(name);
        % If every element is a numeric (in a cell), convert to numeric column
        if all(cellfun(@(x) isnumeric(x), col))
            tmpT.(name) = cell2mat(col);
        else
            tmpT.(name) = string(col);  % convert to string column
        end
    end 
else% one pulse for each tif
    for v = allVars
        name = v{1};
        col = tifStimParamTable.(name);
        if iscell(col)
            if all(cellfun(@(x) isnumeric(x), col))
                tmpT.(name) = cell2mat(col);
            else
                tmpT.(name) = string(col);
            end
        else
            tmpT.(name) = col;
        end
    end
end
%%
if any(contains(tifStimParamTable.Properties.VariableNames,'rawPulse'))
    tmpT = removevars(tmpT,'rawPulse');
end

% tmpT = convertvars(tmpT, @(x) iscell(x) && isnumeric(x), 'double');
% tmpT = convertvars(tmpT, @(x) iscell(x) && any(cellfun(@ischar, x)), 'string');

[stimTable,~,ic] = unique(tmpT);

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
TanmlROI = table(repmat(string(animal),[length(moCorROI) 1]),string({moCorROI.ID})',...%changed from char to string
    'VariableNames',{'animal','roiID'});%modified since 1 roi caused error
TanmlROI = repmat(TanmlROI,[size(stimTable,1) 1]);
stimID = ones(length(moCorROI),1);
roiTstim = repmat(stimTable(1,:),[length(moCorROI) 1]);
for nStim = 2:size(stimTable,1)
    stimID = [stimID; ones(length(moCorROI),1)*nStim];
    roiTstim = [roiTstim; repmat(stimTable(nStim,:),[length(moCorROI) 1])];
end
TanmlROI = [addvars(TanmlROI,stimID) roiTstim];

%get roiF data by tif file
tifRawF = {};
tifMoCorRawF = {};
if FISSA
    tifFissaFroi = {};
    tifSCALEDfissaFroi = {};
end
for i = 1:length(tifFileListStim)
    if ismember('totalPulses', tifStimParamTable.Properties.VariableNames) && tifStimParamTable.totalPulses(i) > 1
        % multiple pulses for each tif
        framesPreTrig = tifFileListStim(i).frameRate*tifStimParamTable{i,'trigDelay'};
        framesPerPulse = tifFileListStim(i).frameRate*tifStimParamTable{i,'ISI'};
        totalPulse = tifStimParamTable.totalPulses(i);
        tifRawF = [tifRawF;mat2cell(tifFileListStim(i).rawFroi(:,framesPreTrig+1:framesPreTrig+framesPerPulse*totalPulse),length(moCorROI),repmat(framesPerPulse,1,totalPulse))'];
        tifMoCorRawF = [tifMoCorRawF;mat2cell(tifFileListStim(i).moCorRawFroi(:,framesPreTrig+1:framesPreTrig+framesPerPulse*totalPulse),length(moCorROI),repmat(framesPerPulse,1,totalPulse))'];
        if FISSA
            tifFissaFroi = [tifFissaFroi,mat2cell(tifFileListStim(i).fissaFroi(:,framesPreTrig+1:framesPreTrig+framesPerPulse*totalPulse),length(moCorROI),repmat(framesPerPulse,1,totalPulse))'];
            tifSCALEDfissaFroi = [tifSCALEDfissaFroi,mat2cell(tifFileListStim(i).SCALEDfissaFroi(:,framesPreTrig+1:framesPreTrig+framesPerPulse*totalPulse),length(moCorROI),repmat(framesPerPulse,1,totalPulse))'];
        end
    else% one pulse for each tif
        tifRawF = [tifRawF; {tifFileListStim(i).rawFroi}'];
        tifMoCorRawF = [tifMoCorRawF; {tifFileListStim(i).moCorRawFroi}'];
        if FISSA
            tifFissaFroi = [tifFissaFroi; {tifFileListStim(i).fissaFroi}'];
            tifSCALEDfissaFroi = [tifSCALEDfissaFroi; {tifFileListStim(i).SCALEDfissaFroi}'];
        end
    end
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
     
    frameRate(nStim) = tifFileListStim(1).frameRate;% assume same frameRate for all tif files
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
TanmlROI = addvars(TanmlROI,FR,'NewVariableNames',{'frameRate'});

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