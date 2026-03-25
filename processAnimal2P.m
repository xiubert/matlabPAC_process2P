%% PIPELINE: NoRMCorre (motion correction) and FISSA
clearvars;close all;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  EDIT HERE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataPath = 'C:\Users\JIC402\OneDrive - University of Pittsburgh\Data\AA0051\2P';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DONE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfolder(dataPath)
    error('Data path not found')
end
dList = dir(dataPath);
try
    animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');
    if isempty(animal)
        error('')
    end
catch
    animal = char(inputdlg('Enter animal ID:','Animal Input',[1 35],{'AA0000'}));
end

%% Get list of tif files for analysis and treatment info
% DO NOT RUN IF RUNNING AGAIN FOR A DIFFERENT CELL TYPE

tifFiles = dir(fullfile(dataPath,[animal '*.tif']));
treatmentName = char(inputdlg('Enter treatment name: (cancel if none)','Treatment Input',[1 55],{'ZX1'}));
treatment = cell(length(tifFiles),1);
if ~isempty(treatmentName)
    preIDX = false(length(tifFiles),1);
    selIDX = listdlg('PromptString','Select pre treatment tif files',...
        'ListString',{tifFiles.name}','SelectionMode','multiple');
    preIDX(selIDX)= true;
    prePost = {['pre' treatmentName],['post' treatmentName]};
    treatment(preIDX)= {prePost{1}};
    treatment(~preIDX)= {prePost{2}};
    clear treatmentName preIDX prePost
else
    treatment(cellfun(@isempty,treatment)) = {'none'};
end

disp('Likely Map Files: ')
disp(string({tifFiles(extractfield(tifFiles,'bytes')>11000000).name})')
% check file size, disp('warning')
FRAmapIDX = false(length(tifFiles),1);
FRAmapSelIDX = listdlg('PromptString',{'Select tif files for BF mapping','(cancel if none)'},...
    'ListString',{tifFiles.name}','SelectionMode','multiple');
FRAmapIDX(FRAmapSelIDX) = true;
clear FRAmapSelIDX
if any(FRAmapIDX)
    treatment(FRAmapIDX) = cellfun(@(c) horzcat(c,' FRAmap'),treatment(FRAmapIDX),'uni',0);
end
%%%%
BPNIDX = false(length(tifFiles),1);
BPNSelIDX = listdlg('PromptString',{'Select tif files for BPN','(cancel if none)'},...
    'ListString',{tifFiles.name}','SelectionMode','multiple');
BPNIDX(BPNSelIDX) = true;
clear BPNSelIDX
if any(BPNIDX)
    treatment(BPNIDX) = cellfun(@(c) horzcat(c,' BPN'),treatment(BPNIDX),'uni',0);
end
%%%%
[tifFiles.treatment] = treatment{:};

save(fullfile(dataPath,[animal '_tifFileLegend.mat']),'tifFiles')

%% Split tifs into condition groups to be motion corrected separately
% DO NOT RUN IF RUNNING AGAIN FOR A DIFFERENT CELL TYPE

%In case of REDO:
% clear Ycon rawCatImg options_nonrigid NoRMCorreParams moCorrImgNonRigid shifts2 template2
% load(fullfile(dataPath,[animal '_tifFileLegend.mat']),'tifFiles')

%tifFiles filter
resp = inputdlg('Enter comma separated list of treatment filters for tif files (eg: preZX1, postZX1). NOTE: CASE SENSITIVE. MUST MATCH TREATMENT IN tifFiles.treatment',...
    'Split motion correction by treatment?',[1 90]);
if ~isempty(resp)
    filters = strsplit(string(resp),', ');
    locLogical = false(length(tifFiles),1);
    for k = 1:length(filters)
        tifList.(filters{k}) = tifFiles(contains({tifFiles.treatment}',filters{k}));
        locLogical(contains({tifFiles.treatment}',filters{k}))=1;
    end
    if sum(locLogical)<length(tifFiles)
        tifList.remaining = tifFiles(~locLogical);
        filters{2} = 'remaining';
    end
else
    tifList.all = tifFiles;
end
save(fullfile(dataPath,[animal '_tifCondSplitLegend.mat']),'tifList')

%% perform motion correction w/ NoRMCorre
% DO NOT RUN IF RUNNING AGAIN FOR A DIFFERENT CELL TYPE

moCorN = fieldnames(tifList);
outputPath = fullfile(dataPath,'NoRMCorred');
if ~(exist(outputPath,'dir')==7)
    mkdir(outputPath)
end

%check if .tifs have multiple channels, if so, split them and rearrange
%files
[img,tmpHeader] = readSCIMtif(fullfile(tifList.(moCorN{1})(1).folder,...
    tifList.(moCorN{1})(1).name));
% %{
if isfield(tmpHeader,'hChannels') && ...
        numel(tmpHeader.hChannels.channelSave)>1 && ...
        isstruct(img)
    splitTifs = cell(length(tifFiles),1);
    for tifM = 1:length(tifFiles)
        splitTifs(tifM) = splitTifChans(fullfile(tifFiles(tifM).folder,tifFiles(tifM).name),2);
    end
    if ~(exist(fullfile(tifFiles(1).folder,'rawMergedTifs'),'dir')==7)
        mkdir(fullfile(tifFiles(1).folder,'rawMergedTifs'))
    end
    for tifM = 1:length(tifFiles)
        movefile(fullfile(tifFiles(tifM).folder,tifFiles(tifM).name),...
            fullfile(tifFiles(1).folder,'rawMergedTifs',tifFiles(tifM).name))
    end
    for tifM = 1:length(tifFiles)
        movefile(splitTifs{tifM},...
            fullfile(tifFiles(tifM).folder,tifFiles(tifM).name))
    end
end
%}

%motion correction
for k = 1:length(moCorN)
    %concatenate tifs
    [Ycon,~] = concatenate_files(tifList.(moCorN{k}));
    rawCattemp = single(Ycon); % convert to single precision
    rawCatImg.(moCorN{k}) = rawCattemp - min(rawCattemp(:));
    clear Ycon
    % non-rigid motion correction (in parallel)
    gcp;
    options_nonrigid = NoRMCorreSetParms('d1',size(rawCatImg.(moCorN{k}),1),...
        'd2',size(rawCatImg.(moCorN{k}),2),...
        'grid_size',[32,32],'mot_uf',4,'bin_width',200,'max_shift',15,...
        'max_dev',3,'us_fac',50,'init_batch',200);
    
    tic; [moCorrImgNonRigid.(moCorN{k}),NoRMCorreParams.(moCorN{k}).shifts,...
        ~,NoRMCorreParams.(moCorN{k}).options_nonrigid] = ...
        normcorre_batch(rawCatImg.(moCorN{k}),options_nonrigid); toc
    
    % generate motion corrected tifs from concatenated motion corrected image data for FISSA
    writeMoCorTifs(tifList.(moCorN{k}),moCorrImgNonRigid.(moCorN{k}),outputPath)
    
    clear options_nonrigid rawCattemp
end
save(fullfile(dataPath,'NoRMCorred',[animal '_NoRMCorreParams.mat']),'NoRMCorreParams')

%% In case of REDO or RESUME after motion correction:
% RUN THIS IF RUNNING AGAIN FOR A DIFFERENT CELL TYPE

load(fullfile(dataPath,[animal '_tifCondSplitLegend.mat']),'tifList')
load(fullfile(dataPath,[animal '_tifFileLegend.mat']),'tifFiles')
% FRAmapIDX = contains({tifFiles.treatment},'map')';

moCorN = fieldnames(tifList);
for k = 1:length(moCorN)
    moCorrImgNonRigid.(moCorN{k}) = loadNoRMCorrNonRigidImgViaTifs(dataPath,tifList.(moCorN{k}));
    try
        [Ycon,~] = concatenate_files(tifList.(moCorN{k}));
    catch
        [tifList.(moCorN{k}).folder] = deal(dataPath);
        [Ycon,~] = concatenate_files(tifList.(moCorN{k}));
    end
    rawCattemp = single(Ycon); % convert to single precision
    rawCatImg.(moCorN{k}) = rawCattemp - min(rawCattemp(:));
    clear rawCattemp
end


%% FOR EACH treatment condition (pre/post): 
%1. RUN THIS CELL to Draw ROI from motion corrected concatenated stack
%2. RUN NEXT CELL to save ROI to file

moCorSeqN = listdlg('PromptString','Select treatment condition (pre/post) for which to draw ROI',...
    'ListString',moCorN);
%outputs ROI to workspace and [animalID]_moCorrROI.mat into dataPath for use w/ FISSA
TIFcatROIgui(moCorrImgNonRigid.(moCorN{moCorSeqN}))

%% Save ROIs:  RUN AFTER DRAWING ROIs FOR EACH TREATMENT CONDITION

nTifs = length(tifList.(moCorN{moCorSeqN}));
tifIDXinAllTifList = ismember({tifFiles.name}',{tifList.(moCorN{moCorSeqN}).name}');

save([dataPath filesep ...
    regexp(dataPath,'[A-Z]{2}\d{4}','match','once') ...
    '_moCorrROI_' moCorN{moCorSeqN} '.mat'],...
    'moCorROI','moCorSeqN','nTifs','tifIDXinAllTifList')
disp([regexp(dataPath,'[A-Z]{2}\d{4}','match','once') ...
    '_moCorrROI_' moCorN{moCorSeqN} '.mat saved to animal directory'])
clear nTifs tifIDXinAllTifList

%% Match ROI pre & post (after confirming all ROI for all tif sequences)
%RUN THIS CELL ONLY AFTER RUNNING ABOVE 2 CELLS FOR EACH TREATMENT CONDITION

%ensure true ROI are same pre and post
intersectROIfiles(dataPath,animal,moCorN,tifList,tifFiles)

%% AFTER ROI ARE SAVED:
%Add rawF, moCorr rawF, nFrames and frameRate to tifFiles struct

for ROIfileN = 1:length(moCorN)
    clear moCorROI
    temp = load(fullfile(dataPath,[animal '_moCorrROI_' moCorN{ROIfileN} '.mat']),...
        'moCorROI');
    moCorROI = temp.moCorROI;
    clear temp
    
    %outputs tifFiles struct w/ rawF and moCorr rawF, adds nFrames and frameRate
    tifList.(moCorN{ROIfileN}) = ...
        moCorRawF2tifList(tifList.(moCorN{ROIfileN}),...
        moCorrImgNonRigid.(moCorN{ROIfileN}),...
        moCorROI,...
        rawCatImg.(moCorN{ROIfileN}));
end
clear moCorrImgNonRigid rawCatImg moCorROI

%Save NoRMCorre params and file names of tifs used
allTifFiles = tifFiles;
try
    save([dataPath filesep animal '_moCorr_Tifs_Params.mat'],'tifList','NoRMCorreParams','allTifFiles','-v7.3')
catch
    save([dataPath filesep animal '_moCorr_Tifs_Params.mat'],'tifList','allTifFiles','-v7.3')
end

%% NEUROPIL CORRECTION via FISSA: FISSAviaMatlab_prePostTreatment.py

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.  edit animalDataPath line accordingly
% 2.  run 'python FISSAviaMatlab_prePostTreatment.py' in python environment w/ fissa installed
%FISSA basically just takes ROIs and folder of motion corrected tifs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Separate map data from tifFiles and FISSAoutput, add FISSAoutput to tifFileList 

%load FISSA output: (again from exp.save_to_matlab() in 'FISSAscript_splitPrePost.py')
%FISSA was given moCorr tiffs and ROIs from concatenated moCorr Tifs
FISSAoutput = load(fullfile(dataPath,'NoRMCorred','FISSAoutput','matlab.mat'));

%in 'result', for a given cell and trial there is a n x numTraceFrames double;
%row 1 is ROI trace, rows 2->n are traces from neuropil regions around ROI
%in 'ROIs', for a given cell and trial there is a n x 1 cell of doubles,
%the doubles are nPoints x 2, col1=Y,col2=X; 1st cell is ROI,  2->n are neuropil regions around ROI

%split map and stim from FISSAoutput
if ~exist('FRAmapIDX')
    FRAmapIDX = contains({tifFiles.treatment}','map');
end
if ~exist('BPNIDX')
    BPNIDX = contains({tifFiles.treatment}','BPN');
end

if any(FRAmapIDX) || any(BPNIDX)
    if length(fieldnames(tifList))>1
        C = struct2cell(tifList);
        tifFiles = vertcat(C{:});
    else
        tifFiles = tifList.all;
    end
    tifFileList.map = tifFiles(FRAmapIDX);
    tifFileList.BPN = tifFiles(BPNIDX);
    tifFileList.stim = tifFiles(~FRAmapIDX & ~BPNIDX);
    
    fID = fieldnames(FISSAoutput);
    trials.all = fieldnames(FISSAoutput.raw.cell0);
    trials.BF = strcat(cellstr(repmat('trial',[sum(FRAmapIDX) 1])),cellstr(string(0:sum(FRAmapIDX)-1))');
    trials.BPN = strcat(cellstr(repmat('trial',[sum(BPNIDX) 1])),cellstr(string(0:sum(BPNIDX)-1))');
    trials.stim = strcat(cellstr(repmat('trial',[sum(~FRAmapIDX & ~BPNIDX) 1])),cellstr(string(0:sum(~FRAmapIDX & ~BPNIDX)-1))');
    for nField = 1:length(fID)
        for nCell = 1:length(fieldnames(FISSAoutput.(fID{nField})))
            tmp = struct2cell(FISSAoutput.(fID{nField}).(['cell' num2str(nCell-1)]));
            FISSAout.map.(fID{nField}).(['cell' num2str(nCell-1)]) = ...
                cell2struct(tmp(FRAmapIDX),trials.BF);
            FISSAout.BPN.(fID{nField}).(['cell' num2str(nCell-1)]) = ...
                cell2struct(tmp(BPNIDX),trials.BPN);
            FISSAout.stim.(fID{nField}).(['cell' num2str(nCell-1)]) = ...
                cell2struct(tmp(~FRAmapIDX & ~BPNIDX),trials.stim);
            clear tmp
        end
    end
    clear FISSAoutput
    FISSAoutput.map = FISSAout.map;
    FISSAoutput.BPN = FISSAout.BPN;
    FISSAoutput.stim = FISSAout.stim;
    clear tmp FISSAout trials fID
elseif isfield(tifList,'all')
    tifFileList.stim = tifList.all;
    tmp = FISSAoutput;
    clear FISSAoutput
    FISSAoutput.stim = tmp;
    clear tmp
else
    if length(fieldnames(tifList))>1
        C = struct2cell(tifList);
        tifFileList.stim = vertcat(C{:});
    else
        tifFileList.stim = tifFiles;
    end
    tmp = FISSAoutput;
    clear FISSAoutput
    FISSAoutput.stim = tmp;
    clear tmp
end

%if updating existing tifFileList w/ FISSAoutput:
% [fName,fPath] = uigetfile('*tifFileList.mat','Locate [ANIMAL]_tifFileList.mat...');
% load(fullfile(fPath,fName))
% animal = regexp(dataPath,'[A-Z]{2}\d{4}','match');
% animal = animal{1};

fissaScaleFactor = str2double(inputdlg('ENTER FACTOR BY WHICH TO SCALE FISSA SUBTRACTION (eg. 0.8): ',...
    'FISSA SCALING FACTOR',[1 80],{'0.8'}));
tifFileList = FISSAoutput2tifFileList(FISSAoutput,tifFileList,fissaScaleFactor);

save(fullfile(dataPath,[animal '_tifFileList.mat']),...
        'dataPath','FISSAoutput','tifFileList','fissaScaleFactor','-v7.3')

%% COMPLETE. YOU ARE NOW LEFT WITH MOTION AND NEUROPIL CORRECTED FLUORESCENCE TRACES FOR EACH ROI FOR EACH TIF

% eg. tifFileList.map.SCALEDfissaFroi is in shape ROI x fluoresence 
% at corresponding .tif frame

%see stimParam2ROI.m if have stim params in _Pulses.mat files
[pulseLegend2P,stimGroupIDX,ROIoutputTables] = stimParam2ROI(dataPath);
%% Clear Workspace first if continuing AFTER GETTING tifFileList.map.SCALEDfissaFroi

if ~exist('dataPath','var')
    dataPath='C:\Users\JIC402\OneDrive - University of Pittsburgh\Data\AA0044';
    animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');
    load(fullfile(dataPath,[animal '_anmlROI_stimTable.mat']))
end
if ~exist('nCell','var')
    nCell=size(tifFileList.stim(1).rawFroi,1);
end
if exist('ROIoutputTables','var')
    tifStimParamTable=ROIoutputTables{2};
    anmlROIbyStim=ROIoutputTables{4};
    stimTable=ROIoutputTables{6};
end
%% FRAmap
%%%%%%%%%%%%%%%%%%%%%%%%%% EDIT IF NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkPTsigSD = 2;
nFramesPostPulse = 2;
plotAllROI = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MapOutputDir = fullfile(dataPath,'BFMap');
if ~isfolder(MapOutputDir)
    mkdir(MapOutputDir)
end 

adHocFRAmap = FRAmap(tifFileList,pkPTsigSD,nFramesPostPulse,'SCALEDfissaFroi');

%% PLOT OUTPUT | Sig Responses

plotFRAmap(adHocFRAmap,'plotAllROI',plotAllROI)

figure;
semilogx(adHocFRAmap.freqList,nanmean(adHocFRAmap.uSigPkResp,1));
[~, maxID] = max(nanmean(adHocFRAmap.uSigPkResp,1));

disp(adHocFRAmap.freqList(maxID))

%% PLOT OUTPUT | All Responses

plotFRAmap(adHocFRAmap,'sigResp',false,'plotAllROI',plotAllROI)

figure;
semilogx(adHocFRAmap.freqList,nanmean(adHocFRAmap.uPkResp,1));
xlabel('Frequency/Hz')
ylabel('dF/F')
[~, maxID] = max(nanmean(adHocFRAmap.uPkResp,1));
disp(adHocFRAmap.freqList(maxID))

%% PLOT RESP to CGC
% calculate average F
anmlROIbyStim.t_total = rowfun(@(F,fr,trigDelay) ...
    {(1:size(F,2))/fr-trigDelay},...
    anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi','frameRate','trigDelay'},...
    'ExtractCellContents',true,'OutputFormat','uniform');
anmlROIbyStim.F_avg = rowfun(@(F) ...
    {mean(F)},...
    anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

%% Plot rawF
anmlROIbyStim.roiID = string(strtrim(cellstr(anmlROIbyStim.roiID)));
roiList = unique(anmlROIbyStim.roiID, 'stable');
ROIperFig = 9;
remROIplotNo = rem(nCell,ROIperFig);
roiFigNo = floor(nCell/ROIperFig)+(remROIplotNo>=1);
dBdeltaList = unique(stimTable.dBdelta);
ndBdelta=height(dBdeltaList);
%initialize subplots for multiple ROI per fig
for roiFigN = 1:roiFigNo
    figure;
    title('raw F responses for each ROI')
    for roiSubPlotN = 1:ROIperFig
        curROIno = roiFigN*ROIperFig - ROIperFig + roiSubPlotN;
        if curROIno <= nCell
            subplot(3,3,roiSubPlotN);
            roi=roiList(curROIno);
            rows=anmlROIbyStim(anmlROIbyStim.roiID==roi,:);
            hold on;
            label=strings(height(rows));
            for r=1:height(rows)
                x=rows.t_total{r};
                y=rows.F_avg{r};
                plot(x,y);
                if rows.dBdelta(r) == dBdeltaList(1)
                    label(r)='Low contrast';
                elseif rows.dBdelta(r) == dBdeltaList(2)
                    label(r)='High contrast';
                end
            end
            xlabel('time/s')
            ylabel('raw F')
            xline(2,'--','pure tone')
            hold off;
            title('ROI'+string(curROIno));
            legend(label(1), label(2))
        end
        clear curROIno
    end
    clear roiSubPlotN
end
%% 

tBaseDRC = [-1.2 0];
tBasePT = [1 2];

%dFF re DRC

anmlROIbyStim.t_dFF_DRC = rowfun(@(t) ...
    {t(find((t>=tBaseDRC(1)),1,'first'):end)},...
    anmlROIbyStim,'InputVariables',{'t_total'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

anmlROIbyStim.dFF_DRC = rowfun(@(F,t) ...
    {dFoFcalc(F,[find((t>=tBaseDRC(1)),1,'first')...
    find((t<=tBaseDRC(2)),1,'last')],1)},...
    anmlROIbyStim,'InputVariables',{'F_avg','t_total'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

% dFF_DRC_high = dFoFcalc(F_high_avg,[find((time_vector>=tBaseDRC(1)),1,'first')...
%      find((time_vector<=tBaseDRC(2)),1,'last')],1);
% dFF_DRC_low = dFoFcalc(F_low_avg,[find((time_vector>=tBaseDRC(1)),1,'first')...
%      find((time_vector<=tBaseDRC(2)),1,'last')],1);
% t_dFF_DRC = time_vector(find((time_vector>=tBaseDRC(1) & time_vector<=tBaseDRC(2)),1,'first'):end);

%initialize subplots for multiple ROI per fig
for roiFigN = 1:roiFigNo
    figure;
    title('dF/F responses for each ROI')
    for roiSubPlotN = 1:ROIperFig
        curROIno = roiFigN*ROIperFig - ROIperFig + roiSubPlotN;
        if curROIno <= nCell
            subplot(3,3,roiSubPlotN);
            roi=roiList(curROIno);
            rows=anmlROIbyStim(anmlROIbyStim.roiID==roi,:);
            hold on;
            label=strings(height(rows));
            for r=1:height(rows)
                x=rows.t_dFF_DRC{r};
                y=rows.dFF_DRC{r};
                plot(x,y);
                if rows.dBdelta(r) == dBdeltaList(1)
                    label(r)='Low contrast';
                elseif rows.dBdelta(r) == dBdeltaList(2)
                    label(r)='High contrast';
                end
            end
            xlabel('time/s')
            ylabel('dF/F')
            xline(2,'--')
            hold off;
            title('ROI'+string(curROIno));
            legend(label(1), label(2),'pure tone')
        end
        clear curROIno
    end
    clear roiSubPlotN
end

% dFF_PT_preDRCf0_high = dFF_DRC_high  - nanmean(dFF_DRC_high(:,t_dFF_DRC>=tBasePT(1) & t_dFF_DRC<=tBasePT(2)),2)
% dFF_PT_preDRCf0_low = dFF_DRC_low  - nanmean(dFF_DRC_low(:,t_dFF_DRC>=tBasePT(1) & t_dFF_DRC<=tBasePT(2)),2)
%
%% dFF re PT re DRCf0
anmlROIbyStim.t_dFF_PT = rowfun(@(t) ...
    {t(find((t>=tBasePT(1)),1,'first'):end)},...
    anmlROIbyStim,'InputVariables',{'t_total'},...
    'ExtractCellContents',true,'OutputFormat','uniform');
anmlROIbyStim.dFF_PT = rowfun(@(F,t) ...
    {dFoFcalc(F,[find((t>=tBasePT(1)),1,'first')...
    find((t<=tBasePT(2)),1,'last')],1)},...
    anmlROIbyStim,'InputVariables',{'F_avg','t_total'},...
    'ExtractCellContents',true,'OutputFormat','uniform');
anmlROIbyStim.dFF_PT_preDRCf0 = rowfun(@(dFF_DRC,t) ...
    {dFF_DRC  - ...
    nanmean(dFF_DRC(:,t>=tBasePT(1) & t<=tBasePT(2)),2)},...
    anmlROIbyStim,'InputVariables',{'dFF_DRC','t_dFF_DRC'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

%% 
%initialize subplots for multiple ROI per fig
for roiFigN = 1:roiFigNo
    figure;
    title('dF/F responses re PT re DRCf0 for each ROI');
    for roiSubPlotN = 1:ROIperFig
        curROIno = roiFigN*ROIperFig - ROIperFig + roiSubPlotN;
        if curROIno <= nCell
            subplot(3,3,roiSubPlotN);
            roi=roiList(curROIno);
            rows=anmlROIbyStim(anmlROIbyStim.roiID==roi,:);
            hold on;
            label=strings(height(rows));
            for r=1:height(rows)
                x=rows.t_dFF_DRC{r};
                y=rows.dFF_PT_preDRCf0{r};
                plot(x,y,'LineWidth', 2);
                if rows.dBdelta(r) == dBdeltaList(1)
                    label(r)='Low contrast';
                elseif rows.dBdelta(r) == dBdeltaList(2)
                    label(r)='High contrast';
                end
            end
            xlabel('time/s')
            ylabel('dF/F')
            xline(2,'--')
            hold off;
            title('ROI'+string(curROIno));
            legend(label(1), label(2),'pure tone')
        end
        clear curROIno
    end
    clear roiSubPlotN
end
%% Plot avg of all ROIs

colors.lohi = [0,0.451000000000000,0.741200000000000;0.851000000000000,0.329400000000000,0.102000000000000];
colors.lohiTrace = [0.729400000000000,0.874500000000000,1;1,0.694100000000000,0.541200000000000];
[groups,idC] = findgroups(anmlROIbyStim.dBdelta);
tmp = splitapply(@(x) {vertcat(x{:})},anmlROIbyStim.dFF_PT,groups);
[G0,idC0] = findgroups(idC);
dFF_PT_avg = splitapply(@(x) {cellfun(@nanmean,x,'uni',0)},tmp,G0);
dFF_PT_sem = splitapply(@(x) {cellfun(@SEMcalc,x,'uni',0)},tmp,G0);

figure;
hold on
for r=1:ndBdelta
    x=anmlROIbyStim(1,:).t_dFF_PT{1,1};
    y=dFF_PT_avg{r,1}{1,1};
    yerr=dFF_PT_sem{r,1}{1,1};
    fillSEMplot(x,y,yerr,colors.lohi(r,:),colors.lohiTrace(r,:));
end


xlabel('time/s')
ylabel('dF/F')
xline(2,'--','pure tone')
xlim([1 5])
hold off;
title('Average across cell');
legend('Low contrast', 'High contrast')
%% Peak dFF
pkPTframeBin = 4;
pkPTsigSD=1.5;

tmp = rowfun(@(dFF,t,PTonset,fr) ...
    pkFcalc(dFF,...
    find(t>=(PTonset+(1/fr)),1,'first'),...
    pkPTframeBin,pkPTsigSD),...
    anmlROIbyStim,'InputVariables',{'dFF_PT_preDRCf0','t_dFF_DRC','PTsOnset','frameRate'},...
    'ExtractCellContents',true,'OutputFormat','cell','OutputVariableNames',{'sigPk','sig','pk'});
anmlROIbyStim.pkPT_sig = tmp(:,1);
anmlROIbyStim.sigPk = tmp(:,2);
anmlROIbyStim.pkPT = tmp(:,3);
%% Save
save(fullfile(dataPath,[animal '_anmlROI_stimTable.mat']),"anmlROIbyStim",'-append');
%% Low vs High per ROI

x = nan(nCell,1);
y = nan(nCell,1);
for i = 1:nCell
    roi=roiList(i);
    rows=anmlROIbyStim(anmlROIbyStim.roiID==roi,:);
    x(i)=cell2mat(rows.pkPT(rows.dBdelta==dBdeltaList(2)));
    y(i)=cell2mat(rows.pkPT(rows.dBdelta==dBdeltaList(1)));
end

% keep only roiIDs with both values positive
valid = (x>0) & (y>0);
x = x(valid); y = y(valid);
roiList_pos = roiList(valid);

% make scatter
figure;
scatter(x,y,45,'filled','MarkerFaceAlpha',0.8);
hold on;
% identity line
lims = [0 1];
plot(lims, lims, '--k', 'LineWidth', 1);
hold off;

xlabel('High contrast');
ylabel('Low contrast');
title('peak dF/F per roi');
% grid on;
axis equal;
xlim(lims); ylim(lims);

%% Bar graph
group=cell(ndBdelta,1);
pkResp_means=NaN(ndBdelta,1);
pkResp_sems=NaN(ndBdelta,1);

for k = 1:ndBdelta
    vals = cell2mat(anmlROIbyStim.pkPT(groups == k));
    vals=vals(valid);
    group{k}=vals;
    pkResp_means(k) = mean(vals,'omitnan');
    pkResp_sems(k)= std(vals)/sqrt(nCell);
end

% Create bar plot and error bars (SEM)
figure;
b = bar(1:ndBdelta, pkResp_means,'FaceColor','flat');
hold on;
errorbar(1:ndBdelta, pkResp_means, pkResp_sems, 'k.', 'LineWidth',1);

% Overlay individual points with jitter
rng(0); % for reproducible jitter
jitterAmount = 0.08; % tweak for spread
for k = 1:ndBdelta
    b.CData(k,:)=colors.lohi(k,:);
    vals = cell2mat(anmlROIbyStim.pkPT(groups == k));
    vals=vals(valid);
    x = (k) + (rand(size(vals)) - 0.5) * 2 * jitterAmount;
    % Plot points
    scatter(x, vals, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
end

% Formatting

xticklabels({'Low contrast','High contrast'});             % works for numeric/categorical/string
ylabel('peak dF/F');
% title('pkPT by dBdelta (individual points and group mean)');
box on;

if ndBdelta == 2
    [h,p,ci,stats] = ttest(group{1}, group{2});
    means=[pkResp_means(1),pkResp_means(2)];
    sems=[pkResp_sems(1),pkResp_sems(2)];
    % Add significance star or text
    yMax = max([means + sems]) ;
    yStar = yMax + 0.5*range([means sems]) ;  % vertical position for star/line
    % Draw bar connecting line
    plot([1 2], [yStar yStar], '-k', 'LineWidth',1);
    % Draw short ticks
    plot([1 1], [yStar-0.1*range([means sems]) yStar], '-k', 'LineWidth',1);
    plot([2 2], [yStar-0.1*range([means sems]) yStar], '-k', 'LineWidth',1);

    if p < 0.001
        sigtxt = '***';
    elseif p < 0.01
        sigtxt = '**';
    elseif p < 0.05
        sigtxt = '*';
    else
        sigtxt = sprintf('p=%.3g', p);
    end
    text(1.5, yStar + 0.03*range([means sems]), sigtxt, 'HorizontalAlignment','center', 'FontSize',14);
end

hold off;
