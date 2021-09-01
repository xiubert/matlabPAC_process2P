%% PIPELINE: NoRMCorre (motion correction) and FISSA
clearvars;close all;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  EDIT HERE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataPath = 'C:\Data\sampleData\AA0211';

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

FRAmapIDX = false(length(tifFiles),1);
FRAmapSelIDX = listdlg('PromptString',{'Select tif files for BF mapping','(cancel if none)'},...
    'ListString',{tifFiles.name}','SelectionMode','multiple');
FRAmapIDX(FRAmapSelIDX) = true;
clear FRAmapSelIDX
if any(FRAmapIDX)
    treatment(FRAmapIDX) = cellfun(@(c) horzcat(c,' FRAmap'),treatment(FRAmapIDX),'uni',0);
end
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
%{
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
%}

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

if any(FRAmapIDX)
    if length(fieldnames(tifList))>1
        C = struct2cell(tifList);
        tifFiles = vertcat(C{:});
    else
        tifFiles = tifList.all;
    end
    tifFileList.map = tifFiles(FRAmapIDX);
    tifFileList.stim = tifFiles(~FRAmapIDX);
    
    fID = fieldnames(FISSAoutput);
    trials.all = fieldnames(FISSAoutput.raw.cell0);
    trials.BF = strcat(cellstr(repmat('trial',[sum(FRAmapIDX) 1])),cellstr(string(0:sum(FRAmapIDX)-1))');
    trials.stim = strcat(cellstr(repmat('trial',[sum(~FRAmapIDX) 1])),cellstr(string(0:sum(~FRAmapIDX)-1))');
    for nField = 1:length(fID)
        for nCell = 1:length(fieldnames(FISSAoutput.(fID{nField})))
            tmp = struct2cell(FISSAoutput.(fID{nField}).(['cell' num2str(nCell-1)]));
            FISSAout.map.(fID{nField}).(['cell' num2str(nCell-1)]) = ...
                cell2struct(tmp(FRAmapIDX),trials.BF);
            FISSAout.stim.(fID{nField}).(['cell' num2str(nCell-1)]) = ...
                cell2struct(tmp(~FRAmapIDX),trials.stim);
            clear tmp
        end
    end
    clear FISSAoutput
    FISSAoutput.map = FISSAout.map;
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