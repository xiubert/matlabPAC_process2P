%% PIPELINE: NoRMCorre (motion correction) and FISSA
%
% Full 2P processing pipeline for a single animal session:
%
%   1. Tif inventory     - list all tifs, assign treatment (pre/post) and
%                          FRA map labels, save tifFiles legend
%   2. Condition split   - group tifs by treatment for independent motion
%                          correction
%   3. Motion correction - split multi-channel tifs if needed, concatenate
%                          per group, run NoRMCorre non-rigid correction,
%                          write corrected tifs for FISSA
%   4. ROI drawing       - interactively draw ROIs on the motion-corrected
%                          stack for each treatment condition
%   5. ROI matching      - reconcile ROIs across pre/post conditions so the
%                          same cells are tracked throughout
%   6. Raw F extraction  - extract rawF and motion-corrected rawF per ROI
%                          per tif, save to disk
%   7. FISSA             - neuropil correction via Python (see section below)
%   8. FISSA parsing     - load FISSA output, separate map vs. stim trials,
%                          apply neuropil scaling, save tifFileList
%   9. Stim alignment    - attach stimulus parameters to corrected traces
%                          via stimParam2ROI
%
% OUTPUT: tifFileList struct with motion- and neuropil-corrected fluorescence
%         traces (SCALEDfissaFroi, shape: ROI x frames) for each tif.
%

clearvars;close all;clc;

dataPath = uigetdir(".","Data path...");

if ~isfolder(dataPath)
    error('Data path not found')
end

try
    animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');
    if isempty(animal)
        error('')
    end
catch
    animal = char(inputdlg('Enter animal ID:','Animal Input',[1 35],{'AA0000'}));
end

%% 1. Get list of tif files for analysis and treatment info
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

%% 2. Split tifs into condition groups to be motion corrected separately
% DO NOT RUN IF RUNNING AGAIN FOR A DIFFERENT CELL TYPE
% Tifs are motion-corrected per group so that pre/post-treatment sessions
% are not corrupted by large inter-group motion differences.

%In case of REDO:
% clear Ycon rawCatImg options_nonrigid NoRMCorreParams moCorrImgNonRigid shifts2 template2
% load(fullfile(dataPath,[animal '_tifFileLegend.mat']),'tifFiles')

% --- Handle alternate-zoom acquisitions -----------------------------------
% Standard acquisitions use one zoom (scanZoomFactor); both 256x256 and the
% 256x128 10 Hz centered crop share it. A tif at a different zoom is a
% different (magnified) field of view that cannot share a motion-correction
% group or reuse the same ROIs. Policy: <=2 such tifs are stray
% misconfigurations -> dropped + warned. >2 -> prompt to keep them as a
% separate group (default No -> drop); headless/-batch runs take the default.
zoomPerTif = zeros(numel(tifFiles),1);
for i = 1:numel(tifFiles)
    [~,hZ] = readSCIMtif(fullfile(tifFiles(i).folder,tifFiles(i).name),'metaOnly');
    zoomPerTif(i) = hZ.hRoiManager.scanZoomFactor;
end
zoomMain = mode(zoomPerTif);
altIDX = zoomPerTif ~= zoomMain;
if any(altIDX)
    altNames = {tifFiles(altIDX).name};
    keepAlt = false;
    if sum(altIDX) > 2 && usejava('desktop') && ~batchStartupOptionUsed
        btn = questdlg(sprintf(['%d tifs are at a non-standard zoom (main zoom = %g):\n%s\n\n' ...
            'Keep them as a separate motion-correction group/set? (No = drop them)'], ...
            sum(altIDX), zoomMain, strjoin(altNames,', ')), ...
            'Alternate-zoom tifs','Yes','No','No');
        keepAlt = strcmp(btn,'Yes');
    end
    if keepAlt
        fprintf('Keeping %d alternate-zoom tif(s) as separate geometry group(s): %s\n', ...
            sum(altIDX), strjoin(altNames,', '));
    else
        warning('processAnimal2P:altZoomDropped',...
            'Dropping %d alternate-zoom tif(s) (main zoom = %g): %s', ...
            sum(altIDX), zoomMain, strjoin(altNames,', '));
        tifFiles = tifFiles(~altIDX);
        save(fullfile(dataPath,[animal '_tifFileLegend.mat']),'tifFiles')
    end
end
clear zoomPerTif zoomMain altIDX altNames keepAlt btn hZ i

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

% Auto-split each condition by scan GEOMETRY so every motion-correction group
% is dimensionally AND optically uniform. 256x256 and 256x128 tifs cannot share
% a group (concatenate_files stacks frames with cat(3,...) using the first
% tif's dims, so mixed pixel dims error); and same-size frames at different
% zoom are different fields of view that must not be co-registered or share
% ROIs. The split key is [H W zoom mFast mSlow]. The main-zoom full frame keeps
% the condition's base name; other groups get a <H>x<W> suffix (disambiguated
% by zoom if needed, e.g. 'postZX1' -> 'postZX1' + 'postZX1_128x256'). Uniform
% conditions are left untouched.
condNames = fieldnames(tifList);
for c = 1:numel(condNames)
    cond = condNames{c};
    grp  = tifList.(cond);
    g    = zeros(numel(grp),5);   % [H W zoom mFast mSlow]
    for i = 1:numel(grp)
        [~,hi] = readSCIMtif(fullfile(grp(i).folder,grp(i).name),'metaOnly');
        rm = hi.hRoiManager;
        g(i,:) = [hi.imHeight hi.imWidth rm.scanZoomFactor ...
                  rm.scanAngleMultiplierFast rm.scanAngleMultiplierSlow];
    end
    [uG,~,iu] = unique(g,'rows','stable');
    if size(uG,1) > 1
        % primary (keeps base name) = largest-area frame at the main zoom
        score = uG(:,1).*uG(:,2);
        score(uG(:,3)~=mode(g(:,3))) = -inf;
        [~,primary] = max(score);
        tifList = rmfield(tifList,cond);
        newNames = strings(size(uG,1),1);
        for s = 1:size(uG,1)
            if s == primary
                subName = cond;
            else
                subName = sprintf('%s_%dx%d',cond,uG(s,1),uG(s,2));
                if isfield(tifList,subName) || any(strcmp(newNames,subName))
                    subName = sprintf('%s_z%g',subName,uG(s,3));   % disambiguate by zoom
                end
                subName = matlab.lang.makeValidName(subName);
            end
            tifList.(subName) = grp(iu==s);
            newNames(s) = subName;
        end
        fprintf('Condition ''%s'' split by geometry -> %s\n', ...
            cond, strjoin(cellstr(newNames)',', '));
    end
end
clear condNames cond grp g hi rm uG iu score primary newNames s subName c i

save(fullfile(dataPath,[animal '_tifCondSplitLegend.mat']),'tifList')

%% 3. Motion correction via NoRMCorre
% DO NOT RUN IF RUNNING AGAIN FOR A DIFFERENT CELL TYPE

moCorN = fieldnames(tifList);
outputPath = fullfile(dataPath,'NoRMCorred');
if ~(exist(outputPath,'dir')==7)
    mkdir(outputPath)
end

% If tifs contain multiple channels (e.g. tdTomato + GCaMP), extract only
% the functional channel so that NoRMCorre and FISSA receive single-channel
% input. Channel 1 = redImg (e.g. tdTomato), channel 2 = greenImg (e.g. GCaMP).
% Original multi-channel tifs are moved to rawMergedTifs/.
[~,tmpHeader] = readSCIMtif(fullfile(tifList.(moCorN{1})(1).folder,...
    tifList.(moCorN{1})(1).name),'metaOnly');
if isfield(tmpHeader,'hChannels') && ...
        numel(tmpHeader.hChannels.channelSave)>1
    funcChan = str2double(inputdlg(...
        'Functional channel to extract (1 = redImg, 2 = greenImg):',...
        'Channel select',[1 60],{'2'}));
    splitTifs = cell(length(tifFiles),1);
    for tifM = 1:length(tifFiles)
        splitTifs(tifM) = splitTifChans(fullfile(tifFiles(tifM).folder,tifFiles(tifM).name),funcChan);
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
%motion correction
for k = 1:length(moCorN)
    %concatenate tifs
    [Ycon,~] = concatenate_files(tifList.(moCorN{k}));
    rawCattemp = single(Ycon); % convert to single precision
    rawCatImg.(moCorN{k}) = rawCattemp - min(rawCattemp(:)); % shift to non-negative for NoRMCorre
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

%% 3b. REDO / RESUME after motion correction
% RUN THIS BLOCK (uncomment %{ %}) INSTEAD OF SECTION 3 when:
%   - motion correction was already run and you are drawing ROIs for a new cell type
%   - picking up from a saved session (moCorrImgNonRigid is not in workspace)
% %{
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

%% 4. Draw ROIs on motion-corrected stack
% Run this section once per treatment condition, then run section 5 to save.
% Repeat for each condition before proceeding to section 6.

moCorSeqN = listdlg('PromptString','Select treatment condition (pre/post) for which to draw ROI',...
    'ListString',moCorN);
%outputs ROI to workspace and [animalID]_moCorrROI.mat into dataPath for use w/ FISSA
TIFcatROIgui(moCorrImgNonRigid.(moCorN{moCorSeqN}))

%% 5. Save ROIs — run after drawing ROIs for each treatment condition

nTifs = length(tifList.(moCorN{moCorSeqN}));
tifIDXinAllTifList = ismember({tifFiles.name}',{tifList.(moCorN{moCorSeqN}).name}');
% moCorr tif basenames for this condition (what writeMoCorTifs wrote to
% NoRMCorred/), in order. Lets the FISSA driver build per-group image lists
% and §9 map FISSA trials back to tifs by name.
moCorTifNames = strrep({tifList.(moCorN{moCorSeqN}).name}','.tif','_NoRMCorre.tif');

save([dataPath filesep ...
    regexp(dataPath,'[A-Z]{2}\d{4}','match','once') ...
    '_moCorrROI_' moCorN{moCorSeqN} '.mat'],...
    'moCorROI','moCorSeqN','nTifs','tifIDXinAllTifList','moCorTifNames')
disp([regexp(dataPath,'[A-Z]{2}\d{4}','match','once') ...
    '_moCorrROI_' moCorN{moCorSeqN} '.mat saved to animal directory'])
clear nTifs tifIDXinAllTifList moCorTifNames

%% 5b. Reuse 256x256 ROIs on any 256x128 (10 Hz spont) condition  [AUTO]
% Runs automatically with no manual lever: any motion-correction condition
% whose tifs are 256x128 is treated as a centered crop of the 256x256 field,
% and its ROIs are remapped from the matching 256x256 condition's moCorrROI
% file (drawn in sections 4-5). Only ROIs fully contained in the crop are
% kept; IDs are preserved. No-op when no 256x128 condition is present.
% remapROItoAcq ERRORS if the geometry is not the expected zoom-matched,
% centered crop (e.g. an erroneous zoom=2), so a misconfigured session fails
% loudly instead of mis-sampling cells.
% PREREQUISITE: draw + save (sections 4-5) the 256x256 condition(s) first.
moCorN = fieldnames(tifList);
isCrop = false(1,numel(moCorN));
for ci = 1:numel(moCorN)
    [~,hCi] = readSCIMtif(fullfile(tifList.(moCorN{ci})(1).folder,...
        tifList.(moCorN{ci})(1).name),'metaOnly');
    isCrop(ci) = hCi.imHeight==128 && hCi.imWidth==256;
end
src256 = moCorN(~isCrop);
for kc = find(isCrop)
    tgtCond = moCorN{kc};
    % resolve the 256x256 source condition: prefer one sharing a treatment
    % token with the crop condition, else the sole 256x256 condition
    match = src256(cellfun(@(s) contains(tgtCond,s) || contains(s,tgtCond), src256));
    if isscalar(match)
        srcCond = match{1};
    elseif isscalar(src256)
        srcCond = src256{1};
    else
        error('processAnimal2P:ambiguousSource',...
            ['Cannot resolve a unique 256x256 ROI source for crop condition ''%s''. '...
             'Candidates: %s. Rename conditions so the source shares a treatment token.'],...
            tgtCond, strjoin(src256,', '));
    end
    srcROIpath = fullfile(dataPath,[animal '_moCorrROI_' srcCond '.mat']);
    if exist(srcROIpath,'file')~=2
        error('processAnimal2P:noSourceROI',...
            'Draw 256x256 ROIs for condition ''%s'' (sections 4-5) before remapping ''%s''.',...
            srcCond, tgtCond);
    end
    remapROIfile(srcROIpath,...
        fullfile(tifList.(srcCond)(1).folder, tifList.(srcCond)(1).name),...
        fullfile(tifList.(tgtCond)(1).folder, tifList.(tgtCond)(1).name),...
        'outPath', fullfile(dataPath,[animal '_moCorrROI_' tgtCond '.mat']),...
        'nTifs', numel(tifList.(tgtCond)),...
        'tifIDXinAllTifList', ismember({tifFiles.name}',{tifList.(tgtCond).name}'),...
        'moCorTifNames', strrep({tifList.(tgtCond).name}','.tif','_NoRMCorre.tif'),...
        'moCorSeqN', kc);
end
clear isCrop ci hCi kc tgtCond srcCond match srcROIpath

%% 6. Match ROIs across treatment conditions
% Run only after sections 4-5 (and 5b for any 256x128 condition) are done.
% intersectROIfiles keeps only ROIs present in every condition so that the
% same set of cells is compared pre and post treatment. 256x128 (spont)
% conditions follow a separate analysis path and are EXCLUDED here so they
% do not reduce the 256x256 stim ROI set (5 Hz and 10 Hz are never pooled).
if exist('src256','var') && ~isempty(src256)
    intersectConds = src256;
else
    intersectConds = fieldnames(tifList);
end
intersectROIfiles(dataPath,animal,intersectConds,tifList,tifFiles)
clear intersectConds src256

%% 7. Extract raw fluorescence per ROI and add to tifList
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

%% 8. Neuropil correction via FISSA  [PYTHON STEP — run outside MATLAB]
%
%   FISSA separates each ROI's signal from contaminating neuropil by
%   decomposing the ROI trace and traces from surrounding neuropil rings.
%
%   Steps:
%     1. Edit the animalDataPath variable in FISSAviaMatlab_prePostTreatment.py
%     2. In a Python environment with FISSA installed, run:
%          python FISSAviaMatlab_prePostTreatment.py
%     3. FISSA reads ROI masks and the motion-corrected tifs from NoRMCorred/
%        and writes output to NoRMCorred/FISSAoutput/matlab.mat
%
%   Resume here (section 9) once matlab.mat exists.

%% 9. Parse FISSA output: separate map vs. stim trials, apply neuropil scaling

fissaDir = fullfile(dataPath,'NoRMCorred','FISSAoutput');

% fissaScaleFactor scales how aggressively neuropil is subtracted:
%   corrected = ROI - scaleFactor * neuropil   (0.8 = common conservative default)
fissaScaleFactor = str2double(inputdlg('ENTER FACTOR BY WHICH TO SCALE FISSA SUBTRACTION (eg. 0.8): ',...
    'FISSA SCALING FACTOR',[1 80],{'0.8'}));

if exist(fullfile(fissaDir,'groups.json'),'file')==2
% ===== GROUPED PATH: per-ROI-count FISSA outputs (mixed 256x256 / 256x128) =====
% FISSA was run once per ROI-count group (different frame sizes cannot share a
% run; see FISSAviaMatlab_prePostTreatment.py). Build tifFileList.stim/.map
% from tifList (map split by treatment), then attach each tif's neuropil traces
% from the matching group, by name, via mergeFISSAgroups. Each tif's trace row
% count follows its own group's ROI count (e.g. 256x256 stim vs 256x128 spont).
C = struct2cell(tifList);
allcond = vertcat(C{:});
isMap = contains({allcond.treatment}','map');
tifFileList = struct();
tifFileList.stim = allcond(~isMap);
if any(isMap)
    tifFileList.map = allcond(isMap);
end
tifFileList = mergeFISSAgroups(tifFileList, fissaDir, fissaScaleFactor);
FISSAoutput = fileread(fullfile(fissaDir,'groups.json'));   % manifest text, for provenance

else
% ===== LEGACY single-output path (unchanged: no 256x128, one ROI count) =====
%FISSA was given moCorr tiffs and ROIs from concatenated moCorr Tifs
FISSAoutput = load(fullfile(fissaDir,'matlab.mat'));

%in 'result', for a given cell and trial there is a n x numTraceFrames double;
%row 1 is ROI trace, rows 2->n are traces from neuropil regions around ROI

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

tifFileList = FISSAoutput2tifFileList(FISSAoutput,tifFileList,fissaScaleFactor);
end

%if updating existing tifFileList w/ FISSAoutput:
% [fName,fPath] = uigetfile('*tifFileList.mat','Locate [ANIMAL]_tifFileList.mat...');
% load(fullfile(fPath,fName))

save(fullfile(dataPath,[animal '_tifFileList.mat']),...
        'dataPath','FISSAoutput','tifFileList','fissaScaleFactor','-v7.3')

%% COMPLETE. YOU ARE NOW LEFT WITH MOTION AND NEUROPIL CORRECTED FLUORESCENCE TRACES FOR EACH ROI FOR EACH TIF
%
%   tifFileList.stim(n).SCALEDfissaFroi  ->  nROI x nFrames for the nth stimulus tif
%   tifFileList.map(n).SCALEDfissaFroi   ->  nROI x nFrames for the nth FRA/BF mapping tif

%% 10. Align stimulus parameters to corrected traces
% Requires _Pulses.mat files co-located with tifs.
[pulseLegend2P,stimGroupIDX,ROIoutputTables] = stimParam2ROI(dataPath);