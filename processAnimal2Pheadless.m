function out = processAnimal2Pheadless(dataPath,varargin)
% processAnimal2Pheadless  Run the whole per-animal pipeline unattended.
%
%   out = processAnimal2Pheadless(dataPath)
%   out = processAnimal2Pheadless(dataPath,'treatmentName','ZX1','preTifs',1:8)
%   out = processAnimal2Pheadless(dataPath,'stages',[4 11])
%   out = processAnimal2Pheadless(dataPath,cfgStruct)
%
%   A second entry point to the pipeline processAnimal2P runs interactively.
%   It is NOT a replacement: processAnimal2P stays the path for hand-drawn
%   ROIs and for anything you want to watch. This runs the same stages with
%   every dialog replaced by a config field (see headlessConfig) and
%   TIFcatROIgui replaced by Cellpose segmentation (see cellposeROIset), so a
%   whole animal -- or a cohort, in a loop -- can go through in one call.
%
%   Stages, numbered to match processAnimal2P's sections:
%      1  tif inventory + treatment/FRAmap labels  -> _tifFileLegend.mat
%      2  condition split (treatment x geometry)   -> _tifCondSplitLegend.mat
%      3  NoRMCorre motion correction              -> NoRMCorred/*.tif
%      4  ROI detection (Cellpose) + save + crop remap -> _moCorrROI_<cond>.mat
%      6  cross-condition ROI matching             (intersectROIfiles)
%      7  raw F extraction                         -> _moCorr_Tifs_Params.mat
%      8  FISSA neuropil correction (Python)       -> FISSAoutput/
%      9  FISSA parsing + neuropil scaling         -> _tifFileList.mat
%     10  stimulus alignment                       -> _anmlROI_<Fam>*_raw.mat
%     11  per-family dF/F + peak responses         -> _anmlROI_<Fam>*.mat
%         (including FRA, which has no _raw stage and is driven straight
%          from the tif inventory)
%   (5 is folded into 4, as saving is not a separate step here.)
%
%   Each stage is skipped when its artefact already exists unless
%   cfg.overwrite is set, so re-running to re-tune segmentation does not
%   repeat motion correction. Use 'stages' to run a subset.
%
%   Inputs
%     dataPath - animal data folder.
%     varargin - name-value pairs, or a single struct; both go to
%                headlessConfig, which documents every field and its default.
%
%   Output (struct)
%     .cfg      the resolved config
%     .stages   per-stage struct array: .n .name .ran .skipped .reason
%               .seconds .error
%     .roi      cellposeROIset's output, when stage 4 ran
%     .stim     processAnimalStimFamilies' output, when stage 11 ran
%     .files    key artefact paths
%
%   A stage that errors stops the run and is reported in .stages; earlier
%   stages have already written their artefacts, so a fixed re-run resumes
%   from where it stopped.
%
%   Example -- a cohort, resuming past motion correction:
%     for a = ["TO0006","TO0007"]
%         processAnimal2Pheadless(fullfile(root,a),'stages',[4 11], ...
%             'treatmentName','ZX1','preTifs','_0003[4-9]_');
%     end
%
%   See also headlessConfig, cellposeROIset, processAnimal2P,
%   stimParam2ROI, processAnimalStimFamilies, runFRA

cfg = headlessConfig(dataPath,varargin{:});

animal   = cfg.animal;
dataPath = cfg.dataPath;
V        = cfg.verbose;

out = struct('cfg',cfg,'stages',struct('n',{},'name',{},'ran',{},...
    'skipped',{},'reason',{},'seconds',{},'error',{}),...
    'roi',[],'stim',[],'files',struct());

F = struct();
F.legend     = fullfile(dataPath,[animal '_tifFileLegend.mat']);
F.condSplit  = fullfile(dataPath,[animal '_tifCondSplitLegend.mat']);
F.moCorrDir  = fullfile(dataPath,'NoRMCorred');
F.ncParams   = fullfile(F.moCorrDir,[animal '_NoRMCorreParams.mat']);
F.moCorrTifs = fullfile(dataPath,[animal '_moCorr_Tifs_Params.mat']);
F.fissaDir   = fullfile(F.moCorrDir,'FISSAoutput');
F.tifFileList= fullfile(dataPath,[animal '_tifFileList.mat']);
out.files = F;

%workspace carried between stages
W = struct('tifFiles',[],'tifList',[],'moCorrImgNonRigid',[],'rawCatImg',[]);

%stage numbers and names; the bodies are dispatched below against the CURRENT
%workspace, not captured in a handle -- each stage feeds the next
stageDefs = { 1,'tif inventory';
              2,'condition split';
              3,'motion correction';
              4,'ROI detection';
              6,'ROI matching';
              7,'raw F extraction';
              8,'FISSA';
              9,'FISSA parsing';
             10,'stim alignment';
             11,'per-stim processing' };

for s = 1:size(stageDefs,1)
    n    = stageDefs{s,1};
    name = stageDefs{s,2};
    if ~ismember(n,cfg.stages) && ~(n==4 && ismember(5,cfg.stages))
        continue
    end

    rec = struct('n',n,'name',name,'ran',false,'skipped',false,...
        'reason','','seconds',0,'error','');
    tS = tic;
    if V
        fprintf('\n===== stage %d: %s =====\n',n,name);
    end
    try
        %stages need each other's outputs, so anything not in W yet is
        %loaded from disk -- that is what makes a partial run resumable
        W = ensureWorkspace(n,cfg,F,W);
        [done,reason] = alreadyDone(n,cfg,F);
        if done && ~cfg.overwrite
            rec.skipped = true;
            rec.reason  = reason;
            if V, fprintf('skipped: %s\n',reason); end
        else
            r = runStage(n,cfg,F,W);
            if isstruct(r)
                fn = fieldnames(r);
                for k = 1:numel(fn)
                    if isfield(W,fn{k}) || any(strcmp(fn{k},{'roi','stim'}))
                        W.(fn{k}) = r.(fn{k});
                    end
                end
                if isfield(r,'roi'),  out.roi  = r.roi;  end
                if isfield(r,'stim'), out.stim = r.stim; end
            end
            rec.ran = true;
        end
    catch ME
        rec.error   = ME.message;
        rec.seconds = toc(tS);
        out.stages(end+1) = rec;
        fprintf(2,'stage %d (%s) FAILED: %s\n',n,name,ME.message);
        rethrow(ME);
    end
    rec.seconds = toc(tS);
    out.stages(end+1) = rec;
    if V && rec.ran
        fprintf('stage %d done in %.1f s\n',n,rec.seconds);
    end
end

if V
    fprintf('\n===== %s complete =====\n',animal);
    for k = 1:numel(out.stages)
        st = out.stages(k);
        if st.skipped
            fprintf('  %2d %-22s skipped (%s)\n',st.n,st.name,st.reason);
        else
            fprintf('  %2d %-22s %.1f s\n',st.n,st.name,st.seconds);
        end
    end
end
end

% ========================================================================
function r = runStage(n,cfg,F,W)
switch n
    case 1,  r = stage1(cfg,F);
    case 2,  r = stage2(cfg,F,W);
    case 3,  r = stage3(cfg,F,W);
    case 4,  r = stage4(cfg,F,W);
    case 6,  r = stage6(cfg,F,W);
    case 7,  r = stage7(cfg,F,W);
    case 8,  r = stage8(cfg,F);
    case 9,  r = stage9(cfg,F,W);
    case 10, r = stage10(cfg,F);
    case 11, r = stage11(cfg,F);
    otherwise
        error('processAnimal2Pheadless:badStage','No stage %d.',n);
end
end

% ========================================================================
function [done,reason] = alreadyDone(n,cfg,F)
done = false; reason = '';
switch n
    case 1, done = isfile(F.legend);     reason = 'tifFileLegend exists';
    case 2, done = isfile(F.condSplit);  reason = 'tifCondSplitLegend exists';
    case 3
        done = isfolder(F.moCorrDir) && ~isempty(dir(fullfile(F.moCorrDir,'*_NoRMCorre.tif')));
        reason = 'NoRMCorred tifs exist';
    case 4
        done = ~isempty(dir(fullfile(cfg.dataPath,[cfg.animal '_moCorrROI_*.mat'])));
        reason = 'moCorrROI files exist';
    case 7, done = isfile(F.moCorrTifs);  reason = 'moCorr_Tifs_Params exists';
    case 8
        done = isfile(fullfile(F.fissaDir,'matlab.mat')) || ...
               isfile(fullfile(F.fissaDir,'groups.json'));
        reason = 'FISSA output exists';
    case 9, done = isfile(F.tifFileList); reason = 'tifFileList exists';
    case 10
        done = ~isempty(dir(fullfile(cfg.dataPath,[cfg.animal '_anmlROI_*_raw.mat'])));
        reason = 'raw stim tables exist';
end
%6 and 11 are cheap and idempotent enough to just re-run
end

% ========================================================================
function W = ensureWorkspace(n,cfg,F,W)
%load whatever this stage needs and does not have yet
if n >= 2 && isempty(W.tifFiles) && isfile(F.legend)
    S = load(F.legend,'tifFiles');
    W.tifFiles = repointFolders(S.tifFiles,cfg.dataPath,'tifFileLegend');
end
if n >= 3 && isempty(W.tifList) && isfile(F.condSplit)
    S = load(F.condSplit,'tifList');
    cn = fieldnames(S.tifList);
    for c = 1:numel(cn)
        S.tifList.(cn{c}) = repointFolders(S.tifList.(cn{c}),cfg.dataPath,...
            sprintf('tifCondSplitLegend.%s',cn{c}));
    end
    W.tifList = S.tifList;
end
%stage 7 is the only one that needs the image stacks in memory
if n == 7 && isempty(W.moCorrImgNonRigid)
    if cfg.verbose
        fprintf('loading motion-corrected and raw stacks per condition...\n');
    end
    condN = fieldnames(W.tifList);
    for k = 1:numel(condN)
        W.moCorrImgNonRigid.(condN{k}) = ...
            loadNoRMCorrNonRigidImgViaTifs(cfg.dataPath,W.tifList.(condN{k}));
        try
            [Ycon,~] = concatenate_files(W.tifList.(condN{k}));
        catch
            [W.tifList.(condN{k}).folder] = deal(cfg.dataPath);
            [Ycon,~] = concatenate_files(W.tifList.(condN{k}));
        end
        raw = single(Ycon);
        W.rawCatImg.(condN{k}) = raw - min(raw(:));
        clear Ycon raw
    end
end
end

% ========================================================================
function f = repointFolders(f,dataPath,what)
%A legend built on another machine carries that machine's paths -- these
%folders routinely read C:\Users\...\OneDrive\Data\<animal>. Every consumer
%joins .folder with .name, so a stale path fails far downstream (in
%moCorRawF2tifList, say) with a Windows path in the message. Repoint to the
%folder the legend itself was found in, which is where the tifs are. A folder
%that does exist is left alone -- that is how a legitimately split layout
%(TO0002's BPN subfolder) keeps working.
if isempty(f) || ~isfield(f,'folder')
    return
end
stale = ~arrayfun(@(x) isfolder(x.folder), f);
if any(stale)
    warning('processAnimal2Pheadless:stalePaths',...
        '%s: %d of %d tif folder(s) do not exist here (e.g. %s); repointing to %s',...
        what,nnz(stale),numel(f),f(find(stale,1)).folder,dataPath);
    [f(stale).folder] = deal(dataPath);
end
end

% ========================================================================
function r = stage1(cfg,F)
%tif inventory + treatment and FRAmap labels
tifFiles = dir(fullfile(cfg.dataPath,cfg.tifPattern));
tifFiles = tifFiles(~[tifFiles.isdir]);
if isempty(tifFiles)
    error('processAnimal2Pheadless:noTifs',...
        'No tifs matching %s in %s',cfg.tifPattern,cfg.dataPath);
end

treatment = cell(numel(tifFiles),1);
if isempty(cfg.treatmentName)
    treatment(:) = {'none'};
else
    treatment(cfg.preIDX)  = {['pre'  cfg.treatmentName]};
    treatment(~cfg.preIDX) = {['post' cfg.treatmentName]};
end
if any(cfg.mapIDX)
    treatment(cfg.mapIDX) = cellfun(@(c) [c ' FRAmap'],treatment(cfg.mapIDX),'uni',0);
end
[tifFiles.treatment] = treatment{:};

save(F.legend,'tifFiles')
if cfg.verbose
    u = unique(treatment);
    for k = 1:numel(u)
        fprintf('  %-24s %d tifs\n',u{k},nnz(strcmp(treatment,u{k})));
    end
end
r = struct('tifFiles',tifFiles);
end

% ========================================================================
function r = stage2(cfg,F,W)
%condition split: alternate zoom policy, treatment filters, geometry split
tifFiles = W.tifFiles;

%--- alternate-zoom tifs ---
zoomPerTif = zeros(numel(tifFiles),1);
for i = 1:numel(tifFiles)
    [~,hZ] = readSCIMtif(fullfile(tifFiles(i).folder,tifFiles(i).name),'metaOnly');
    zoomPerTif(i) = hZ.hRoiManager.scanZoomFactor;
end
zoomMain = mode(zoomPerTif);
altIDX = zoomPerTif ~= zoomMain;
if any(altIDX)
    altNames = {tifFiles(altIDX).name};
    switch lower(cfg.altZoomPolicy)
        case 'error'
            error('processAnimal2Pheadless:altZoom',...
                '%d tif(s) at a non-standard zoom (main %g): %s',...
                sum(altIDX),zoomMain,strjoin(altNames,', '));
        case 'keep'
            fprintf('Keeping %d alternate-zoom tif(s) as their own geometry group: %s\n',...
                sum(altIDX),strjoin(altNames,', '));
        otherwise
            warning('processAnimal2Pheadless:altZoomDropped',...
                'Dropping %d alternate-zoom tif(s) (main zoom = %g): %s',...
                sum(altIDX),zoomMain,strjoin(altNames,', '));
            tifFiles = tifFiles(~altIDX);
            save(F.legend,'tifFiles')
    end
end

%--- treatment filters ---
filters = cfg.condFilters;
if isempty(filters)
    tr = unique(strrep({tifFiles.treatment}',' FRAmap',''));
    if isscalar(tr) && strcmp(tr{1},'none')
        filters = {};
    else
        filters = tr(:)';
    end
end
tifList = struct();
if isempty(filters)
    tifList.all = tifFiles;
else
    filters = cellstr(filters);
    locLogical = false(numel(tifFiles),1);
    for k = 1:numel(filters)
        hit = contains({tifFiles.treatment}',filters{k});
        tifList.(matlab.lang.makeValidName(filters{k})) = tifFiles(hit);
        locLogical(hit) = true;
    end
    if nnz(locLogical) < numel(tifFiles)
        tifList.remaining = tifFiles(~locLogical);
    end
end

%--- geometry sub-split (same rule as processAnimal2P §2) ---
condNames = fieldnames(tifList);
for c = 1:numel(condNames)
    cond = condNames{c};
    grp  = tifList.(cond);
    if isempty(grp), tifList = rmfield(tifList,cond); continue, end
    g = zeros(numel(grp),5);
    for i = 1:numel(grp)
        [~,hi] = readSCIMtif(fullfile(grp(i).folder,grp(i).name),'metaOnly');
        rm = hi.hRoiManager;
        g(i,:) = [hi.imHeight hi.imWidth rm.scanZoomFactor ...
                  rm.scanAngleMultiplierFast rm.scanAngleMultiplierSlow];
    end
    [uG,~,iu] = unique(g,'rows','stable');
    if size(uG,1) > 1
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
                    subName = sprintf('%s_z%g',subName,uG(s,3));
                end
                subName = matlab.lang.makeValidName(subName);
            end
            tifList.(subName) = grp(iu==s);
            newNames(s) = subName;
        end
        fprintf('Condition ''%s'' split by geometry -> %s\n',...
            cond,strjoin(cellstr(newNames)',', '));
    end
end

save(F.condSplit,'tifList')
if cfg.verbose
    cn = fieldnames(tifList);
    for k = 1:numel(cn)
        fprintf('  %-24s %d tifs\n',cn{k},numel(tifList.(cn{k})));
    end
end
r = struct('tifFiles',tifFiles,'tifList',tifList);
end

% ========================================================================
function r = stage3(cfg,F,W)
%NoRMCorre, per condition group
tifFiles = W.tifFiles;
tifList  = W.tifList;
moCorN   = fieldnames(tifList);
if ~isfolder(F.moCorrDir), mkdir(F.moCorrDir); end

%--- split multi-channel tifs down to the functional channel ---
[~,tmpHeader] = readSCIMtif(fullfile(tifList.(moCorN{1})(1).folder,...
    tifList.(moCorN{1})(1).name),'metaOnly');
if isfield(tmpHeader,'hChannels') && numel(tmpHeader.hChannels.channelSave) > 1
    if cfg.verbose
        fprintf('extracting functional channel %d from %d multi-channel tifs\n',...
            cfg.funcChan,numel(tifFiles));
    end
    splitTifs = cell(numel(tifFiles),1);
    for m = 1:numel(tifFiles)
        splitTifs(m) = splitTifChans(fullfile(tifFiles(m).folder,tifFiles(m).name),cfg.funcChan);
    end
    rawDir = fullfile(tifFiles(1).folder,'rawMergedTifs');
    if ~isfolder(rawDir), mkdir(rawDir); end
    for m = 1:numel(tifFiles)
        movefile(fullfile(tifFiles(m).folder,tifFiles(m).name),...
            fullfile(rawDir,tifFiles(m).name))
    end
    for m = 1:numel(tifFiles)
        movefile(splitTifs{m},fullfile(tifFiles(m).folder,tifFiles(m).name))
    end
end

moCorrImgNonRigid = struct();
rawCatImg = struct();
NoRMCorreParams = struct();
for k = 1:numel(moCorN)
    if cfg.verbose
        fprintf('  %s: concatenating %d tifs...\n',moCorN{k},numel(tifList.(moCorN{k})));
    end
    [Ycon,~] = concatenate_files(tifList.(moCorN{k}));
    raw = single(Ycon);
    rawCatImg.(moCorN{k}) = raw - min(raw(:));
    clear Ycon raw

    nc = cfg.normcorre;
    opts_nr = NoRMCorreSetParms('d1',size(rawCatImg.(moCorN{k}),1),...
        'd2',size(rawCatImg.(moCorN{k}),2),...
        'grid_size',nc.grid_size,'mot_uf',nc.mot_uf,'bin_width',nc.bin_width,...
        'max_shift',nc.max_shift,'max_dev',nc.max_dev,'us_fac',nc.us_fac,...
        'init_batch',nc.init_batch);

    if cfg.verbose, fprintf('  %s: NoRMCorre...\n',moCorN{k}); end
    [moCorrImgNonRigid.(moCorN{k}),NoRMCorreParams.(moCorN{k}).shifts,~,...
        NoRMCorreParams.(moCorN{k}).options_nonrigid] = ...
        normcorre_batch(rawCatImg.(moCorN{k}),opts_nr);

    writeMoCorTifs(tifList.(moCorN{k}),moCorrImgNonRigid.(moCorN{k}),F.moCorrDir)
end
save(F.ncParams,'NoRMCorreParams')

r = struct('moCorrImgNonRigid',moCorrImgNonRigid,'rawCatImg',rawCatImg);
end

% ========================================================================
function r = stage4(cfg,F,W)
%ROI detection via Cellpose, then the 256x128 crop remap (§5b)
tifList  = W.tifList;
tifFiles = W.tifFiles;
condN    = fieldnames(tifList);

%which conditions are 256x128 centred crops? those reuse a 256x256 set
isCrop = false(1,numel(condN));
for c = 1:numel(condN)
    [~,hC] = readSCIMtif(fullfile(tifList.(condN{c})(1).folder,...
        tifList.(condN{c})(1).name),'metaOnly');
    isCrop(c) = hC.imHeight==128 && hC.imWidth==256;
end
fullConds = condN(~isCrop);
if isempty(fullConds)
    error('processAnimal2Pheadless:noFullFrameCond',...
        'Every condition is a 256x128 crop; there is no full-frame condition to segment.');
end

roiArgs = cfg.roi;
roiNV = struct2nv(roiArgs);
roiOut = cellposeROIset(cfg.dataPath,cfg.animal,tifList,tifFiles,...
    'conds',fullConds,roiNV{:},'verbose',cfg.verbose);

%--- §5b: derive each crop condition from its 256x256 source ---
for c = find(isCrop)
    tgtCond = condN{c};
    match = fullConds(cellfun(@(s) contains(tgtCond,s) || contains(s,tgtCond),fullConds));
    if isscalar(match)
        srcCond = match{1};
    elseif isscalar(fullConds)
        srcCond = fullConds{1};
    else
        error('processAnimal2Pheadless:ambiguousSource',...
            ['Cannot resolve a unique 256x256 ROI source for crop condition ''%s''. '...
             'Candidates: %s.'],tgtCond,strjoin(fullConds',', '));
    end
    srcROIpath = fullfile(cfg.dataPath,[cfg.animal '_moCorrROI_' srcCond '.mat']);
    remapROIfile(srcROIpath,...
        fullfile(tifList.(srcCond)(1).folder,tifList.(srcCond)(1).name),...
        fullfile(tifList.(tgtCond)(1).folder,tifList.(tgtCond)(1).name),...
        'outPath',fullfile(cfg.dataPath,[cfg.animal '_moCorrROI_' tgtCond '.mat']),...
        'nTifs',numel(tifList.(tgtCond)),...
        'tifIDXinAllTifList',ismember({tifFiles.name}',{tifList.(tgtCond).name}'),...
        'moCorTifNames',strrep({tifList.(tgtCond).name}','.tif','_NoRMCorre.tif'),...
        'moCorSeqN',c);
end

r = struct('roi',roiOut);
end

% ========================================================================
function r = stage6(cfg,F,W)
%keep only ROIs present in every full-frame condition. With a shared
%Cellpose set this is normally a no-op, which is the point: it verifies the
%invariant instead of being what establishes it.
tifList = W.tifList;
condN   = fieldnames(tifList);
isCrop  = false(1,numel(condN));
for c = 1:numel(condN)
    [~,hC] = readSCIMtif(fullfile(tifList.(condN{c})(1).folder,...
        tifList.(condN{c})(1).name),'metaOnly');
    isCrop(c) = hC.imHeight==128 && hC.imWidth==256;
end
intersectConds = condN(~isCrop);

%legacy ROI files (drawn before the grouped FISSA driver existed) carry only
%the ROIs; FISSA needs moCorTifNames and errors without it
ensureROIfileMeta(cfg.dataPath,cfg.animal,tifList,W.tifFiles);

before = zeros(1,numel(intersectConds));
for c = 1:numel(intersectConds)
    S = load(fullfile(cfg.dataPath,[cfg.animal '_moCorrROI_' intersectConds{c} '.mat']),'moCorROI');
    before(c) = numel(S.moCorROI);
end
intersectROIfiles(cfg.dataPath,cfg.animal,intersectConds,tifList,W.tifFiles)
for c = 1:numel(intersectConds)
    S = load(fullfile(cfg.dataPath,[cfg.animal '_moCorrROI_' intersectConds{c} '.mat']),'moCorROI');
    if cfg.verbose
        fprintf('  %-24s %d -> %d ROIs\n',intersectConds{c},before(c),numel(S.moCorROI));
    end
end
r = struct();
end

% ========================================================================
function r = stage7(cfg,F,W)
%raw F per ROI per tif, from both the corrected and uncorrected stacks
tifList = W.tifList;
condN   = fieldnames(tifList);
for c = 1:numel(condN)
    S = load(fullfile(cfg.dataPath,[cfg.animal '_moCorrROI_' condN{c} '.mat']),'moCorROI');
    if cfg.verbose
        fprintf('  %s: %d ROIs x %d tifs\n',condN{c},numel(S.moCorROI),numel(tifList.(condN{c})));
    end
    tifList.(condN{c}) = moCorRawF2tifList(tifList.(condN{c}),...
        W.moCorrImgNonRigid.(condN{c}),S.moCorROI,W.rawCatImg.(condN{c}));
end

allTifFiles = W.tifFiles; %#ok<NASGU>
NoRMCorreParams = [];
if isfile(F.ncParams)
    Sp = load(F.ncParams,'NoRMCorreParams');
    NoRMCorreParams = Sp.NoRMCorreParams; %#ok<NASGU>
end
save(F.moCorrTifs,'tifList','NoRMCorreParams','allTifFiles','-v7.3')

r = struct('tifList',tifList,'moCorrImgNonRigid',[],'rawCatImg',[]);
end

% ========================================================================
function r = stage8(cfg,F)
%the Python neuropil step
if isempty(cfg.fissaCmd)
    error('processAnimal2Pheadless:noFissaCmd',...
        ['FISSA output is missing and fissaCmd is empty. Run '...
         'FISSAviaMatlab_prePostTreatment.py on %s, or set fissaCmd.'],cfg.dataPath);
end
cmd = sprintf(cfg.fissaCmd,cfg.dataPath);
%MATLAB points LD_LIBRARY_PATH at its own libs, which breaks system python
cmd = ['env -u LD_LIBRARY_PATH -u LD_PRELOAD ' cmd];
if cfg.verbose, fprintf('%s\n',cmd); end
[st,outTxt] = system(cmd);
if st ~= 0
    error('processAnimal2Pheadless:fissaFailed',...
        'FISSA failed (exit %d).\ncommand: %s\noutput:\n%s',st,cmd,outTxt);
end
if ~isfile(fullfile(F.fissaDir,'matlab.mat')) && ...
        ~isfile(fullfile(F.fissaDir,'groups.json'))
    error('processAnimal2Pheadless:noFissaOutput',...
        'FISSA reported success but wrote no output to %s.\n%s',F.fissaDir,outTxt);
end
r = struct();
end

% ========================================================================
function r = stage9(cfg,F,W)
%parse FISSA output into tifFileList, applying the neuropil scale factor
%must be the tifList stage 7 wrote: FISSAoutput2tifFileList reads
%moCorRawFroi off each entry, and the condition-split legend's copy has no
%fluorescence in it
if isfile(F.moCorrTifs)
    S = load(F.moCorrTifs,'tifList');
    tifList = S.tifList;
else
    tifList = W.tifList;
end
if isempty(tifList) || ~isfield(tifList.(subsref(fieldnames(tifList),...
        substruct('{}',{1}))),'moCorRawFroi')
    error('processAnimal2Pheadless:noRawF',...
        ['Stage 9 needs the raw F extracted in stage 7 (%s). Run stage 7 '...
         'first.'],F.moCorrTifs);
end
fissaScaleFactor = cfg.fissaScaleFactor;

if isfile(fullfile(F.fissaDir,'groups.json'))
    %grouped path: one FISSA run per ROI count
    C = struct2cell(tifList);
    allcond = vertcat(C{:});
    isMap = contains({allcond.treatment}','map');
    tifFileList = struct();
    tifFileList.stim = allcond(~isMap);
    if any(isMap)
        tifFileList.map = allcond(isMap);
    end
    tifFileList = mergeFISSAgroups(tifFileList,F.fissaDir,fissaScaleFactor);
    FISSAoutput = fileread(fullfile(F.fissaDir,'groups.json'));
else
    FISSAoutput = load(fullfile(F.fissaDir,'matlab.mat'));
    C = struct2cell(tifList);
    allcond = vertcat(C{:});
    FRAmapIDX = contains({allcond.treatment}','map');

    tifFileList = struct();
    if any(FRAmapIDX)
        tifFileList.map  = allcond(FRAmapIDX);
        tifFileList.stim = allcond(~FRAmapIDX);

        fID = fieldnames(FISSAoutput);
        trials.BF   = strcat({'trial'},cellstr(string(0:sum(FRAmapIDX)-1))');
        trials.stim = strcat({'trial'},cellstr(string(0:sum(~FRAmapIDX)-1))');
        for nField = 1:numel(fID)
            cells = fieldnames(FISSAoutput.(fID{nField}));
            for nCell = 1:numel(cells)
                tmp = struct2cell(FISSAoutput.(fID{nField}).(['cell' num2str(nCell-1)]));
                FISSAout.map.(fID{nField}).(['cell' num2str(nCell-1)])  = ...
                    cell2struct(tmp(FRAmapIDX),trials.BF);
                FISSAout.stim.(fID{nField}).(['cell' num2str(nCell-1)]) = ...
                    cell2struct(tmp(~FRAmapIDX),trials.stim);
            end
        end
        clear FISSAoutput
        FISSAoutput.map  = FISSAout.map;
        FISSAoutput.stim = FISSAout.stim;
    else
        tifFileList.stim = allcond;
        tmp = FISSAoutput;
        clear FISSAoutput
        FISSAoutput.stim = tmp;
    end
    tifFileList = FISSAoutput2tifFileList(FISSAoutput,tifFileList,fissaScaleFactor);
end

dataPath = cfg.dataPath; %#ok<NASGU>
save(F.tifFileList,'dataPath','FISSAoutput','tifFileList','fissaScaleFactor','-v7.3')
if cfg.verbose
    fprintf('  %d stim tifs',numel(tifFileList.stim));
    if isfield(tifFileList,'map'), fprintf(', %d map tifs',numel(tifFileList.map)); end
    fprintf('; traces are %d ROIs x %d frames\n',...
        size(tifFileList.stim(1).SCALEDfissaFroi,1),...
        size(tifFileList.stim(1).SCALEDfissaFroi,2));
end
r = struct();
end

% ========================================================================
function r = stage10(cfg,~)
%FRA is not handled here: it has no _raw stage, and processAnimalStimFamilies
%(stage 11) drives it straight from the tif inventory along with every other
%family.
[~,~,~] = stimParam2ROI(cfg.dataPath,'excludeNeg',cfg.excludeNeg);
r = struct();
end

% ========================================================================
function r = stage11(cfg,~)
if ~cfg.runStimFamilies
    r = struct(); return
end
args = {'showPlots',false,'scriptVars',cfg.stimScriptVars,'verbose',cfg.verbose};
if ~cfg.runFRAmap
    fams = stimGroupSpec();
    args = [args,{'families',fams(~strcmp(fams,'FRA'))}];
end
stim = processAnimalStimFamilies(cfg.dataPath,args{:});
r = struct('stim',stim);
end

% ========================================================================
function nv = struct2nv(s)
f = fieldnames(s);
nv = cell(1,2*numel(f));
for k = 1:numel(f)
    nv{2*k-1} = f{k};
    nv{2*k}   = s.(f{k});
end
end
