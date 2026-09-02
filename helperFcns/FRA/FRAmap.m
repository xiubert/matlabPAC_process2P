function [FRAoutput,PTfreq,PTdBampl] = FRAmap(tifFileList,...
    pkPTsigSD,nFramesPostPulse,varargin)
% TRFmap: output TRF map for an animal.
%   [FRAoutput,PTfreq,PTdBampl] = TRFmap(
%           tifFileList, --> list of BF map tif files w/ F signal
%           pkPTsigSD, --> SD threshold for a significant PT response
%           nFramesPostPulse, --> frames after pulse ends for peak search
%           (eg if pulse is 400 ms and FR is 5, total frames post pulse
%           onset will be 2 + nFramesPostPulse
%           varargin --> optionally define Fsource string
%
%   FRAoutput: structure with FRAmap, BF, and dPrime
%   PTfreq: list of pure-tone frequencies in order of presentation
%   PTdBampl: list of pure-tone dB amplitudes in order of presentation
%
%   See also plotFRAmap.

% Get pulse params for each tif
p = inputParser;
addRequired(p,'tifFileList',@isstruct)
addRequired(p,'pkPTsigSD',@isscalar)
addRequired(p,'nFramesPostPulse',@isscalar)
addOptional(p,'FsourceString','SCALEDfissaFroi',@ischar);

parse(p,tifFileList,pkPTsigSD,nFramesPostPulse,varargin{:});

tifFileList = p.Results.tifFileList;
pkPTsigSD = p.Results.pkPTsigSD;
nFramesPostPulse = p.Results.nFramesPostPulse;
FsourceString = p.Results.FsourceString;

tifDir = tifFileList.map(1).folder;
[PTfreq,PTdBampl,PTonsetInPulse,msPTpulseLen,paramS] = deal(cell(length(tifFileList.map),1));
for nTif = 1:length(tifFileList.map)
    pulseFile = strrep(tifFileList.map(nTif).name,'.tif','_Pulses.mat');
    try
        S = load(fullfile(tifDir,pulseFile));
    catch
        %prompting is only possible in an interactive session; under -batch
        %the dialog throws an opaque error, so fail with the missing path
        if ~usejava('desktop')
            error(['Cannot find pulse file:\n  %s\n' ...
                'tifFileList.map(1).folder may be stale (e.g. a path from ' ...
                'another machine). Set it to the folder holding the map ' ...
                'tifs and their _Pulses.mat files.'],fullfile(tifDir,pulseFile));
        end
        disp('Can''t find pulses associated with .tifs, locate dir...')
        tifDir = uigetdir(pwd,...
            'Locate directory containing TRF map tifs...');
        S = load(fullfile(tifDir,pulseFile));
    end
    
    if isfield(S.params, 'pulseFrameNo')
        paramS{nTif} = rmfield(S.params, 'pulseFrameNo');
    else
        paramS{nTif} = S.params;
    end
    if length(S.pulse)>S.params.totalPulses
        error('Too many pulses in pulse file, may have been appended in error')
    elseif length(S.pulse)<S.params.totalPulses
        disp(['Warning: pulse number lower than noted in pulse params (tif file no. ' num2str(nTif) ')'])
    end
    [PTfreq{nTif},PTdBampl{nTif},PTonsetInPulse{nTif},msPTpulseLen{nTif},~,~,~] = extractMapPulseParams(S.pulse);

    clear S
end

PTfreq = cell2mat(PTfreq);
PTdBampl = cell2mat(PTdBampl);
PTonsetInPulse = cell2mat(PTonsetInPulse);
msPTpulseLen = cell2mat(msPTpulseLen);

paramS = struct2table(vertcat(paramS{:}));
paramS = paramS(:,{'stimDelay','totalPulses','ISI','PulseHiTime','PulseLoTime'});
param = unique(paramS);

if size(param,1)~=1
    error('mapping stim parameters differ between .tif files')
else
    clear paramS
end

%% calculate dFF + peak dFF for each trace for each cell
nPulsePerFile = length(PTfreq)/length(tifFileList.map);
nCell = size(tifFileList.map(1).(FsourceString),1);
%round BEFORE unique: rates differing only by float noise must collapse to
%one value, else framesPerPulse below becomes non-scalar and every frame
%index downstream is silently wrong
frameRates = extractfield(tifFileList.map,'frameRate');
fs = unique(round(frameRates));
if ~isscalar(fs)
    error('frame rates differ between map tifs: %s',mat2str(unique(frameRates)))
end
framesPerPulse = param.ISI*fs;

%PT onset frame: frame k spans t = (k-1)/fs, so the first frame at or after
%onset is find(t>=onset), NOT onset*fs (which lands one frame early).
%Matches the idiom used by processBPN2P / processCGC.
tLocalPulse = (0:framesPerPulse-1)/fs;
PTonsetIDX = arrayfun(@(o) find(tLocalPulse >= o - 1e-9,1,'first'),PTonsetInPulse);
if isempty(PTonsetIDX) || min(PTonsetIDX) < 3
    error(['PT onset lands at frame %d; need >=2 strictly pre-onset frames ' ...
        'for a baseline SD. Check onset times / frame rate.'],min(PTonsetIDX))
end

%for PT relative traces: frames before onset is kept consistent such that
%PT occurs at same x for every dFFptRel trace
%with onset of 0.6s, dFF starts from ABSidx 1, pulse on ABSidx 3
%with onset of 1s, dFF starts from ABSidx 3, pulse on ABSidx 5
%pulse occurs at ptRelIDX 3 for all then regardless of onset
maxFramesAfterOnset = framesPerPulse-max(PTonsetIDX);
maxFramesBeforeOnset = min(PTonsetIDX);

F = cat(3,tifFileList(:).map.(FsourceString));

tAbs = 0:1/fs:(size(F,2)/fs)-1/fs;
tAbs = tAbs(param.stimDelay*fs+1:param.stimDelay*fs+framesPerPulse*nPulsePerFile);
tAbsTracePulse = reshape(tAbs,framesPerPulse,nPulsePerFile);
tAbsTracePulse = repmat(tAbsTracePulse,[1 length(tifFileList.map)]);

F = F(:,param.stimDelay*fs+1:param.stimDelay*fs+framesPerPulse*nPulsePerFile,:);
F = reshape(F,nCell,framesPerPulse,nPulsePerFile*length(tifFileList.map)); %now nCell x pulseFrames x nPulse

%% obtain dFF trace for each pulse, output dFFptRel --> dFF for each pulse
%onset-aligned so PT onset sits at column maxFramesBeforeOnset for every
%pulse, regardless of that pulse's own onset time. Peaks are NOT taken here:
%significance is tested once per (cell,freq,dB) on the trial-averaged trace
%below, per pkFcalc's contract.
peakDFFbyTrial = zeros(nCell,length(PTfreq));
[rawFptRel,dFFptRel] = deal(zeros(nCell,maxFramesAfterOnset+maxFramesBeforeOnset,length(PTfreq)));
tPTrel = zeros(maxFramesAfterOnset+maxFramesBeforeOnset,length(PTfreq));

onsetCol = maxFramesBeforeOnset;    %PT onset column in the aligned traces

for pulseN = 1:length(PTfreq)
    tPTrel(:,pulseN) = tAbsTracePulse(PTonsetIDX(pulseN)-(maxFramesBeforeOnset-1):PTonsetIDX(pulseN)+maxFramesAfterOnset,pulseN);

    %F0 over STRICTLY pre-onset frames; the onset frame already carries
    %signal (a 400 ms tone spans the whole onset frame at 5 Hz)
    dFF = dFoFcalc(F(:,:,pulseN),[1 PTonsetIDX(pulseN)-1],1);

    rawFptRel(:,:,pulseN) = F(:,PTonsetIDX(pulseN)-(maxFramesBeforeOnset-1):PTonsetIDX(pulseN)+maxFramesAfterOnset,pulseN);
    dFFptRel(:,:,pulseN) = dFF(:,PTonsetIDX(pulseN)-(maxFramesBeforeOnset-1):PTonsetIDX(pulseN)+maxFramesAfterOnset);

    %retained for diagnostics only; nothing downstream tests these
    peakDFFbyTrial(:,pulseN) = max(dFFptRel(:,onsetCol:onsetCol+ ...
        nFrameWindowFor(msPTpulseLen(pulseN),fs,nFramesPostPulse)-1,pulseN),[],2);
    clear dFF
end

%% organize by freq/dB
[uFreq,~,fqIDX] = unique(PTfreq);
[uDB,~,uDBidx] = unique(PTdBampl);
dBFreqMap = cell(length(uDB),length(uFreq));
dBFreqMap(:) = {struct('tPTrel',[],'rawFptRel',[],'dFFptRel',[],'dFFavg',[],...
    'pkDFF',[],'sigPkDFF',[],'pkDFFbyTrial',[],'msPulseLen',[])};

%uDB is ascending as rows inc so fill cells according to that, then flip UD
%(70 in last row)
for pNo = 1:length(PTfreq)
    fqCol = fqIDX(pNo);
    dBrow = uDBidx(pNo);

    dBFreqMap{dBrow,fqCol}.tPTrel = cat(2,dBFreqMap{dBrow,fqCol}.tPTrel,tPTrel(:,pNo));
    dBFreqMap{dBrow,fqCol}.rawFptRel = cat(3,dBFreqMap{dBrow,fqCol}.rawFptRel,rawFptRel(:,:,pNo));
    dBFreqMap{dBrow,fqCol}.dFFptRel = cat(3,dBFreqMap{dBrow,fqCol}.dFFptRel,dFFptRel(:,:,pNo));

    dBFreqMap{dBrow,fqCol}.pkDFFbyTrial = cat(2,dBFreqMap{dBrow,fqCol}.pkDFFbyTrial,peakDFFbyTrial(:,pNo));
    dBFreqMap{dBrow,fqCol}.msPulseLen = cat(2,dBFreqMap{dBrow,fqCol}.msPulseLen,msPTpulseLen(pNo));
    clear fqCol dBrow
end

%% trial-averaged peak + significance, one test per (cell,freq,dB)
%Averaging over trials first is what pkFcalc requires: a single presentation
%gives only a handful of pre-onset frames, and averaging only the trials that
%passed a per-trial test would condition the reported response on being large.
%Trials within a cell may come from different onset groups; dFFptRel is
%onset-aligned, so that is exactly what the averaging is for.
%
%BASELINE: strictly pre-onset (1:onsetCol-1). The onset frame carries real
%signal, so pkFcalc's default 1:frameStart would inflate the threshold in
%proportion to the response being tested.
for i = 1:numel(dBFreqMap)
    c = dBFreqMap{i};
    if isempty(c.dFFptRel)
        continue    %missing freq/dB pair; filled in below
    end
    uLen = unique(c.msPulseLen);
    if ~isscalar(uLen)
        error('pulse lengths differ within one freq/dB condition: %s ms',mat2str(uLen))
    end
    nFrameWindow = nFrameWindowFor(uLen,fs,nFramesPostPulse);

    dBFreqMap{i}.dFFavg = mean(c.dFFptRel,3,'omitnan');
    [~,sig,pk] = pkFcalc(dBFreqMap{i}.dFFavg,onsetCol,nFrameWindow,...
        pkPTsigSD,1:onsetCol-1);
    dBFreqMap{i}.pkDFF = pk;            %nCell x 1, peak of the trial average
    dBFreqMap{i}.sigPkDFF = sig;        %nCell x 1 LOGICAL mask
end

%% fill missing frequency / level pairs
%built explicitly from the known nCell / nFrames rather than probing a
%sibling cell for its field sizes: that probe searched the already-filled
%array and so picked the empty pair itself, yielding NaN([0 0]) instead of
%NaN(nCell,1) and breaking the reshape below.
nFrameAligned = maxFramesAfterOnset+maxFramesBeforeOnset;
emptyIDX = cellfun(@(c) isempty(c.dFFptRel),dBFreqMap);
if any(emptyIDX,'all')
    disp('There are missing frequency / level pairs')
    emptyStruct = struct('tPTrel',NaN(nFrameAligned,0),...
        'rawFptRel',NaN(nCell,nFrameAligned,0),...
        'dFFptRel',NaN(nCell,nFrameAligned,0),...
        'dFFavg',NaN(nCell,nFrameAligned),...
        'pkDFF',NaN(nCell,1),...
        'sigPkDFF',false(nCell,1),...
        'pkDFFbyTrial',NaN(nCell,0),...
        'msPulseLen',NaN(1,0));
    dBFreqMap(emptyIDX) = deal({emptyStruct});
end

%% mean pk responses across cells
%sigPkOrNaN keeps the PEAK where significant and NaN elsewhere. The old form
%multiplied pkDFF by sigPkDFF, which held peak VALUES rather than a 0/1 mask,
%so the "significant response" maps were peak SQUARED.
sigPkOrNaN = @(c) local_sigPkOrNaN(c);

uPkResp = cellfun(@(c) mean(c.pkDFF,'omitnan'),dBFreqMap);
uSigPkResp = cellfun(@(c) mean(sigPkOrNaN(c),'omitnan'),dBFreqMap);

%% Calculate BF for each cell from peak response
%takes average of cell peak responses first, then max across
%frequency/level
dBFreqLin = numel(dBFreqMap);

%reshape maps the linear column index to sub2ind([nDB nFreq]) ordering, which
%is what anmlFRA2BF and anmlFRA2dPrime assume
sigPkCell = reshape(cell2mat(cellfun(sigPkOrNaN,dBFreqMap,'uni',0)),[nCell dBFreqLin]);
pkCell = reshape(cell2mat(cellfun(@(c) c.pkDFF,dBFreqMap,'uni',0)),[nCell dBFreqLin]);


%% send relevant vars to output
FRAoutput.sourceF = FsourceString;
FRAoutput.tifDir = tifDir;
FRAoutput.dBFreqMap = dBFreqMap;
FRAoutput.dBlist = uDB;
FRAoutput.freqList = uFreq;
FRAoutput.uPkResp = uPkResp;
FRAoutput.uSigPkResp = uSigPkResp;
FRAoutput.CellPkRespLinDBfreq = pkCell;
FRAoutput.CellSigPkLinDBfreq = sigPkCell;

%provenance: lets a stale saved _FRAmap.mat be told apart from one built by
%the trial-averaged pipeline
FRAoutput.params = struct(...
    'pkPTsigSD',pkPTsigSD,...
    'nFramesPostPulse',nFramesPostPulse,...
    'FsourceString',FsourceString,...
    'onsetCol',onsetCol,...
    'sigMethod','trialAveraged',...
    'baseline','preOnsetExclusive',...
    'onsetConvention','firstFrameAtOrAfterOnset');

% BF
FRAoutput.BFuDB = anmlFRA2BF(FRAoutput);
% dPrime
FRAoutput.dPrime = anmlFRA2dPrime(FRAoutput);

end %function

function n = nFrameWindowFor(msPulseLen,fs,nFramesPostPulse)
%frames spanned by the tone plus the post-pulse tail. ceil, not truncation:
%a bare colon floors a fractional window and silently drops a frame.
n = ceil(msPulseLen/1000*fs) + nFramesPostPulse;
end

function p = local_sigPkOrNaN(c)
%peak of the trial-averaged trace where significant, NaN otherwise
p = c.pkDFF;
p(~c.sigPkDFF) = NaN;
end

