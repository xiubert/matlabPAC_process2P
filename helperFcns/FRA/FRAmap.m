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
addParameter(p,'nBaselineFrames',[],@(x) isempty(x)||(isscalar(x)&&x>=2));
addParameter(p,'baselineDecayMargin',2,@(x) isscalar(x)&&x>=0);

parse(p,tifFileList,pkPTsigSD,nFramesPostPulse,varargin{:});

tifFileList = p.Results.tifFileList;
pkPTsigSD = p.Results.pkPTsigSD;
nFramesPostPulse = p.Results.nFramesPostPulse;
FsourceString = p.Results.FsourceString;
nBaselineFrames = p.Results.nBaselineFrames;
baselineDecayMargin = p.Results.baselineDecayMargin;

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
%PT occurs at same x for every dFFptRel trace, regardless of that pulse's
%own onset time
maxFramesAfterOnset = framesPerPulse-max(PTonsetIDX);

%Keep the FULL per-tif array. The baseline below reaches BACKWARDS out of a
%pulse window into the preceding silence, so cropping to the stim period
%first (as this used to) would throw away exactly the frames it needs.
F = cat(3,tifFileList(:).map.(FsourceString));
nTif = length(tifFileList.map);
tAbsFull = (0:size(F,2)-1)/fs;

%Absolute within-tif frame of each pulse's tone onset. Pulses are contiguous
%segments of one recording, so absolute indexing is what lets a baseline span
%a pulse boundary.
nFrameWindowPulse = ceil(msPTpulseLen/1000*fs) + nFramesPostPulse;
[absOnset,tifOfPulse] = deal(zeros(length(PTfreq),1));
for pulseN = 1:length(PTfreq)
    tifOfPulse(pulseN) = ceil(pulseN/nPulsePerFile);
    pLocal = pulseN - (tifOfPulse(pulseN)-1)*nPulsePerFile;
    absOnset(pulseN) = param.stimDelay*fs + (pLocal-1)*framesPerPulse + PTonsetIDX(pulseN);
end

%% baseline length
%The ISI (%g s here) is much longer than the response, so each pulse ends in
%genuine silence. Reaching back across the pulse boundary buys a far better
%SD estimate than the 2-3 pre-onset frames inside the window: measured on
%TO0003 it lifts the real-vs-sham discrimination from 1.42x to 1.53x, and it
%lowers the false-positive rate (sham 0.33 -> 0.25).
%
%The limit is the PREVIOUS tone's response, which must stay out of the
%baseline. For the first pulse of a tif the limit is instead the start of the
%recording, and the stimDelay period before it is silent.
if isempty(nBaselineFrames)
    nBase = Inf;
    for pulseN = 1:length(PTfreq)
        pLocal = pulseN - (tifOfPulse(pulseN)-1)*nPulsePerFile;
        if pLocal == 1
            lim = absOnset(pulseN)-1;                       % into the stim delay
        else
            prevOnset = absOnset(pulseN-1);
            lim = absOnset(pulseN) - (prevOnset + nFrameWindowPulse(pulseN-1) + baselineDecayMargin);
        end
        nBase = min(nBase,lim);
    end
    %never shorter than the original within-pulse baseline
    nBase = max(nBase, min(PTonsetIDX)-1);
else
    nBase = nBaselineFrames;
end
if nBase < 2
    error('FRAmap:baselineTooShort', ...
        'baseline resolved to %d frames; need >=2 for an SD',nBase);
end
if any(absOnset-nBase < 1)
    error('FRAmap:baselineBeforeRecording', ...
        'a %d-frame baseline runs before the start of the recording',nBase);
end

%onset sits at this column of every aligned trace
onsetCol = nBase + 1;

%% obtain dFF trace for each pulse, output dFFptRel --> dFF for each pulse
%Each pulse is extracted as ONE contiguous segment: nBase silent frames, the
%onset frame, then maxFramesAfterOnset. The segment is onset-aligned by
%construction, so PT onset sits at column onsetCol for every pulse regardless
%of its own onset time, and the baseline that F0 and the significance test
%both use is the leading nBase columns.
%
%Peaks are NOT taken here: significance is tested once per (cell,freq,dB) on
%the trial-averaged trace below, per pkFcalc's contract.
alignedLen = nBase + 1 + maxFramesAfterOnset;
peakDFFbyTrial = zeros(nCell,length(PTfreq));
[rawFptRel,dFFptRel] = deal(zeros(nCell,alignedLen,length(PTfreq)));
tPTrel = zeros(alignedLen,length(PTfreq));

for pulseN = 1:length(PTfreq)
    span = absOnset(pulseN)-nBase : absOnset(pulseN)+maxFramesAfterOnset;
    seg  = F(:,span,tifOfPulse(pulseN));

    %F0 over the extended STRICTLY pre-onset baseline. The onset frame is
    %excluded: a 400 ms tone spans the whole onset frame at 5 Hz, so it
    %already carries signal.
    F0 = mean(seg(:,1:nBase),2);

    tPTrel(:,pulseN)      = tAbsFull(span)';
    rawFptRel(:,:,pulseN) = seg;
    dFFptRel(:,:,pulseN)  = (seg - F0)./F0;

    %retained for diagnostics only; nothing downstream tests these
    peakDFFbyTrial(:,pulseN) = max(dFFptRel(:,onsetCol:onsetCol+ ...
        nFrameWindowPulse(pulseN)-1,pulseN),[],2);
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
nFrameAligned = alignedLen;
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
    'nBaselineFrames',nBase,...
    'baselineDecayMargin',baselineDecayMargin,...
    'alignedLen',alignedLen,...
    'sigMethod','trialAveraged',...
    'baseline','preOnsetExtended',...
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

