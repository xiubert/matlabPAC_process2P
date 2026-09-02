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

%% obtain dFF trace for each pulse, extract peak, output dFFptRel --> dFF for each pulse w/ PT onset at IDX 3 (first post pulse frame is 4)
[peakDFF,sigPkDff] = deal(zeros(nCell,length(PTfreq)));
[rawFptRel,dFFptRel] = deal(zeros(nCell,maxFramesAfterOnset+maxFramesBeforeOnset,length(PTfreq)));
tPTrel = deal(zeros(maxFramesAfterOnset+maxFramesBeforeOnset,length(PTfreq)));

for pulseN = 1:length(PTfreq)
    tPTrel(:,pulseN) = tAbsTracePulse(PTonsetIDX(pulseN)-(maxFramesBeforeOnset-1):PTonsetIDX(pulseN)+maxFramesAfterOnset,pulseN);
    
    %F0 over STRICTLY pre-onset frames; the onset frame already carries
    %signal (a 400 ms tone spans the whole onset frame at 5 Hz)
    dFF = dFoFcalc(F(:,:,pulseN),[1 PTonsetIDX(pulseN)-1],1);
    [sigPkDff_tmp,sig_tmp,peakDFF(:,pulseN),~,~] = ...        
        pkFcalc(dFF,PTonsetIDX(pulseN),...
        msPTpulseLen(pulseN)/1000*fs+nFramesPostPulse,pkPTsigSD);
    sigPkDff(sig_tmp,pulseN) = sigPkDff_tmp;
%     [sigPkDff(:,pulseN),~,peakDFF(:,pulseN),~,~] = ...        
%         pkFcalc(dFF,PTonsetIDX(pulseN),...
%         msPTpulseLen(pulseN)/1000*fs+nFramesPostPulse,pkPTsigSD);
%     
    rawFptRel(:,:,pulseN) = F(:,PTonsetIDX(pulseN)-(maxFramesBeforeOnset-1):PTonsetIDX(pulseN)+maxFramesAfterOnset,pulseN);
    dFFptRel(:,:,pulseN) = dFF(:,PTonsetIDX(pulseN)-(maxFramesBeforeOnset-1):PTonsetIDX(pulseN)+maxFramesAfterOnset);
    clear dFF sig
end

%% organize by freq/dB
[uFreq,~,fqIDX] = unique(PTfreq);
[uDB,~,uDBidx] = unique(PTdBampl);
dBFreqMap = cell(length(uDB),length(uFreq));
dBFreqMap(:) = {struct('tPTrel',[],'rawFptRel',[],'dFFptRel',[],'pkDFF',[],'sigPkDFF',[])};

%uDB is ascending as rows inc so fill cells according to that, then flip UD
%(70 in last row)
for pNo = 1:length(PTfreq)
    fqCol = fqIDX(pNo);
    dBrow = uDBidx(pNo);
    
    dBFreqMap{dBrow,fqCol}.tPTrel = cat(2,dBFreqMap{dBrow,fqCol}.tPTrel,tPTrel(:,pNo));
    dBFreqMap{dBrow,fqCol}.rawFptRel = cat(3,dBFreqMap{dBrow,fqCol}.rawFptRel,rawFptRel(:,:,pNo));
    dBFreqMap{dBrow,fqCol}.dFFptRel = cat(3,dBFreqMap{dBrow,fqCol}.dFFptRel,dFFptRel(:,:,pNo));
    
    dBFreqMap{dBrow,fqCol}.pkDFF = cat(2,dBFreqMap{dBrow,fqCol}.pkDFF,peakDFF(:,pNo));
    dBFreqMap{dBrow,fqCol}.sigPkDFF = cat(2,dBFreqMap{dBrow,fqCol}.sigPkDFF,sigPkDff(:,pNo));
    clear fqCol dBrow
end

%% mean pk responses across cells
uPkResp = cellfun(@nanmean,cellfun(@(c) nanmean(c.pkDFF,2),dBFreqMap,'uni',0));
uSigPkResp = cellfun(@(c) nanmean(nanmean(zero2nan(c.pkDFF.*c.sigPkDFF),2)),dBFreqMap);

%% Calculate BF for each cell from peak response
%takes average of cell peak responses first, then max across
%frequency/level
dBFreqLin = numel(dBFreqMap);

tmp = cellfun(@(c) nanmean(zero2nan(c.pkDFF.*c.sigPkDFF),2),dBFreqMap,'uni',0);
%deal with missing frequency / level pairs
if any(cellfun(@isempty,tmp),'all')
    disp('There are missing frequency / level pairs')
    emptyIDX = cellfun(@isempty,tmp);
    tmp(emptyIDX) = deal({NaN(nCell,1)});
    tmpS = dBFreqMap{find(~cellfun(@isempty,tmp),1)};
    emptyStruct = cell2struct(cellfun(@(c) NaN(c),...
        varfun(@(v) size(v{1}),struct2table(tmpS,'AsArray',1),...
        'OutputFormat','cell'),'uni',0),fieldnames(tmpS),2);
    dBFreqMap(emptyIDX) = deal({emptyStruct});
    clear tmpS
end
sigPkCell = reshape(cell2mat(tmp),[nCell dBFreqLin]);
pkCell = reshape(cell2mat(cellfun(@(c) nanmean(c.pkDFF,2),dBFreqMap,'uni',0)),[nCell dBFreqLin]);


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

% BF
FRAoutput.BFuDB = anmlFRA2BF(FRAoutput);
% dPrime
FRAoutput.dPrime = anmlFRA2dPrime(FRAoutput);

