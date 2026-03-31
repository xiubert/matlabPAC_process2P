% No longer in use. See BPN2.m


% function [BPNoutput,BPNdBampl] = BPN(tifFileList,...
%     pkBPNsigSD,nFramesPostPulse,varargin)
% % BPN: output BPN for an animal.
% %   [BPNoutput,BPNdBampl] = BPN(
% %           tifFileList, --> list of BPN tif files w/ F signal
% %           pkBPNsigSD, --> SD threshold for a significant BPN response
% %           nFramesPostPulse, --> frames after pulse ends for peak search
% %           (eg if pulse is 400 ms and FR is 5, total frames post pulse
% %           onset will be 2 + nFramesPostPulse
% %           varargin --> optionally define Fsource string
% %
% %   output: structure with ?
% %   BPNdBampl: list of BPN dB amplitudes in order of presentation
% 
% % Get pulse params for each tif
% p = inputParser;
% addRequired(p,'tifFileList',@isstruct)
% addRequired(p,'pkBPNsigSD',@isscalar)
% addRequired(p,'nFramesPostPulse',@isscalar)
% addOptional(p,'FsourceString','SCALEDfissaFroi',@ischar);
% 
% parse(p,tifFileList,pkPTsigSD,nFramesPostPulse,varargin{:});
% 
% tifFileList = p.Results.tifFileList;
% pkBPNsigSD = p.Results.pkBPNsigSD;
% nFramesPostPulse = p.Results.nFramesPostPulse;
% FsourceString = p.Results.FsourceString;
%% test from here
pkBPNsigSD=2;
nFramesPostPulse=4;
FsourceString='SCALEDfissaFroi';
%
tifDir = tifFileList.BPN(1).folder;
[BPNfreqcutoff,BPNdBampl,BPNonsetInPulse,msBPNpulseLen,paramS] = deal(cell(length(tifFileList.BPN),1));
for nTif = 1:length(tifFileList.BPN)
    try
        S = load(fullfile(tifDir,strrep(tifFileList.BPN(nTif).name,'.tif','_Pulses.mat')));
    catch
        disp('Can''t find pulses associated with .tifs, locate dir...')
        tifDir = uigetdir('D:\Data',...
            'Locate directory containing tifs...');
        S = load(fullfile(tifDir,strrep(tifFileList.BPN(nTif).name,'.tif','_Pulses.mat')));
    end
    paramS{nTif} = S.params;
    if length(S.pulse)~=S.params.totalPulses
        error(['pulse number lower or higher than noted in pulse params (tif file no. ' num2str(nTif) ')'])
    end
    [BPNfreqcutoff{nTif},BPNdBampl{nTif},BPNonsetInPulse{nTif},msBPNpulseLen{nTif},~] = extractBPNPulseParams(S.pulse);

    clear S
end

BPNdBampl = cell2mat(BPNdBampl);
BPNonsetInPulse = cell2mat(BPNonsetInPulse);
msBPNpulseLen = cell2mat(msBPNpulseLen);

paramS = struct2table(vertcat(paramS{:}));
paramS = paramS(:,{'stimDelay','totalPulses','ISI','PulseHiTime','PulseLoTime'});
param = unique(paramS);

if size(param,1)~=1
    error('stim parameters differ between .tif files')
else
    clear paramS
end

%% calculate dFF + peak dFF for each trace for each cell
nPulsePerFile = length(BPNdBampl)/length(tifFileList.BPN);%% try to extract from elsewhere
nCell = size(tifFileList.BPN(1).(FsourceString),1);
fs = unique(extractfield(tifFileList.BPN,'frameRate'));
framesPerPulse = param.ISI*fs;
BPNonsetIDX = BPNonsetInPulse*fs;

%frames before onset is kept consistent such that BPN occurs at same x for every dFFptRel trace
%with onset of 0.6s, dFF starts from ABSidx 1, pulse on ABSidx 3
%with onset of 1s, dFF starts from ABSidx 3, pulse on ABSidx 5

F = cat(3,tifFileList(:).BPN.(FsourceString));

tAbs = 0:1/fs:(size(F,2)/fs)-1/fs;
tAbs = tAbs(param.stimDelay*fs+1:param.stimDelay*fs+framesPerPulse*nPulsePerFile);
tAbsTracePulse = reshape(tAbs,framesPerPulse,nPulsePerFile);
tAbsTracePulse = repmat(tAbsTracePulse,[1 length(tifFileList.BPN)]);

F = F(:,param.stimDelay*fs+1:param.stimDelay*fs+framesPerPulse*nPulsePerFile,:);
F = reshape(F,nCell,framesPerPulse,nPulsePerFile*length(tifFileList.BPN)); %now nCell x pulseFrames x nPulse

%% obtain dFF trace for each pulse, extract peak, output dFFptRel --> dFF for each pulse w/ BPN onset at ABSidx 5 (first post pulse frame is 6)
[peakDFF,sigPkDff] = deal(nan(nCell,length(BPNdBampl)));
[rawFBPN,dFFBPN] = deal(zeros(nCell,framesPerPulse,length(BPNdBampl)));

for pulseN = 1:length(BPNdBampl) 
    dFF = dFoFcalc(F(:,:,pulseN),[1 BPNonsetIDX(pulseN)],1);
    [sigPkDff_tmp,sig_tmp,peakDFF(:,pulseN),~,~] = ...        
        pkFcalc(dFF,BPNonsetIDX(pulseN),...
        msBPNpulseLen(pulseN)/1000*fs+nFramesPostPulse,pkBPNsigSD);
    sigPkDff(sig_tmp,pulseN) = sigPkDff_tmp;    
    rawFBPN(:,:,pulseN) = F(:,:,pulseN);
    dFFBPN(:,:,pulseN) = dFF(:,:);
    clear dFF sig
end

%% organize by dB
tBPN=tAbsTracePulse;
n=length(BPNdBampl);
% [uDB,~] = unique(BPNdBampl);
animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');
animal_tmp = repmat(animal,[n*nCell 1]);
roiID = repelem((1:nCell).', n);
MeasurementID = repmat((1:n).', nCell, 1);
SoundLevel = repmat(BPNdBampl, nCell, 1);
Peakdff_tmp = reshape(peakDFF.', [], 1);
sigPkDff_tmp = reshape(sigPkDff.', [], 1);
% t_tmp = reshape(tBPN', [], 20);
t_cell = cell(nCell*n,1);
dff_cell = cell(nCell*n,1);
rawf_cell = cell(nCell*n,1);
% Tlen = size(tBPN,1);
for j = 1:n
    tj = tBPN(:,j);
    for i = 1:nCell
        rowIdx = (i-1)*n + j;
        t_cell{rowIdx} = tj;
        rawf_cell{rowIdx} = rawFBPN(i,:,j).';
        dff_cell{rowIdx} = dFFBPN(i,:,j).';
    end
end
Ttbl = table(animal_tmp,roiID, MeasurementID, SoundLevel, t_cell, rawf_cell, dff_cell, Peakdff_tmp, sigPkDff_tmp, ...
    'VariableNames', {'animal','roiID','MeasurementID','dB','time','rawF','dFF','PeakdFF','sigPeakdFF'});

%% 
G = findgroups(Ttbl.roiID, Ttbl.dB);   % group by Cell and Sound level
% Aggregate measurement IDs (numeric vector per group)
MeasurementIDs = splitapply(@(x) {x.'}, Ttbl.MeasurementID, G);

% Aggregate peak values (numeric vector per group)
PeakList = splitapply(@(x) {x.'}, Ttbl.PeakdFF, G);
sigPeakList = splitapply(@(x) {x.'}, Ttbl.sigPeakdFF, G);
% Mean peak per group (scalar)
MeanPeak = splitapply(@(x) {mean(x,'omitnan')}, Ttbl.PeakdFF, G);
MeansigPeak = splitapply(@(x) {mean(x,'omitnan')}, Ttbl.sigPeakdFF, G);
% Aggregate Time and DFF: return cell arrays of cell arrays (each entry is a cell array of vectors)
timeList = splitapply(@(c) {c}, Ttbl.time, G);   % gives cell array where each element is an n_k x 1 cell array of time vectors
dFFList  = splitapply(@(c) {c}, Ttbl.dFF,  G);

% Also collect the grouping keys (CellID and SoundLevel)
[roiID_group, SoundLevel_group] = splitapply(@(a,b) deal(a(1), b(1)), Ttbl.roiID, Ttbl.dB, G);
% splitapply returned cell arrays of scalars, convert to numeric vectors:
% roiID_group = cell2mat(grouproiID).';
% SoundLevel_group = cell2mat(groupSound).';

% Build grouped table
GroupedTbl = table(animal_tmp(1:6*nCell,:), roiID_group, SoundLevel_group, MeasurementIDs, PeakList, MeanPeak, timeList, dFFList, sigPeakList, MeansigPeak,...
    'VariableNames', {'animal','roiID','dB','MeasurementIDs','Peakdff','MeanPeak','time','dFF','sigPeakList', 'MeansigPeak'});
save(fullfile(dataPath,[animal '_anmlROI_BPNTable.mat']),'GroupedTbl','-v7.3');
%% send relevant vars to output
% BPNoutput.sourceF = FsourceString;
% BPNoutput.tifDir = tifDir;
% BPNoutput.dBMap = dBMap;
% BPNoutput.dBlist = uDB;
% BPNoutput.uPkResp = uPkResp;
% BPNoutput.uSigPkResp = uSigPkResp;
% BPNoutput.CellPkRespLinDB = pkCell;

% sigPkCell = reshape(cell2mat(MeansigPeak), nCell,6);
% BPNoutput.CellSigPkLinDB = sigPkCell;
% BPNoutput.dPrime = dPrime_ampOnly(BPNoutput);

