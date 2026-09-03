function [sigPkResp,sig,pkResp,pkIDXsig,pkIDX] = ...
    pkFcalc(Ftrace,frameStart,nFrameWindow,pkPTsigSD,baseIDX)
% PKFCALC  Find peak fluorescence response and flag significant responses.
%
%   [sigPkResp, sig, pkResp, pkIDXsig, pkIDX] = pkFcalc(Ftrace, frameStart, nFrameWindow, pkPTsigSD)
%   [...] = pkFcalc(..., baseIDX)
%
%   Searches for the peak value of each trace within the window
%   [frameStart, frameStart + nFrameWindow - 1].  A response is considered
%   significant if its peak exceeds the baseline mean plus
%   pkPTsigSD * baseline SD, taken over the frames baseIDX.
%
%   NOTE: Ftrace should contain per-cell average traces (e.g. trial-averaged
%   dF/F), not individual pixel or trial traces.  Significance thresholding
%   based on baseline SD is not meaningful for single-trial data.
%
%   BASELINE WINDOW: baseIDX defaults to 1:frameStart, i.e. the onset frame
%   sits in BOTH the baseline and the peak-search window.  That overlap is
%   deliberate for BPN (processBPN2P) and CGC (processCGC) — see the design
%   note at processCGC.m, which wants the baseline to run up to and including
%   PT onset — and those callers use the 4-argument form.
%
%   FRA (FRAmap) passes an explicit STRICTLY PRE-ONSET window,
%   baseIDX = 1:frameStart-1.  There the onset frame already carries signal
%   (a 400 ms tone spans the whole 200 ms onset frame at 5 Hz), so including
%   it inflates the threshold in proportion to the very response being
%   tested, and rejects clearly-responding cells.
%
%   Inputs:
%     Ftrace        - fluorescence traces; nTraces x nFrames matrix, or a
%                     cell array of such matrices (traces along dimension 2)
%     frameStart    - first frame of the post-stimulus search window (1-based)
%     nFrameWindow  - number of frames in the search window
%     pkPTsigSD     - significance threshold in units of baseline SD
%     baseIDX       - (optional) frame indices defining the baseline window.
%                     Defaults to 1:frameStart (legacy behaviour).
%
%   Outputs:
%     sigPkResp - peak responses for significant traces only (subset of pkResp)
%     sig       - nTraces x 1 logical vector; true where peak >= baseline + pkPTsigSD*SD
%     pkResp    - nTraces x 1 peak response within the search window (all traces)
%     pkIDXsig  - frame indices of peak within the search window for significant traces
%     pkIDX     - frame indices of peak within the search window for all traces
%
%   See also dFoFcalc, zero2nan, FRAmap

if nargin < 5 || isempty(baseIDX)
    baseIDX = 1:frameStart;
end

%frame span for peak search
maxSpan(1) = frameStart;
maxSpan(2) = maxSpan(1) + nFrameWindow-1;
maxSpan = maxSpan(1):maxSpan(2);

if iscell(Ftrace)
    [pkResp,pkIDX] = cellfun(@(c) max(c(:,maxSpan),[],2),Ftrace,'uni',0);
    
    %mean of F trace over baseIDX plus pkPTsigSD*std of F trace over baseIDX
    %respective to each trace
    baselinePlusSD = cellfun(@(c) zero2nan(nanmean(c(:,baseIDX),2)+(pkPTsigSD.*...
        nanstd(c(:,baseIDX),0,2))),Ftrace,'uni',0);
    
    x = cell2struct([pkResp baselinePlusSD pkIDX],{'pk','baselinePlusSD','pkIDX'},2);
    sig = arrayfun(@(r) r.pk>=r.baselinePlusSD,x,'uni',0);
    sigPkResp = arrayfun(@(r) r.pk(r.pk>=r.baselinePlusSD),x,'uni',0);
    pkIDXsig = arrayfun(@(r) r.pkIDX(r.pk>=r.baselinePlusSD),x,'uni',0);
else
    [pkResp,pkIDX] = max(Ftrace(:,maxSpan),[],2);
    baselinePlusSD = zero2nan(nanmean(Ftrace(:,baseIDX),2)+(pkPTsigSD.*...
        nanstd(Ftrace(:,baseIDX),0,2)));
    sig = pkResp>=baselinePlusSD;
    sigPkResp = pkResp(sig);
    pkIDXsig = pkIDX(sig);
end
