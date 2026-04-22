function [sigPkResp,sig,pkResp,pkIDXsig,pkIDX] = ...
    pkFcalc(Ftrace,frameStart,nFrameWindow,pkPTsigSD)
% PKFCALC  Find peak fluorescence response and flag significant responses.
%
%   [sigPkResp, sig, pkResp, pkIDXsig, pkIDX] = pkFcalc(Ftrace, frameStart, nFrameWindow, pkPTsigSD)
%
%   Searches for the peak value of each trace within the window
%   [frameStart, frameStart + nFrameWindow - 1].  A response is considered
%   significant if its peak exceeds the pre-stimulus baseline mean plus
%   pkPTsigSD * baseline SD (computed over frames 1:frameStart).
%
%   NOTE: Ftrace should contain per-cell average traces (e.g. trial-averaged
%   dF/F), not individual pixel or trial traces.  Significance thresholding
%   based on baseline SD is not meaningful for single-trial data.
%
%   Inputs:
%     Ftrace        - fluorescence traces; nTraces x nFrames matrix, or a
%                     cell array of such matrices (traces along dimension 2)
%     frameStart    - first frame of the post-stimulus search window (1-based)
%     nFrameWindow  - number of frames in the search window
%     pkPTsigSD     - significance threshold in units of baseline SD
%
%   Outputs:
%     sigPkResp - peak responses for significant traces only (subset of pkResp)
%     sig       - nTraces x 1 logical vector; true where peak >= baseline + pkPTsigSD*SD
%     pkResp    - nTraces x 1 peak response within the search window (all traces)
%     pkIDXsig  - frame indices of peak within the search window for significant traces
%     pkIDX     - frame indices of peak within the search window for all traces
%
%   See also dFoFcalc, zero2nan

%frame span for peak search
maxSpan(1) = frameStart;
maxSpan(2) = maxSpan(1) + nFrameWindow-1;
maxSpan = maxSpan(1):maxSpan(2);

if iscell(Ftrace)
    [pkResp,pkIDX] = cellfun(@(c) max(c(:,maxSpan),[],2),Ftrace,'uni',0);
    
    %mean of F trace up to frameStart plus pkPTsigSD*std of F trace up to frameStart
    %respective to each trace
    baselinePlusSD = cellfun(@(c) zero2nan(nanmean(c(:,1:frameStart),2)+(pkPTsigSD.*...
        nanstd(c(:,1:frameStart),0,2))),Ftrace,'uni',0);
    
    x = cell2struct([pkResp baselinePlusSD pkIDX],{'pk','baselinePlusSD','pkIDX'},2);
    sig = arrayfun(@(r) r.pk>=r.baselinePlusSD,x,'uni',0);
    sigPkResp = arrayfun(@(r) r.pk(r.pk>=r.baselinePlusSD),x,'uni',0);
    pkIDXsig = arrayfun(@(r) r.pkIDX(r.pk>=r.baselinePlusSD),x,'uni',0);
else
    [pkResp,pkIDX] = max(Ftrace(:,maxSpan),[],2);
    baselinePlusSD = zero2nan(nanmean(Ftrace(:,1:frameStart),2)+(pkPTsigSD.*...
        nanstd(Ftrace(:,1:frameStart),0,2)));
    sig = pkResp>=baselinePlusSD;
    sigPkResp = pkResp(sig);
    pkIDXsig = pkIDX(sig);
end