function [sigPkResp,sig,pkResp,pkIDXsig,pkIDX] = ...
    pkFcalc(Ftrace,frameStart,nFrameWindow,pkPTsigSD)
% pkFcalc output peak response of input Ftrace within nFrameWindow from frameStart.
%   
%   [sigPkResp,sig,pkResp,pkIDXsig,pkIDX] = pkFcalc(
%           Ftrace, --> fluorescence traces (cell,vector or matrix), 
%                       assumes traces are along dimension 1
%           frameStart, --> start frame for peak search
%           nFrameWindow, --> frames from start frame for peak search
%           pkPTsigSD --> SD threshold for significant responses
%
%   OUTPUTS:
%       sigPkResp --> max responses that are pkPTsigSD above mean
%                       Ftrace before frameStart
%       sig --> logical array of significant peak responses
%       pkResp --> max response within search window
%       pkIDXsig --> index of Ftrace for significant peak responses
%       pkIDX --> index of Ftrace at which peak response occurs
%
%   See also dFoFcalc.

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