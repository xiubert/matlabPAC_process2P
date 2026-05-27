function dFoF = dFoFcalc(rawF,baseIDX,traceDim)
% DFOFCALC  Compute dF/F from raw fluorescence traces.
%
%   dFoF = dFoFcalc(rawF, baseIDX)
%   dFoF = dFoFcalc(rawF, baseIDX, traceDim)
%
%   Computes (F - F0) / F0 for each trace, where F0 is the mean of rawF
%   over the baseline window baseIDX(1):baseIDX(2).  The output trace
%   begins at baseIDX(1) — pre-baseline frames are discarded.
%
%   Inputs:
%     rawF      - raw fluorescence; vector, nTraces x nFrames matrix, or
%                 cell array of such matrices
%     baseIDX   - [first last] frame indices defining the F0 baseline window
%     traceDim  - dimension along which individual traces run:
%                   1 (default) → traces are rows, frames are columns
%                   2           → frames are rows, traces are columns
%
%   Output:
%     dFoF - dF/F in the same format as rawF, starting at frame baseIDX(1)
%
%   See also baseIDXfromPTonset, pkFcalc

if nargin==2
    traceDim=1;
end

if traceDim==1
    if ~iscell(rawF)
        dFoF = (rawF(:,baseIDX(1):end)-mean(rawF(:,baseIDX(1):baseIDX(2)),2))./...
            mean(rawF(:,baseIDX(1):baseIDX(2)),2);
    else
        dFoF = cellfun(@(c) (c(:,baseIDX(1):end)-mean(c(:,baseIDX(1):baseIDX(2)),2))./...
            mean(c(:,baseIDX(1):baseIDX(2)),2),rawF,'uni',0);
    end
    
elseif traceDim==2
    if ~iscell(rawF)
        dFoF = (rawF(baseIDX(1):end,:)-mean(rawF(baseIDX(1):baseIDX(2),:),1))./...
            mean(rawF(baseIDX(1):baseIDX(2),:),1);
    else
        dFoF = cellfun(@(c) (c(baseIDX(1):end,:)-mean(c(baseIDX(1):baseIDX(2),:),1))./...
            mean(c(baseIDX(1):baseIDX(2),:),1),rawF,'uni',0);
    end
end

end %function


