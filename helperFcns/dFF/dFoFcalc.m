function dFoF = dFoFcalc(rawF,baseIDX,traceDim)
% dFoFcalc output deltaF/F for a vector, matrix, or cell array of raw F.
%   dFoF = dFoFcalc(
%           rawF, --> raw fluorescence as vector, matrix, or cell array
%           baseIDX, --> [first last] indices of rawF to be used for F0 calculation
%           traceDim --> dimension along traces
%
%   Calculates deltaF/F from rawF input with F0 as mean of rawF within
%   baseIDX indices. dF/F is calculated for each trace. dF/F trace starts
%   at first baseIDX
%
%   See also baseIDXfromPTonset.m, pkFcalc.m

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


