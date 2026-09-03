%denseAmplTo3ampl
function TfraAnml = denseAmplTo3ampl(TfraAnml)
%REMOVE IF FUTURE MAPS WILL BE > 3 LVLs
if any(cellfun(@length,{TfraAnml.dBlist})>3)
    IDXdense = find(cellfun(@length,{TfraAnml.dBlist})>3);
    tmpMap = TfraAnml(IDXdense).dBFreqMap(ismember(TfraAnml(IDXdense).dBlist,[30 50 70]),:);
    %peak where significant, NaN elsewhere -- matching FRAmap. The old form
    %multiplied pkDFF by sigPkDFF, which held peak VALUES rather than a 0/1
    %mask, so this map was peak SQUARED.
    tmp = cellfun(@(c) local_sigPkOrNaN(c),tmpMap,'uni',0);
    TfraAnml(IDXdense).CellSigPkLinDBfreq = reshape(cell2mat(tmp),[length(TfraAnml(IDXdense).BFuDB) numel(tmp)]);
    %pkDFF.*pkDFF squared the ALL-response map too, disagreeing with FRAmap
    %and so silently rescaling any animal with >3 levels relative to the rest
    %of the cohort
    tmp = cellfun(@(c) c.pkDFF,tmpMap,'uni',0);
    TfraAnml(IDXdense).CellPkRespLinDBfreq = reshape(cell2mat(tmp),[length(TfraAnml(IDXdense).BFuDB) numel(tmp)]);
    
    TfraAnml(IDXdense).uPkResp = TfraAnml(IDXdense).uPkResp(ismember(TfraAnml(IDXdense).dBlist,[30 50 70]),:);
    TfraAnml(IDXdense).uSigPkResp = TfraAnml(IDXdense).uSigPkResp(ismember(TfraAnml(IDXdense).dBlist,[30 50 70]),:);
    
    TfraAnml(IDXdense).dBFreqMap = tmpMap;
    TfraAnml(IDXdense).dBlist = [30; 50; 70];
end
end %function

function p = local_sigPkOrNaN(c)
%peak of the trial-averaged trace where significant, NaN otherwise
p = c.pkDFF;
p(~c.sigPkDFF) = NaN;
end
