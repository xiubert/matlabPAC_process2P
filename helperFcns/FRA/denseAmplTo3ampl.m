%denseAmplTo3ampl
function TfraAnml = denseAmplTo3ampl(TfraAnml)
%REMOVE IF FUTURE MAPS WILL BE > 3 LVLs
if any(cellfun(@length,{TfraAnml.dBlist})>3)
    IDXdense = find(cellfun(@length,{TfraAnml.dBlist})>3);
    tmpMap = TfraAnml(IDXdense).dBFreqMap(ismember(TfraAnml(IDXdense).dBlist,[30 50 70]),:);
    tmp = cellfun(@(c) nanmean(zero2nan(c.pkDFF.*c.sigPkDFF),2),tmpMap,'uni',0);
    TfraAnml(IDXdense).CellSigPkLinDBfreq = reshape(cell2mat(tmp),[length(TfraAnml(IDXdense).BFuDB) numel(tmp)]);
    tmp = cellfun(@(c) nanmean(zero2nan(c.pkDFF.*c.pkDFF),2),tmpMap,'uni',0);
    TfraAnml(IDXdense).CellPkRespLinDBfreq = reshape(cell2mat(tmp),[length(TfraAnml(IDXdense).BFuDB) numel(tmp)]);
    
    TfraAnml(IDXdense).uPkResp = TfraAnml(IDXdense).uPkResp(ismember(TfraAnml(IDXdense).dBlist,[30 50 70]),:);
    TfraAnml(IDXdense).uSigPkResp = TfraAnml(IDXdense).uSigPkResp(ismember(TfraAnml(IDXdense).dBlist,[30 50 70]),:);
    
    TfraAnml(IDXdense).dBFreqMap = tmpMap;
    TfraAnml(IDXdense).dBlist = [30; 50; 70];
end