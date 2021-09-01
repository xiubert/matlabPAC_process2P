function tifFileList = FISSAoutput2tifFileList(FISSAoutput,tifFileList,fissaScaleFactor)
tifListFields = fieldnames(tifFileList);
for k = 1:length(tifListFields)
    if ~isempty(tifFileList.(tifListFields{k}))
        
        %get first row from each trial in FISSAoutput for each cell
        %first row corresponds to cell fluorescence remaining rows are
        %fluorescence from neuropil regions
        tmp = struct2cell(structfun(@(F) structfun(@(f) f(1,:), F, 'UniformOutput', false), ...
            FISSAoutput.(tifListFields{k}).result))';
        tmp1 = arrayfun(@(colIDX) vertcat(tmp{:,colIDX}), 1:size(tmp,2),'uni',0)';
        [tifFileList.(tifListFields{k}).fissaFroi] = tmp1{:};
        
        %scale fissa subtraction by fissaScaleFactor
        tmp2 = arrayfun(@(r) r.moCorRawFroi - fissaScaleFactor.*(r.moCorRawFroi-r.fissaFroi),...
            tifFileList.(tifListFields{k}), 'UniformOutput', false);
        [tifFileList.(tifListFields{k}).SCALEDfissaFroi] = tmp2{:};
        
        clear tmp*
    end
end

