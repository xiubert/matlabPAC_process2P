function [] = intersectROIfiles(dataPath,animal,filters,tifList,tifFiles)

for ROIfileN = 1:length(filters)
    preROI{ROIfileN} = load(fullfile(dataPath,[animal '_moCorrROI_' filters{ROIfileN} '.mat']),'moCorROI');
    ROIids{ROIfileN} = string({preROI{ROIfileN}.moCorROI.ID}');  
end

allROI = vertcat(ROIids{:});
[C,~,ib] = unique(allROI,'stable');
temp = accumarray(ib,str2double(allROI),[],@numel)==length(filters);
ROIinALL = C(temp);
clear C ib allROI temp ROIfileN ROIids

for ROIfileN = 1:length(filters)
    ROIkeep{ROIfileN} = ismember(string({preROI{ROIfileN}.moCorROI.ID}'),ROIinALL);
    if ~all(ROIkeep{ROIfileN})
        nTifs = length(tifList.(filters{ROIfileN}));
        tifIDXinAllTifList = ismember({tifFiles.name}',{tifList.(filters{ROIfileN}).name}');
        matchNroiS{ROIfileN} = preROI{ROIfileN}.moCorROI(ROIkeep{ROIfileN});
        movefile(fullfile(dataPath,[animal '_moCorrROI_' filters{ROIfileN} '.mat']),...
            fullfile(dataPath,[animal '_OLDmoCorrROI_' filters{ROIfileN} '.mat']))
        moCorROI = matchNroiS{ROIfileN};
        save(fullfile(dataPath,[animal '_moCorrROI_' filters{ROIfileN} '.mat']),...
            'moCorROI','nTifs','tifIDXinAllTifList');
        clear moCorROI
    end
end