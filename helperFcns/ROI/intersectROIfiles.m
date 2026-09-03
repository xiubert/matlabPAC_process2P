function [] = intersectROIfiles(roiDir,animal,filters,tifList,tifFiles)
% intersectROIfiles  Keep only the ROIs present in every condition.
%
%   intersectROIfiles(roiDir,animal,filters,tifList,tifFiles)
%
%   roiDir is the folder holding <animal>_moCorrROI_<cond>.mat -- the animal
%   folder in the flat layout, or analysis/<run>/ otherwise. It is NOT
%   necessarily the animal folder, which is why it is not called dataPath.

for ROIfileN = 1:length(filters)
    preROI{ROIfileN} = load(fullfile(roiDir,[animal '_moCorrROI_' filters{ROIfileN} '.mat']),'moCorROI');
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
        roiFile = fullfile(roiDir,[animal '_moCorrROI_' filters{ROIfileN} '.mat']);
        % keep everything else the file carried. moCorTifNames in particular is
        % what the FISSA driver reads to build each group's image list, and
        % cellposeParams/roiParams are the provenance of a segmented set --
        % rewriting only three variables silently dropped them.
        S = load(roiFile);
        movefile(roiFile,...
            fullfile(roiDir,[animal '_OLDmoCorrROI_' filters{ROIfileN} '.mat']))
        S.moCorROI           = matchNroiS{ROIfileN};
        S.nTifs              = nTifs;
        S.tifIDXinAllTifList = tifIDXinAllTifList;
        save(roiFile,'-struct','S');
        clear S roiFile
    end
end