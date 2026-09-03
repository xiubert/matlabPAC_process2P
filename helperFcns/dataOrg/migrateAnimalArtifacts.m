function moved = migrateAnimalArtifacts(dataPath,opts)
% migrateAnimalArtifacts  Move a flat animal folder's artifacts into analysis/<run>/.
%
%   moved = migrateAnimalArtifacts(dataPath)                 % dry run
%   moved = migrateAnimalArtifacts(dataPath,'run','legacy','apply',true)
%
%   Optional. The pipeline reads a flat folder perfectly well, so nothing has
%   to be migrated -- this exists for when you want one layout everywhere
%   rather than two.
%
%   Moves only the PER-RUN artifacts (see animalPaths): the ROI files, the
%   FISSA output, _moCorr_Tifs_Params, _tifFileList, _pulseLegend2P,
%   _stimGroupIDX, the _anmlROI_* tables, _FRAmap, and the QC images. It does
%   NOT touch the raw acquisition, the two legends, or NoRMCorred/*.tif --
%   those are shared by every analysis of the animal and stay where they are.
%
%   Inputs
%     dataPath  animal data folder.
%
%   Name-value
%     'run'     folder name under analysis/. Default 'legacy'.
%     'apply'   actually move. Default FALSE -- the call reports what it would
%               do and changes nothing, so you can look before committing.
%     'animal'  animal ID. Default: derived from dataPath.
%
%   Output
%     moved     table of what was (or would be) moved: name, from, to.
%
%   Back up the folder first. This moves files rather than copying them, so
%   the flat layout is gone afterwards; anything still expecting it (an old
%   script with a hardcoded path, say) will stop finding its inputs.
%
%   See also animalPaths, processAnimal2Pheadless

arguments
    dataPath (1,:) char
    opts.run    (1,:) char = 'legacy'
    opts.apply  (1,1) logical = false
    opts.animal (1,:) char = ''
end

animal = opts.animal;
if isempty(animal); animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once'); end
if isempty(animal)
    error('migrateAnimalArtifacts:noAnimal', ...
        'Could not derive an animal ID from %s; pass ''animal''.',dataPath);
end

P = animalPaths(dataPath,'run',opts.run);
if strcmp(P.artifacts,dataPath)
    error('migrateAnimalArtifacts:sameFolder', ...
        'Resolved artifact folder is the animal folder itself; nothing to do.');
end

pat = {[animal '_moCorrROI_*.mat'],[animal '_OLDmoCorrROI_*.mat'], ...
       [animal '_moCorr_Tifs_Params.mat'],[animal '_tifFileList.mat'], ...
       [animal '_anmlROI_*.mat'],[animal '_FRAmap.mat'], ...
       [animal '_pulseLegend2P.mat'],[animal '_stimGroupIDX.mat'], ...
       [animal '_ROIoverlay_*.png'],[animal '_cellposeMean.tif']};

names = strings(0,1); froms = strings(0,1); tos = strings(0,1);
for k = 1:numel(pat)
    d = dir(fullfile(dataPath,pat{k}));
    for j = 1:numel(d)
        names(end+1,1) = string(d(j).name);                       %#ok<AGROW>
        froms(end+1,1) = string(fullfile(d(j).folder,d(j).name)); %#ok<AGROW>
        tos(end+1,1)   = string(fullfile(P.artifacts,d(j).name)); %#ok<AGROW>
    end
end

%the cellposeMean written by older runs sits under NoRMCorred/
oldMean = fullfile(P.moCorrDir,[animal '_cellposeMean.tif']);
if isfile(oldMean)
    names(end+1,1) = string([animal '_cellposeMean.tif']);
    froms(end+1,1) = string(oldMean);
    tos(end+1,1)   = string(fullfile(P.artifacts,[animal '_cellposeMean.tif']));
end

legacyFissa = fullfile(P.moCorrDir,'FISSAoutput');
hasFissa = isfolder(legacyFissa);
if hasFissa
    names(end+1,1) = "FISSAoutput/";
    froms(end+1,1) = string(legacyFissa);
    tos(end+1,1)   = string(P.fissaDir);
end

moved = table(names,froms,tos,'VariableNames',{'name','from','to'});

if isempty(names)
    fprintf('migrateAnimalArtifacts: nothing to move in %s\n',dataPath);
    return
end

if ~opts.apply
    fprintf(['migrateAnimalArtifacts: DRY RUN -- %d item(s) would move to %s\n' ...
             'Re-run with ''apply'',true to do it.\n'],numel(names),P.artifacts);
    disp(moved(:,{'name'}))
    return
end

if ~isfolder(P.artifacts); mkdir(P.artifacts); end
for k = 1:numel(names)
    if names(k) == "FISSAoutput/"
        if isfolder(P.fissaDir); rmdir(P.fissaDir,'s'); end
        movefile(char(froms(k)),char(tos(k)));
    else
        movefile(char(froms(k)),char(tos(k)));
    end
end
fprintf('migrateAnimalArtifacts: moved %d item(s) to %s\n',numel(names),P.artifacts);
end
