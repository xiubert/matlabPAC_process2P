function FRAout = runFRA(dataPath,opts)
% runFRA  Non-interactive Frequency Response Area mapping for one animal.
%
%   FRAout = runFRA(dataPath)
%   FRAout = runFRA(dataPath,'pkPTsigSD',2,'showPlots',false)
%
%   The headless counterpart to the processFRA script. processFRA resolves
%   dataPath through uigetdir and plots as it goes, and it cannot be driven
%   by processAnimalStimFamilies -- FRA is registered in stimGroupSpec with an
%   empty suffixRaw, so that dispatcher's isfile gate can never fire. This
%   does the same work as a function: load the animal's tifFileList, call
%   FRAmap, save <animal>_FRAmap.mat, and (optionally) the per-animal table
%   FRAmap2table writes next to it.
%
%   The analysis is unchanged from processFRA: significance is tested once
%   per (ROI, freq, dB) on the trial-averaged onset-aligned trace against a
%   strictly pre-onset baseline.
%
%   Inputs
%     dataPath - animal data folder holding <animal>_tifFileList.mat with a
%                .map field (the FRA/BF mapping tifs).
%
%   Name-value
%     'animal'           default: [A-Z]{2}\d{4} from dataPath.
%     'pkPTsigSD'        significance threshold in baseline SDs. Default 2.
%     'nFramesPostPulse' frames added to the peak window. Default 2.
%     'FsourceString'    trace field to use. Default 'SCALEDfissaFroi'.
%     'showPlots'        draw plotFRAmap. Default false.
%     'writeTable'       also write <animal>_anmlROI_FRAtable.mat via
%                        FRAmap2table, when that function is on the path (it
%                        lives on the refactorFRA branch). Default true.
%     'overwrite'        redo when <animal>_FRAmap.mat exists. Default true.
%
%   Output
%     FRAout  the FRAmap struct that was saved.
%
%   See also processFRA, FRAmap, FRAmap2table, processAnimal2Pheadless,
%   plotFRAmap

arguments
    dataPath (1,:) char
    opts.animal           (1,:) char = ''
    opts.pkPTsigSD        (1,1) double = 2
    opts.nFramesPostPulse (1,1) double = 2
    opts.FsourceString    (1,:) char = 'SCALEDfissaFroi'
    opts.showPlots        (1,1) logical = false
    opts.writeTable       (1,1) logical = true
    opts.overwrite        (1,1) logical = true
end

if ~isfolder(dataPath)
    error('runFRA:noDataPath','Not a folder: %s',dataPath);
end
animal = opts.animal;
if isempty(animal)
    animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');
end
if isempty(animal)
    error('runFRA:noAnimal','Could not derive an animal ID from %s',dataPath);
end

outFile = fullfile(dataPath,[animal '_FRAmap.mat']);
if isfile(outFile) && ~opts.overwrite
    S = load(outFile,'FRAmap');
    FRAout = S.FRAmap;
    fprintf('runFRA: %s already exists; not overwriting.\n',outFile);
    return
end

listFile = fullfile(dataPath,[animal '_tifFileList.mat']);
if ~isfile(listFile)
    error('runFRA:noTifFileList','Not found: %s',listFile);
end
S = load(listFile,'tifFileList');
if ~isfield(S,'tifFileList') || ~isfield(S.tifFileList,'map') || isempty(S.tifFileList.map)
    error('runFRA:noMapTifs',...
        ['%s has no .map tifs, so there is no FRA to compute. Mark the mapping '...
         'tifs in stage 1 (mapTifs) if this animal has them.'],listFile);
end
tifFileList = S.tifFileList;

%FRAmap is both a function and, by convention downstream, the variable name
%of its output; keep them apart here so the call is unambiguous inside a
%function workspace
FRAout = FRAmap(tifFileList,opts.pkPTsigSD,opts.nFramesPostPulse,opts.FsourceString);

%same two variables processFRA saves, under the same names
out.dataPath = dataPath;
out.FRAmap   = FRAout;
save(outFile,'-struct','out','-v7.3');
fprintf('runFRA: wrote %s (%d cells)\n',outFile,size(FRAout.uSigPkResp,1));

if opts.writeTable
    if exist('FRAmap2table','file') ~= 2
        warning('runFRA:noTableFcn',...
            ['FRAmap2table is not on the path (it lives on the refactorFRA '...
             'branch), so no per-animal FRA table was written. The '...
             '_FRAmap.mat is complete.']);
    else
        try
            T = FRAmap2table(FRAout); %#ok<NASGU>
            tblFile = fullfile(dataPath,[animal '_anmlROI_FRAtable.mat']);
            save(tblFile,'T','dataPath','-v7.3');
            fprintf('runFRA: wrote %s\n',tblFile);
        catch ME
            warning('runFRA:tableFailed',...
                'FRAmap2table failed (%s); the _FRAmap.mat is still valid.',ME.message);
        end
    end
end

if opts.showPlots
    plotFRAmap(FRAout,'plotAllROI',true);
end
end
