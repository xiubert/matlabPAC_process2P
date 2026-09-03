function FRAout = runFRA(dataPath,opts)
% runFRA  Non-interactive Frequency Response Area mapping for one animal.
%
%   FRAout = runFRA(dataPath)
%   FRAout = runFRA(dataPath,'pkPTsigSD',2,'showPlots',false)
%
%   A function-shaped way to run one animal's FRA: load its tifFileList, call
%   FRAmap, save <animal>_FRAmap.mat, and write the per-animal table via
%   FRAmap2table. The analysis is identical to the processFRA script; this
%   just takes arguments and returns a value instead of reading the caller's
%   workspace and plotting.
%
%   For a whole-animal run you normally do NOT need this:
%   processAnimalStimFamilies drives FRA along with every other family, and
%   processAnimal2Pheadless stage 11 goes through that. Use runFRA to
%   (re)compute FRA on its own, e.g. with a different pkPTsigSD.
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
%                        FRAmap2table -- the file aggregateStimGroup reads.
%                        Default true.
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

%tifFileList.map(n).folder is whatever path the inventory was built on, often
%a drive letter from another machine or a folder that has since moved. FRAmap
%resolves _Pulses.mat relative to it, so repoint it at dataPath -- where the
%tif list itself lives, alongside the tifs and their pulse files -- whenever
%the stored path is not present here. A stored path that does exist is left
%alone. Same rule as processFRA.
if ~isfolder(tifFileList.map(1).folder)
    [tifFileList.map.folder] = deal(dataPath);
end

%FRAmap is both a function and, by convention downstream, the variable name
%of its output; keep them apart here so the call is unambiguous inside a
%function workspace
FRAout = FRAmap(tifFileList,opts.pkPTsigSD,opts.nFramesPostPulse,opts.FsourceString);

%same two variables processFRA saves, under the same names
out.dataPath = dataPath;
out.FRAmap   = FRAout;
save(outFile,'-struct','out','-v7.3');
%uSigPkResp is nDB x nFreq, not per cell -- the per-cell array is
%CellSigPkLinDBfreq, whose rows are cells
if isfield(FRAout,'CellSigPkLinDBfreq')
    nCell = size(FRAout.CellSigPkLinDBfreq,1);
else
    nCell = NaN;
end
fprintf('runFRA: wrote %s (%d cells)\n',outFile,nCell);

if opts.writeTable
    if exist('FRAmap2table','file') ~= 2
        warning('runFRA:noTableFcn',...
            ['FRAmap2table is not on the path, so no per-animal FRA table '...
             'was written. The _FRAmap.mat is complete.']);
    else
        try
            T = FRAmap2table(FRAout,animal); %#ok<NASGU>
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
