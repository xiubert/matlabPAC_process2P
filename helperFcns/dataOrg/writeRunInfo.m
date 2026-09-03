function infoFile = writeRunInfo(P,info)
% writeRunInfo  Describe an analysis run, in its folder and in the animal log.
%
%   infoFile = writeRunInfo(P)
%   infoFile = writeRunInfo(P,info)
%
%   A run folder full of .mat files does not say what it is. This drops a
%   runInfo.json beside them and appends a line to <animal>_analysisRuns.log
%   in the animal folder, so "what analyses exist for this animal, and where
%   are they?" is answerable by looking rather than by loading MATLAB files.
%
%   The log lives at the ANIMAL level on purpose: that is the folder someone
%   opens first, and it is the one place a per-run folder cannot point to
%   itself from.
%
%   Inputs
%     P     an animalPaths struct.
%     info  optional struct of extra fields to record -- roiSource, stages,
%           families, nROI, parameters, whatever the caller knows. Anything
%           already set by this function is not overwritten.
%
%   Output
%     infoFile  path to the runInfo.json written ('' in the flat layout,
%               which has no run folder to describe).
%
%   See also animalPaths, listRuns, processAnimal2Pheadless

arguments
    P    (1,1) struct
    info (1,1) struct = struct()
end

infoFile = '';
if P.isFlat
    return      % nothing to describe: artifacts are loose in the animal folder
end
if ~isfolder(P.artifacts); mkdir(P.artifacts); end

rec = info;
rec.run       = P.run;
rec.animal    = regexp(P.root,'[A-Z]{2}\d{4}','match','once');
rec.animalDir = P.root;
rec.artifacts = P.artifacts;
rec.fissaDir  = P.fissaDir;
if ~isfield(rec,'created') || isempty(rec.created)
    rec.created = char(datetime('now','Format','uuuu-MM-dd HH:mm:ss'));
end
rec.updated = char(datetime('now','Format','uuuu-MM-dd HH:mm:ss'));
if ~isfield(rec,'repoSHA') || isempty(rec.repoSHA)
    rec.repoSHA = repoSHA();
end

infoFile = fullfile(P.artifacts,'runInfo.json');
% keep whatever an earlier stage recorded; a later stage should add to the
% description of a run, not replace it
if isfile(infoFile)
    try
        old = jsondecode(fileread(infoFile));
        f = fieldnames(old);
        for k = 1:numel(f)
            if ~isfield(rec,f{k}) || isempty(rec.(f{k}))
                rec.(f{k}) = old.(f{k});
            end
        end
        if isfield(old,'created') && ~isempty(old.created)
            rec.created = old.created;   % first write wins
        end
    catch
    end
end

fid = fopen(infoFile,'w');
if fid < 0
    warning('writeRunInfo:cannotWrite','Could not write %s',infoFile);
    infoFile = '';
    return
end
fprintf(fid,'%s',jsonencode(rec,'PrettyPrint',true));
fclose(fid);

appendAnimalLog(P,rec);
end

% ------------------------------------------------------------------------
function appendAnimalLog(P,rec)
%one line per run in the animal folder, so the top level of an animal says
%which analyses exist and where each one lives
animal = rec.animal;
if isempty(animal); animal = 'animal'; end
logFile = fullfile(P.root,[animal '_analysisRuns.log']);

rel = strrep(P.artifacts,[P.root filesep],'');
line = sprintf('%s  %-28s %-14s %s', rec.updated, rec.run, ...
    getdef(rec,'roiSource','?'), rel);

prev = '';
if isfile(logFile); prev = fileread(logFile); end
%one entry per run: replace this run's line rather than appending a duplicate
keep = {};
if ~isempty(prev)
    for l = strsplit(strtrim(prev),newline)
        t = strtrim(l{1});
        if isempty(t); continue, end
        if startsWith(t,'#'); continue, end   % the header is rewritten below
        if contains(l{1},sprintf(' %s ',rec.run)); continue, end
        keep{end+1} = l{1}; %#ok<AGROW>
    end
end
keep{end+1} = line;

fid = fopen(logFile,'w');
if fid < 0; return, end
fprintf(fid,'# analyses of %s -- one line per run; artifacts are under analysis/\n',animal);
fprintf(fid,'# updated              run                          roiSource      location\n');
fprintf(fid,'%s\n',keep{:});
fclose(fid);
end
% ------------------------------------------------------------------------
function v = getdef(s,f,d)
if isfield(s,f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
% ------------------------------------------------------------------------
function sha = repoSHA()
sha = '';
try
    here = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    [st,out] = system(sprintf('git -C "%s" rev-parse --short HEAD 2>/dev/null',here));
    if st == 0; sha = strtrim(out); end
catch
end
end
