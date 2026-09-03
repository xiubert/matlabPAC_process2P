function T = listRuns(dataPath,opts)
% listRuns  What analyses exist for an animal, and where.
%
%   T = listRuns(dataPath)
%   listRuns(dataPath)              % prints the table
%
%   Answers "which analyses have been run on this animal?" without opening a
%   single .mat file. Reads the runInfo.json each run folder carries, and
%   falls back to inspecting the folder for runs that predate it (or were
%   written by hand).
%
%   Inputs
%     dataPath  animal data folder.
%
%   Name-value
%     'verbose'  print the table when nargout is 0. Default true.
%
%   Output (table, newest first)
%     run        folder name under analysis/
%     roiSource  'handDrawn' | 'cellpose' | '?' -- read from runInfo.json, or
%                inferred from the ROI file when absent
%     created    when the run was first written
%     nROI       ROI count, when the ROI file records one
%     families   which stim families have a processed table
%     path       the run folder
%
%   A flat animal folder (artifacts loose, no analysis/) returns an empty
%   table -- there are no named runs to list, which is itself the answer.
%
%   See also animalPaths, writeRunInfo, migrateAnimalArtifacts

arguments
    dataPath (1,:) char
    opts.verbose (1,1) logical = true
end

if ~isfolder(dataPath)
    error('listRuns:noDataPath','Not a folder: %s',dataPath);
end
animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');

d = dir(fullfile(dataPath,'analysis','*'));
d = d([d.isdir] & ~ismember({d.name},{'.','..'}));

run = strings(0,1); roiSource = strings(0,1); created = strings(0,1);
nROI = nan(0,1); families = strings(0,1); pth = strings(0,1);

for k = 1:numel(d)
    A = fullfile(d(k).folder,d(k).name);
    rec = struct();
    j = fullfile(A,'runInfo.json');
    if isfile(j)
        try, rec = jsondecode(fileread(j)); catch, end
    end

    run(end+1,1)     = string(d(k).name);                       %#ok<AGROW>
    created(end+1,1) = string(getdef(rec,'created', ...
        datestr(d(k).datenum,'yyyy-mm-dd HH:MM:SS'))); %#ok<AGROW,DATST>
    pth(end+1,1)     = string(A);                               %#ok<AGROW>

    [src,n] = roiFacts(A,animal);
    roiSource(end+1,1) = string(getdef(rec,'roiSource',src));   %#ok<AGROW>
    nROI(end+1,1)      = n;                                     %#ok<AGROW>
    families(end+1,1)  = string(strjoin(familiesPresent(A,animal),','));  %#ok<AGROW>
end

T = table(run,roiSource,created,nROI,families,pth);
if ~isempty(T)
    [~,ord] = sort(T.created,'descend');
    T = T(ord,:);
end

if nargout == 0 && opts.verbose
    if isempty(T)
        fprintf(['No analysis/<run>/ folders in %s.\n' ...
            'Either nothing has been run, or this folder still uses the flat\n' ...
            'layout -- see migrateAnimalArtifacts.\n'],dataPath);
    else
        fprintf('Analyses of %s:\n',animal);
        disp(T(:,{'run','roiSource','created','nROI','families'}));
    end
    clear T
end
end

% ------------------------------------------------------------------------
function [src,n] = roiFacts(A,animal)
%read the ROI set directly when runInfo.json is absent or silent
src = '?'; n = NaN;
f = dir(fullfile(A,[animal '_moCorrROI_*.mat']));
if isempty(f); return, end
try
    S = load(fullfile(f(1).folder,f(1).name));
    if isfield(S,'moCorROI'); n = numel(S.moCorROI); end
    %report what the ROIs ARE, not which function wrote them
    if isfield(S,'cellposeParams') || (isfield(S,'roiParams') && ...
            isfield(S.roiParams,'source') && contains(S.roiParams.source,'cellpose'))
        src = 'cellpose';
    else
        src = 'handDrawn';
    end
catch
end
end
% ------------------------------------------------------------------------
function fam = familiesPresent(A,animal)
fam = {};
for f = stimGroupSpec()
    spec = stimGroupSpec(f{1});
    if isfile(fullfile(A,[animal spec.suffixProcessed]))
        fam{end+1} = f{1}; %#ok<AGROW>
    end
end
end
% ------------------------------------------------------------------------
function v = getdef(s,f,d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
