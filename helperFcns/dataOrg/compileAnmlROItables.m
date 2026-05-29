function T = compileAnmlROItables(family,varargin)
% compileAnmlROItables  Concatenate every animal's anmlROIbyStim table
%                       under a cohort directory into one long table,
%                       for a single stim family.
%
%   T = compileAnmlROItables(family)               % uigetdir prompt
%   T = compileAnmlROItables(family, parentPath)   % char path
%   T = compileAnmlROItables(family, params)       % struct: parentPath, tableDir
%
%   Scans parentPath for animal-ID subfolders (pattern '[A-Z]{2}\d{4}',
%   e.g. AA0067), loads each animal's per-stim-family table, applies a
%   universal column reduction, and vertically concatenates the result
%   across animals into a single cohort-level table.
%
%   Inputs:
%     family       - one of:
%                      'CGC'        (PT-in-contrast)
%                      'BPN'        (band-pass noise)
%                      'Spont'      (spontaneous)
%                      'dContrast'  (contrast change)
%                    Validated against the lookup below; unknown values
%                    raise an error listing the supported set.
%     varargin{1}  - (optional) either:
%                      char folder path (uses default tableDir = '.')
%                    or
%                      params struct with fields:
%                        parentPath - cohort directory
%                        tableDir   - subdirectory under each animal
%                                     folder where the .mat lives.
%                                     Configurable; defaults to '.',
%                                     i.e. the animal folder itself.
%
%   Family -> (filename suffix, .mat variable name) lookup:
%     'CGC'        '_anmlROI_CGCstimTable.mat'      anmlROIbyStim
%     'BPN'        '_anmlROI_BPNstimTable.mat'      anmlROIbyStim
%     'Spont'      '_anmlROI_SpontstimTable.mat'    anmlROIbyStim
%     'dContrast'  '_anmlROI_dContrastTable.mat'    anmlROIdContrast
%   Source: the four save(...) blocks in stimParam2ROI.m. 'dContrast'
%   tables save their table under a different variable name; the load
%   step is family-aware so callers don't have to be.
%
%   Expected per-animal file layout:
%     <parentPath>/<animal>/<tableDir>/<animal><filename suffix>
%
%   Per-animal column reduction (applied before concatenation):
%     - Adds  F  = SCALEDfissaFroi (canonical fluorescence column;
%                  assigned in a clearly delimited block below — change
%                  the source there if needed for non-FISSA data).
%     - Drops every column whose name contains 'Froi'.
%     - Drops {'Pulse','stimID'} (per-animal-local identifiers). Drop is
%       defensive (intersect with existing columns first).
%     - Family-specific stim parameter columns are PRESERVED (e.g.
%       BPN/Spont per-pulse params, CGC dBrange/PTampl/DRC columns).
%
%   Output:
%     T - long-form table: one row per (animal, roiID, unique-stim) for
%         the selected family, with stim params + frameRate + F
%         (canonical fluorescence). T.animal is cast to string;
%         T.roiID is cast to numeric.
%
%   Notes:
%     - F is hard-coded to SCALEDfissaFroi; tables without that column
%       will error (intentional — non-FISSA data should be flagged, not
%       silently filled).
%     - Animals whose <tableDir> or per-family .mat file is missing are
%       skipped with a console notice.

% --- family -> file/varname lookup ---
familyMap = struct(...
    'CGC',       struct('suffix','_anmlROI_CGCstimTable.mat',   'varname','anmlROIbyStim'), ...
    'BPN',       struct('suffix','_anmlROI_BPNstimTable.mat',   'varname','anmlROIbyStim'), ...
    'Spont',     struct('suffix','_anmlROI_SpontstimTable.mat', 'varname','anmlROIbyStim'), ...
    'dContrast', struct('suffix','_anmlROI_dContrastTable.mat', 'varname','anmlROIdContrast'));

if ~isfield(familyMap, family)
    error('compileAnmlROItables:UnknownFamily', ...
        'Unknown family "%s". Supported: %s', family, strjoin(fieldnames(familyMap),', '));
end
suffix  = familyMap.(family).suffix;
varname = familyMap.(family).varname;

% --- resolve parentPath / tableDir ---
defParams.tableDir = '.';
switch numel(varargin)
    case 0
        pDir = uigetdir('D:\Data\','Select cohort directory');
        params = defParams;
        params.parentPath = pDir;
    case 1
        if ischar(varargin{1}) && isfolder(varargin{1})
            params = defParams;
            params.parentPath = varargin{1};
        elseif isstruct(varargin{1})
            params = varargin{1};
            if ~isfield(params,'tableDir')
                params.tableDir = defParams.tableDir;
            end
        else
            error('compileAnmlROItables:BadArg', ...
                'Second argument must be a folder path (char) or a params struct.');
        end
end

% --- list animal-ID subfolders ---
aDir = dir(params.parentPath);
aDir = aDir([aDir.isdir]);
aDir = aDir(cellfun(@any, regexp({aDir.name}', '[A-Z]{2}\d{4}')));

% --- per-animal load + reduce ---
perAnimal = cell(numel(aDir),1);
for aNum = 1:numel(aDir)
    animal = aDir(aNum).name;
    disp(['Compiling stim/response data for: ' animal]);

    animalDir = fullfile(params.parentPath, animal, params.tableDir);
    if ~isfolder(animalDir)
        disp([params.tableDir ' directory does not exist for ' animal ...
              ' — not loaded into table'])
        continue
    end

    matFile = fullfile(animalDir, [animal suffix]);
    if ~isfile(matFile)
        disp([animal suffix ' missing in ' animalDir ' — not loaded into table'])
        continue
    end

    aData = load(matFile, varname);
    Ta = aData.(varname);

    %=== Canonical fluorescence column ====================================
    % F defaults to SCALEDfissaFroi (FISSA-scaled, motion-corrected dF/F
    % source used downstream). To swap sources for a non-FISSA dataset or
    % to A/B a different signal, change the right-hand side here. Any
    % column whose name contains 'Froi' is dropped in the next step, so
    % the assignment must happen before that.
    Ta.F = Ta.SCALEDfissaFroi;
    %======================================================================

    % drop *Froi columns (raw, motion-corrected, FISSA, scaled-FISSA)
    froiVars = Ta.Properties.VariableNames(contains(Ta.Properties.VariableNames,'Froi'));
    if ~isempty(froiVars)
        Ta = removevars(Ta, froiVars);
    end

    % drop per-animal-local identifiers (defensive intersection)
    dropVars = intersect({'Pulse','stimID'}, Ta.Properties.VariableNames);
    if ~isempty(dropVars)
        Ta = removevars(Ta, dropVars);
    end

    perAnimal{aNum} = Ta;
end

% --- concatenate ---
perAnimal = perAnimal(~cellfun(@isempty, perAnimal));
if isempty(perAnimal)
    T = table();
    return
end
T = vertcat(perAnimal{:});
T.animal = string(T.animal);
T.roiID  = str2double(T.roiID);
end
