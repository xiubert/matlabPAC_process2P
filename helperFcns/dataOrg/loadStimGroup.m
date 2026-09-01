function [T,info,rep] = loadStimGroup(src,varargin)
% loadStimGroup  Load a group file, validate it, and check its provenance.
%
%   [T,info,rep] = loadStimGroup(groupFile)
%   [...] = loadStimGroup(groupFile,'family','CGC','strict',true)
%
%   The entry point every group plotter should use instead of a bare `load`.
%   Loads the table, runs validateStimGroup on it, and cross-checks the
%   embedded manifest against the manifest .json on disk.
%
%   Inputs:
%     src - path to a <Family>_Group<g>.mat file.
%
%   Name/Value:
%     'family'  - 'BPN' | 'CGC'. Inferred from the filename by default.
%     'strict'  - error rather than warn when validation fails. Default true.
%     'verbose' - print the validation summary. Default false.
%
%   Outputs:
%     T    - the group table.
%     info - groupInfo struct, or [] for a legacy unstamped file.
%     rep  - validateStimGroup report.
%
%   Divergence between the embedded manifest and the on-disk .json means the
%   group file was built from a manifest that has since been edited, i.e. the
%   file no longer reflects its stated membership. That is exactly the drift
%   behind BPN_GroupB = {TO0009,TO0010,TO0011} while CGC_GroupB additionally
%   contained TO0012, with nothing recording which was intended. It warns
%   rather than errors, since the group file itself is still internally valid.
%
%   See also aggregateStimGroup, validateStimGroup, stimGroupSpec.

p = inputParser;
addRequired(p,'src',@(x) ischar(x)||isstring(x));
addParameter(p,'family','',@(x) ischar(x)||isstring(x));
addParameter(p,'strict',true,@islogical);
addParameter(p,'verbose',false,@islogical);
parse(p,src,varargin{:});
family  = char(p.Results.family);
strict  = p.Results.strict;
verbose = p.Results.verbose;

src = char(src);
if ~isfile(src)
    error('loadStimGroup:notFound','Group file not found: %s', src);
end

S = load(src);
if isempty(family)
    rep0 = validateStimGroup(src,'verbose',false);
    family = rep0.family;
end
spec = stimGroupSpec(family);

if ~isfield(S, spec.varname)
    error('loadStimGroup:missingVar','%s has no variable "%s".', src, spec.varname);
end
T = S.(spec.varname);

if isfield(S,'groupInfo'); info = S.groupInfo; else; info = []; end

rep = validateStimGroup(T, family, 'verbose',verbose);
% validateStimGroup only detects a stamp when handed a path; it was handed the
% table, so carry the authoritative answer over from the load above.
rep.hasGroupInfo = ~isempty(info);
if ~rep.ok
    msg = sprintf('%s failed validation:\n%s', src, ...
        formatProblems(rep.problems(rep.problems.severity=="error",:)));
    if strict
        error('loadStimGroup:invalid','%s', msg);
    else
        warning('loadStimGroup:invalid','%s', msg);
    end
end

%% ---- manifest cross-check ----
if ~isempty(info) && isfield(info,'manifest')
    [d,nm] = fileparts(src);
    jsonPath = fullfile(d,[nm '_manifest.json']);
    if isfile(jsonPath)
        try
            onDisk = jsondecode(fileread(jsonPath));
            a1 = sort(string(getfielddef(info.manifest,'animals',strings(0,1))));
            a2 = sort(string(getfielddef(onDisk,'animals',strings(0,1))));
            if ~isequal(a1(:),a2(:))
                warning('loadStimGroup:manifestDrift', ...
                    ['%s: embedded manifest animals {%s} differ from %s {%s}. ' ...
                     'The group file predates an edit to the manifest; ' ...
                     're-run aggregateStimGroup to match.'], ...
                    nm, strjoin(a1,','), [nm '_manifest.json'], strjoin(a2,','));
            end
        catch ME
            warning('loadStimGroup:manifestUnreadable', ...
                'Could not compare against %s: %s', jsonPath, ME.message);
        end
    end
end

%% ---- convention staleness ----
if ~isempty(info) && isfield(info,'convention')
    if ~isequaln(info.convention, spec.convention)
        warning('loadStimGroup:staleConvention', ...
            ['%s was built under a different analysis convention than ' ...
             'stimGroupSpec(''%s'') currently defines. Re-process the animals ' ...
             'and re-aggregate, or the group is not comparable with others.'], ...
            nm, family);
    end
end
end

%% ---- helpers ----
function v = getfielddef(s,f,d)
if isstruct(s) && isfield(s,f); v = s.(f); else; v = d; end
end

function s = formatProblems(P)
s = '';
for i = 1:height(P)
    s = [s sprintf('  [%s] %s: %s\n', P.severity(i), P.code(i), P.message(i))]; %#ok<AGROW>
end
end
