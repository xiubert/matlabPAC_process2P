function rep = validateStimGroup(src,varargin)
% validateStimGroup  Check a group table against its family's contract.
%
%   rep = validateStimGroup(matPath)                 % family from filename
%   rep = validateStimGroup(matPath,family)
%   rep = validateStimGroup(T,family)
%   rep = validateStimGroup(...,'refVars',vars,'verbose',false)
%
%   Reports problems; never modifies anything. Run it on a group file before
%   plotting, and inside aggregateStimGroup on every load.
%
%   The checks are the generalised gate half of etc/fixBPNgroupTraceLength.m
%   and etc/fixCGCgroupPTconvention.m -- the conditions that, unchecked, gave
%   ragged BPN traces and two incompatible CGC conventions living under
%   identical filenames.
%
%   Inputs:
%     src    - path to a group .mat, or a table.
%     family - 'BPN' | 'CGC'. Inferred from the filename when src is a path.
%
%   Name/Value:
%     'refVars' - canonical column order to compare against. Order mismatch is
%                 a warning (it blocks table vertcat but corrupts nothing).
%     'verbose' - print a summary. Default true.
%
%   Output (struct):
%     .ok        true when no 'error'-severity problems were found
%     .family    resolved family
%     .source    path or '<table>'
%     .nRows .nCells .nAnimals
%     .hasGroupInfo  whether the file carried a provenance stamp
%     .problems  table: severity, code, message, nAffected
%
%   Severities:
%     error   - the table is not safely usable (ragged axes, duplicate cells,
%               missing required columns)
%     warning - usable but non-conforming (column order, mixed frame rates)
%     info    - expected variation worth surfacing (unstamped legacy file,
%               ragged raw F columns, which is normal since raw traces sit on
%               t_total, not the dF/F axis)
%
%   See also stimGroupSpec, aggregateStimGroup, groupN.

p = inputParser;
addRequired(p,'src');
% the optional positional must not swallow a parameter name
addOptional(p,'family','',@(x) (ischar(x)||isstring(x)) && ...
    ~any(strcmpi(char(x),{'refVars','refTimeAxis','verbose'})));
addParameter(p,'refVars',{},@(x) iscellstr(x)||isstring(x)||isempty(x)); %#ok<ISCLSTR>
addParameter(p,'refTimeAxis',[],@(x) isempty(x)||isnumeric(x));
addParameter(p,'verbose',true,@islogical);
parse(p,src,varargin{:});
family  = char(p.Results.family);
refVars = p.Results.refVars;
refT    = p.Results.refTimeAxis;
verbose = p.Results.verbose;
if ~isempty(refVars); refVars = cellstr(refVars); end

sev = strings(0,1); code = strings(0,1); msg = strings(0,1); naff = zeros(0,1);
    function add(s,c,m,n)
        sev(end+1,1)=s; code(end+1,1)=c; msg(end+1,1)=string(m); naff(end+1,1)=n;
    end

%% ---- resolve source ----
hasGroupInfo = false;
if istable(src)
    T = src; srcName = '<table>';
else
    srcName = char(src);
    if ~isfile(srcName)
        error('validateStimGroup:notFound','File not found: %s', srcName);
    end
    if isempty(family); family = inferFamily(srcName); end
    S = load(srcName);
    spec0 = stimGroupSpec(family);
    if ~isfield(S, spec0.varname)
        error('validateStimGroup:missingVar','%s has no variable "%s".', srcName, spec0.varname);
    end
    T = S.(spec0.varname);
    hasGroupInfo = isfield(S,'groupInfo');
end
if isempty(family)
    error('validateStimGroup:noFamily','family could not be inferred; pass it explicitly.');
end
spec = stimGroupSpec(family);
vars = T.Properties.VariableNames;

rep.family = family;
rep.source = srcName;
rep.hasGroupInfo = hasGroupInfo;
rep.nRows = height(T);

if ~hasGroupInfo && ~istable(src)
    add("info","noGroupInfo", ...
        "no groupInfo stamp (legacy file); convention/provenance unknown",0);
end

%% ---- 1. required columns ----
missing = setdiff(spec.requiredVars, vars);
if ~isempty(missing)
    add("error","missingVars", ...
        sprintf('missing required column(s): %s', strjoin(missing,', ')), numel(missing));
    % without the core columns the remaining checks are meaningless
    rep = finish(rep,T,spec,sev,code,msg,naff,verbose);
    return
end

%% ---- 2. cell identity ----
idv = [spec.idVars, spec.stimVars];
idv = idv(ismember(idv,vars));
K = T(:,idv);
[~,ia] = unique(K,'rows','stable');
nDup = height(T) - numel(ia);
if nDup > 0
    add("error","duplicateCellStim", ...
        sprintf('%d duplicate (%s) row(s); cell identity is ambiguous', ...
        nDup, strjoin(idv,', ')), nDup);
end

%% ---- 3./4. dF/F time axis ----
tw = cellfun(@numel, T.(spec.timeVar));
if numel(unique(tw)) > 1
    add("error","raggedTimeAxis", ...
        sprintf('%s widths differ across rows: %s', spec.timeVar, mat2str(unique(tw)')), ...
        sum(tw ~= mode(tw)));
else
    tAll = cell2mat(cellfun(@(c) reshape(c,1,[]), T.(spec.timeVar), 'uni',0));
    dev = max(abs(tAll - tAll(1,:)),[],'all');
    if dev > 1e-9
        add("error","timeAxisMismatch", ...
            sprintf('%s values differ across rows (max dev %.3g)', spec.timeVar, dev), ...
            sum(any(abs(tAll-tAll(1,:))>1e-9,2)));
    end
    % cross-group check: single-file validation cannot see that two groups sit
    % on different axes (pre-fix BPN_GroupC was internally uniform at 20
    % frames while GroupB was uniform at 15), so pooling callers pass refTimeAxis.
    if ~isempty(refT)
        thisT = tAll(1,:);
        if numel(thisT) ~= numel(refT)
            add("warning","timeAxisDiffersFromRef", ...
                sprintf('%s is %d frames, reference is %d; groups cannot be pooled', ...
                spec.timeVar, numel(thisT), numel(refT)), 1);
        elseif max(abs(thisT - reshape(refT,1,[]))) > 1e-9
            add("warning","timeAxisDiffersFromRef", ...
                sprintf('%s values differ from reference (max dev %.3g)', ...
                spec.timeVar, max(abs(thisT - reshape(refT,1,[])))), 1);
        end
    end
end

%% ---- 5. trace widths vs the time axis ----
for k = 1:numel(spec.traceVars)
    v = spec.traceVars{k};
    if ~ismember(v,vars); continue; end
    w = cellfun(@(c) size(c,2), T.(v));
    bad = w ~= tw;
    if any(bad)
        add("error","traceTimeMismatch", ...
            sprintf('%s width does not match %s on %d row(s) (%s vs %s)', ...
            v, spec.timeVar, sum(bad), mat2str(unique(w)'), mat2str(unique(tw)')), sum(bad));
    end
end

%% ---- 6. peak / significance columns are per-cell scalars ----
for k = 1:numel(spec.peakVars)
    v = spec.peakVars{k};
    if ~ismember(v,vars); continue; end
    col = T.(v);
    if iscell(col)
        n = cellfun(@numel, col);
        if strcmp(v,spec.peakVar) || strcmp(v,spec.sigVar)
            if any(n ~= 1)
                add("warning","nonScalarPeak", ...
                    sprintf('%s is not scalar on %d row(s) (sizes %s); expected one value per cell', ...
                    v, sum(n~=1), mat2str(unique(n)')), sum(n~=1));
            end
        end
    end
end

%% ---- 7. raw F columns (informational) ----
for k = 1:numel(spec.rawTraceVars)
    v = spec.rawTraceVars{k};
    if ~ismember(v,vars); continue; end
    w = cellfun(@(c) size(c,2), T.(v));
    if numel(unique(w)) > 1
        add("info","raggedRawTrace", ...
            sprintf('%s widths vary (%s) - expected, raw F is on the t_total axis', ...
            v, mat2str(unique(w)')), numel(unique(w)));
    end
end

%% ---- 8. frame rate ----
if ismember('frameRate',vars)
    fr = unique(T.frameRate);
    if numel(fr) > 1
        add("warning","mixedFrameRate", ...
            sprintf('mixed frame rates in one group: %s Hz', mat2str(fr')), numel(fr));
    end
end

%% ---- 9. column order ----
if ~isempty(refVars)
    if ~isempty(setxor(vars,refVars))
        add("warning","variableSetDiffers", ...
            sprintf('variable set differs from reference (extra: %s | missing: %s)', ...
            strjoin(setdiff(vars,refVars),', '), strjoin(setdiff(refVars,vars),', ')), 1);
    elseif ~isequal(vars,refVars)
        add("warning","columnOrderDiffers", ...
            sprintf('%d of %d column positions differ from reference; blocks table vertcat', ...
            sum(~strcmp(vars,refVars)), numel(refVars)), sum(~strcmp(vars,refVars)));
    end
end

rep = finish(rep,T,spec,sev,code,msg,naff,verbose);
end

%% ================= helpers =================
function rep = finish(rep,T,spec,sev,code,msg,naff,verbose)
rep.problems = table(sev,code,msg,naff, ...
    'VariableNames',{'severity','code','message','nAffected'});
rep.ok = ~any(sev == "error");

idv = spec.idVars(ismember(spec.idVars, T.Properties.VariableNames));
if numel(idv) == numel(spec.idVars)
    N = groupN(T,'idVars',idv);
    rep.nCells = N.nCells; rep.nAnimals = N.nAnimals;
else
    rep.nCells = NaN; rep.nAnimals = NaN;
end

if verbose
    [~,nm,ex] = fileparts(rep.source);
    if isempty(nm); label = rep.source; else; label = [nm ex]; end
    fprintf('%-24s %-4s | %4d rows, %3d cells, %2d animals | %s\n', ...
        label, rep.family, rep.nRows, rep.nCells, rep.nAnimals, ...
        statusWord(rep.ok, sev));
    for i = 1:height(rep.problems)
        fprintf('    [%-7s] %-22s %s\n', rep.problems.severity(i), ...
            rep.problems.code(i), rep.problems.message(i));
    end
end
end

function w = statusWord(ok,sev)
if ~ok
    w = sprintf('FAIL (%d error)', sum(sev=="error"));
elseif any(sev=="warning")
    w = sprintf('ok, %d warning(s)', sum(sev=="warning"));
else
    w = 'ok';
end
end

function family = inferFamily(pathStr)
[~,nm] = fileparts(pathStr);
fams = stimGroupSpec();
hit = fams(cellfun(@(f) startsWith(nm,[f '_']), fams));
if isempty(hit)
    error('validateStimGroup:cannotInferFamily', ...
        'Cannot infer family from "%s"; pass it explicitly. Known: %s', nm, strjoin(fams,', '));
end
family = hit{1};
end
