function results = runTests(varargin)
%runTests  Run the matlabPAC_process2P test suite.
%
%   runTests                         every suite (unit, then integration)
%   runTests('unit')                 only tests that need no recorded data
%   runTests('integration')          only tests that read the lab data drive
%   runTests(...,'Filter','FRA')     only tests whose name contains 'FRA'
%   runTests(...,'Isolate',true)     one MATLAB process per test
%   runTests(...,'Verbose',true)     stream each test's output as it runs
%   results = runTests(...)          struct array: .name .suite .status .secs .msg .out
%
%   From the shell:
%       matlab -batch "addpath('tests'); runTests"
%       matlab -batch "addpath('tests'); runTests('unit')"
%
%   Suites
%   ------
%     unit/            synthetic fixtures only; runs on any machine with the
%                      repo, in seconds. Break one of these and the cause is
%                      in the repo, not in the data.
%     integration/     reads recorded datasets under testConfig().dataRoot. A
%                      test whose fixture is missing reports SKIP, not FAIL,
%                      so the suite still means something without that drive.
%     investigations/  NOT run here. Diagnostic scripts with no assertions,
%                      kept for provenance; see tests/README.md.
%
%   Each test runs inside a disposable function workspace, so a test written
%   as a script cannot clobber the runner's own variables -- which is exactly
%   what goes wrong when tests are driven from a loop in a script. Figures are
%   forced invisible and closed between tests. 'Isolate' spends a MATLAB
%   startup per test (~20 s) to buy full state isolation, which is worth it
%   when chasing a test that only fails in company.
%
%   These tests print a lot. By default each one's output is captured and the
%   report stays one line per test; a failing test's captured output is
%   replayed in full underneath the summary, which is when you actually want
%   it. 'Verbose' streams instead, for watching a long test make progress.
%
%   Exit status: with no output requested the function errors at the end if
%   anything FAILed, so `matlab -batch` returns non-zero for CI.
%
%   See also testConfig, requireFixture, testSandbox.

p = inputParser;
% The suite validator has to REJECT non-suite strings, not just accept any
% char. With a permissive validator inputParser greedily binds the optional
% positional to the first argument, so runTests('Filter','FRA') consumed
% 'Filter' as the suite and then complained that 'FRA' had no value.
addOptional(p,'suite','all',@isSuiteArg);
addParameter(p,'Filter','',@(x) ischar(x)||isstring(x));
addParameter(p,'Isolate',false,@islogical);
addParameter(p,'Verbose',false,@islogical);
parse(p,varargin{:});
suite   = lower(char(p.Results.suite));
filt    = char(p.Results.Filter);
isolate = p.Results.Isolate;
verbose = p.Results.Verbose;

cfg = testConfig();
switch suite
    case 'unit',        suites = {'unit'};
    case 'integration', suites = {'integration'};
    otherwise,          suites = {'unit','integration'};   % 'all'
end

%% ---- collect ----
files = {}; suiteOf = {};
for s = 1:numel(suites)
    d = dir(fullfile(cfg.testsRoot,suites{s},'*.m'));
    for i = 1:numel(d)
        if ~isempty(filt) && ~contains(d(i).name,filt,'IgnoreCase',true); continue; end
        files{end+1}   = fullfile(d(i).folder,d(i).name); %#ok<AGROW>
        suiteOf{end+1} = suites{s};                       %#ok<AGROW>
    end
end
if isempty(files)
    error('runTests:noTests','No tests matched (suite=%s, filter=%s).',suite,filt);
end

fprintf('\n=== matlabPAC_process2P test suite ===\n');
fprintf('repo      %s\n', cfg.repoRoot);
if isfolder(cfg.dataRoot)
    fprintf('data root %s\n', cfg.dataRoot);
else
    fprintf('data root %s  (MISSING -- integration tests will skip)\n', cfg.dataRoot);
end
fprintf('%d test(s) across: %s\n\n', numel(files), strjoin(unique(suiteOf),', '));

vis0 = get(0,'DefaultFigureVisible');
restoreVis = onCleanup(@() set(0,'DefaultFigureVisible',vis0)); %#ok<NASGU>
set(0,'DefaultFigureVisible','off');

%% ---- run ----
results = struct('name',{},'suite',{},'status',{},'secs',{},'msg',{},'out',{},'file',{});
for i = 1:numel(files)
    [~,name] = fileparts(files{i});
    if verbose; fprintf('---- %s\n', name); end
    t0 = tic;
    if isolate
        [status,msg,out] = runIsolated(files{i}, cfg);
    else
        [status,msg,out] = runInProcess(files{i}, verbose);
    end
    secs = toc(t0);
    close all force
    results(end+1) = struct('name',name,'suite',suiteOf{i},'status',status, ...
        'secs',secs,'msg',msg,'out',out,'file',files{i}); %#ok<AGROW>
    fprintf('  %-4s  %-38s %6.1fs', status, name, secs);
    if ~strcmp(status,'PASS'); fprintf('  %s', firstLine(msg)); end
    fprintf('\n');
end

%% ---- summary ----
nPass = sum(strcmp({results.status},'PASS'));
nSkip = sum(strcmp({results.status},'SKIP'));
nFail = sum(strcmp({results.status},'FAIL'));
fprintf('\n=== %d passed, %d skipped, %d failed (%.0f s) ===\n', ...
    nPass, nSkip, nFail, sum([results.secs]));
if nSkip
    fprintf('skipped (fixture not on this machine):\n');
    for r = results(strcmp({results.status},'SKIP'))
        fprintf('  %-38s %s\n', r.name, firstLine(r.msg));
    end
end
if nFail
    for r = results(strcmp({results.status},'FAIL'))
        fprintf('\n--------------- FAILED: %s ---------------\n', r.name);
        fprintf('%s\n', r.msg);
        if ~verbose && ~isempty(strtrim(r.out))
            fprintf('---- output ----\n%s\n', r.out);
        end
        fprintf('file: %s\n', r.file);
    end
    if nargout == 0
        error('runTests:failures','%d test(s) failed.', nFail);
    end
end
if nargout == 0; clear results; end
end

%% ================= argument validation =================
function tf = isSuiteArg(x)
% True for a suite name, false for one of the name-value NAMES (so
% inputParser falls through and binds it as a parameter instead), and a clear
% error for anything else -- which is a mistyped suite.
tf = false;
if ~(ischar(x) || isstring(x)); return; end
x = char(x);
if any(strcmpi(x,{'all','unit','integration'})); tf = true; return; end
if any(strcmpi(x,{'Filter','Isolate','Verbose'})); return; end
error('runTests:unknownSuite', ...
    'Unknown suite "%s". Use ''unit'', ''integration'' or ''all''.', x);
end

%% ================= runners =================
function [status,msg,out] = runInProcess(f, verbose)
% The test executes inside THIS function's workspace, so a script-style test
% can only clobber the locals below, never the caller's bookkeeping.
out = ''; caught = [];
if verbose
    caught = tryInvoke(f);
else
    % The try/catch has to live INSIDE the evalc'd expression. evalc discards
    % what it captured if the expression throws -- which is the one case where
    % the output is worth keeping -- so tryInvoke returns the exception rather
    % than raising it.
    out = evalc('caught = tryInvoke(f);');
end
[status,msg] = classify(caught);
end

function ME = tryInvoke(f)
ME = [];
try
    invokeTest(f);
catch ME %#ok<CTCH,NASGU>
end
end

function [status,msg] = classify(ME)
if isempty(ME)
    status = 'PASS'; msg = '';
elseif strcmp(ME.identifier,'matlabPAC:test:fixtureMissing')
    status = 'SKIP'; msg = ME.message;
else
    status = 'FAIL'; msg = formatErr(ME);
end
end

function invokeTest(f)
[~,name] = fileparts(f);
isFcn = true;
try
    nargin(name);      % errors for a script file
catch
    isFcn = false;
end
if isFcn
    feval(name);
else
    run(f);
end
end

function [status,msg,out] = runIsolated(f, cfg)
% Drive one test in its own MATLAB. The driver is written to a temp .m file
% rather than embedded in -batch, so no quoting survives a trip through the
% shell.
[~,name] = fileparts(f);
drv = [tempname '.m'];
fid = fopen(drv,'w');
fprintf(fid,'addpath(''%s'');\n', cfg.testsRoot);
fprintf(fid,'testConfig();\n');
fprintf(fid,'set(0,''DefaultFigureVisible'',''off'');\n');
fprintf(fid,'try\n');
fprintf(fid,'  isF = true; try, nargin(''%s''); catch, isF = false; end\n', name);
fprintf(fid,'  if isF, %s; else, run(''%s''); end\n', name, f);
fprintf(fid,'  disp(''__TEST_OK__'');\n');
fprintf(fid,'catch ME\n');
fprintf(fid,'  if strcmp(ME.identifier,''matlabPAC:test:fixtureMissing'')\n');
fprintf(fid,'    fprintf(''__TEST_SKIP__%%s\\n'',ME.message);\n');
fprintf(fid,'  else\n');
fprintf(fid,'    fprintf(''__TEST_FAIL__%%s\\n'',ME.message);\n');
fprintf(fid,'  end\n');
fprintf(fid,'end\n');
fclose(fid);
c = onCleanup(@() delete(drv)); %#ok<NASGU>

exe = fullfile(matlabroot,'bin','matlab');
[~,out] = system(sprintf('"%s" -batch "run(''%s'')" 2>&1', exe, drv));

if contains(out,'__TEST_OK__')
    status = 'PASS'; msg = '';
elseif contains(out,'__TEST_SKIP__')
    status = 'SKIP'; msg = afterTag(out,'__TEST_SKIP__');
else
    status = 'FAIL'; msg = afterTag(out,'__TEST_FAIL__');
    if isempty(strtrim(msg)); msg = firstLine(out); end
end
end

%% ================= helpers =================
function s = afterTag(out,tag)
k = strfind(out,tag);
if isempty(k); s = ''; return; end
s = firstLine(out(k(1)+numel(tag):end));
end

function s = firstLine(s)
s = strtrim(s);
nl = find(s==newline,1,'first');
if ~isempty(nl); s = strtrim(s(1:nl-1)); end
end

function s = formatErr(ME)
s = ME.message;
if ~isempty(ME.stack)
    s = sprintf('%s  [%s:%d]', s, ME.stack(1).name, ME.stack(1).line);
end
end
