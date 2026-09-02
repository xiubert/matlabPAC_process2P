function testAggregateStimGroup()
% testAggregateStimGroup  Phase 2: aggregator round-trip + provenance.
%
%   matlab -batch "addpath('tests'); runTests('Filter','testAggregateStimGroup')"
%
% Acceptance criterion from etc/harmonization_plan.md: regenerating a group
% must reproduce the CURRENT table exactly (isequaln), adding only groupInfo.
%
% The real per-animal TOMT tables are not on this machine (only the aggregate
% directory is), so each current group file is split back into per-animal
% tables in the expected on-disk layout, then re-aggregated and compared.
% That exercises resolve -> load -> validate -> canonicalize -> concat ->
% stamp against real data. It does NOT exercise path resolution against the
% researcher's actual cohort tree.

addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();
d = requireFixture(cfg.aggregateDir,'TOMT aggregate group files');
[tmp,tmpCleanup] = testSandbox('aggregateStimGroup'); %#ok<ASGLU>
cohortRoot = fullfile(tmp,'cohort'); outDir = fullfile(tmp,'out');
mkdir(cohortRoot); mkdir(outDir);

nPass = 0; nFail = 0;

%% ============ 1. round-trip every group ============
fprintf('\n=== 1. round-trip: aggregate reproduces the current table ===\n');
for fam = {'BPN','CGC'}
    spec = stimGroupSpec(fam{1});
    for g = {'A','B','C','D'}
        f = fullfile(d,sprintf(spec.groupPattern,g{1}));
        T0 = getfield(load(f,spec.varname),spec.varname); %#ok<GFLD>

        animals = unique(string(T0.animal),'stable');   % preserve row order
        for a = 1:numel(animals)
            aDir = fullfile(cohortRoot,char(animals(a)));
            if ~isfolder(aDir); mkdir(aDir); end
            S = struct(spec.varname, T0(string(T0.animal)==animals(a),:));
            save(fullfile(aDir,[char(animals(a)) spec.suffixProcessed]),'-struct','S','-v7.3');
        end

        man = struct('group',g{1},'family',fam{1},'animals',animals, ...
            'cohortRoot',cohortRoot,'outDir',outDir);
        gi = aggregateStimGroup(man,'verbose',false);

        T1 = getfield(load(gi.outFile,spec.varname),spec.varname); %#ok<GFLD>
        check(isequaln(T0,T1), sprintf('%s_Group%s round-trip is exact (isequaln)',fam{1},g{1}));
        check(gi.nCells==height(unique(T0(:,{'animal','roiID'}))) && ...
              gi.nAnimals==numel(animals), ...
              sprintf('%s_Group%s groupInfo counts: %d cells / %d animals', ...
              fam{1},g{1},gi.nCells,gi.nAnimals));
    end
end

%% ============ 2. provenance stamp ============
fprintf('\n=== 2. provenance ===\n');
spec = stimGroupSpec('CGC');
gf = fullfile(outDir,sprintf(spec.groupPattern,'D'));
[T,info,rep] = loadStimGroup(gf);
check(~isempty(info) && info.schemaVersion==1, 'groupInfo present with schemaVersion');
check(isequaln(info.convention, spec.convention), 'convention stamped from stimGroupSpec');
check(isequal(sort(info.animals(:)),sort(["TO0006";"TO0007"])), 'animals recorded');
check(height(info.sourceFiles)==2 && all(info.sourceFiles.bytes>0), 'source files recorded');
check(height(info.perAnimal)==2, 'per-animal breakdown recorded');
check(isequal(numel(info.timeAxis), numel(T.(spec.timeVar){1})), 'time axis stamped');
check(rep.ok && rep.hasGroupInfo, 'loadStimGroup validates and sees the stamp');
check(isfile(fullfile(outDir,'CGC_GroupD_manifest.json')), 'manifest .json written alongside');

%% ============ 3. manifest drift is caught ============
fprintf('\n=== 3. manifest drift ===\n');
jf = fullfile(outDir,'CGC_GroupD_manifest.json');
m = jsondecode(fileread(jf));
m.animals = [string(m.animals(:)); "TO0099"];      % edit the manifest after the fact
fid = fopen(jf,'w'); fprintf(fid,'%s',jsonencode(m,'PrettyPrint',true)); fclose(fid);
lastwarn('');
[~,~,~] = loadStimGroup(gf);
[wmsg,wid] = lastwarn;
check(strcmp(wid,'loadStimGroup:manifestDrift'), ...
    sprintf('drift warned: %s', firstLine(wmsg)));

%% ============ 4. single-animal group ============
fprintf('\n=== 4. single-animal group ===\n');
man1 = struct('group','solo','family','CGC','animals',"TO0007", ...
    'cohortRoot',cohortRoot,'outDir',outDir);
gi1 = aggregateStimGroup(man1,'verbose',false);
check(gi1.nAnimals==1 && gi1.nCells==27, ...
    sprintf('single-animal group aggregates: %d cells / %d animal',gi1.nCells,gi1.nAnimals));
[T1,~,r1] = loadStimGroup(gi1.outFile);
N1 = groupN(T1);
check(r1.ok && N1.singleAnimal, 'single-animal group validates and is flagged by groupN');

%% ============ 5. explicit file list ============
fprintf('\n=== 5. explicit files ===\n');
manF = struct('group','files','family','BPN','outDir',outDir, ...
    'files',{{fullfile(cohortRoot,'TO0006','TO0006_anmlROI_BPNstimTable.mat'), ...
              fullfile(cohortRoot,'TO0007','TO0007_anmlROI_BPNstimTable.mat')}});
giF = aggregateStimGroup(manF,'verbose',false);
check(giF.nAnimals==2 && isequal(sort(giF.animals(:)),sort(["TO0006";"TO0007"])), ...
    'files-only manifest resolves animals from filenames');

%% ============ 6. failure modes ============
fprintf('\n=== 6. failure modes ===\n');
% missing animal file
manX = struct('group','X','family','BPN','animals',"TO9999", ...
    'cohortRoot',cohortRoot,'outDir',outDir);
check(throws(@() aggregateStimGroup(manX,'verbose',false)), 'missing per-animal file raises');

% an invalid per-animal table must block aggregation
badDir = fullfile(cohortRoot,'TO0006');
good = fullfile(badDir,'TO0006_anmlROI_BPNstimTable.mat');
bak = [good '.bak']; copyfile(good,bak);
Sb = load(good,'anmlROIbyStim'); Tb = Sb.anmlROIbyStim;
Tb.dFF_avg{1} = Tb.dFF_avg{1}(1:5);                 % ragged: 5 frames vs t_dFF 15
Sbad = struct('anmlROIbyStim',Tb); %#ok<NASGU>
save(good,'-struct','Sbad','-v7.3');
manB = struct('group','bad','family','BPN','animals',["TO0006","TO0007"], ...
    'cohortRoot',cohortRoot,'outDir',outDir);
check(throws(@() aggregateStimGroup(manB,'verbose',false)), ...
    'invalid per-animal table blocks aggregation');
movefile(bak,good);


% the common mistake: processAnimal2P run, process* NOT run
rawOnly = fullfile(cohortRoot,'TO0006');
procPath = fullfile(rawOnly,'TO0006_anmlROI_CGCstimTable.mat');
specC = stimGroupSpec('CGC');
Tp = getfield(load(procPath,'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>
Traw = removevars(Tp, intersect(specC.derivedVars,Tp.Properties.VariableNames));
Sraw = struct('anmlROIbyStim',Traw); %#ok<NASGU>
save(fullfile(rawOnly,'TO0006_anmlROI_CGCstimTable_raw.mat'),'-struct','Sraw','-v7.3');
movefile(procPath,[procPath '.hidden']);
manNP = struct('group','np','family','CGC','animals',"TO0006", ...
    'cohortRoot',cohortRoot,'outDir',outDir);
errId = ''; errMsg = '';
try
    aggregateStimGroup(manNP,'verbose',false);
catch ME
    errId = ME.identifier; errMsg = ME.message;
end
check(strcmp(errId,'aggregateStimGroup:notProcessed'), ...
    'unprocessed animal (only _raw present) is refused, not silently aggregated');
check(contains(errMsg,'processCGC'), ...
    'the error names the script to run (processCGC)');
movefile([procPath '.hidden'],procPath);
delete(fullfile(rawOnly,'TO0006_anmlROI_CGCstimTable_raw.mat'));

% dryRun writes nothing
outSolo = fullfile(outDir,'BPN_Groupdry.mat');
if isfile(outSolo); delete(outSolo); end
manD = struct('group','dry','family','BPN','animals',["TO0006","TO0007"], ...
    'cohortRoot',cohortRoot,'outDir',outDir);
aggregateStimGroup(manD,'dryRun',true,'verbose',false);
check(~isfile(outSolo), 'dryRun writes nothing');

% legacy file (no groupInfo) loads with info = []
[~,infoL,repL] = loadStimGroup(fullfile(d,'BPN_GroupA.mat'));
check(isempty(infoL) && repL.ok && ~repL.hasGroupInfo, ...
    'legacy unstamped group file still loads and validates');


%% ============ summary ============
fprintf('\n===============================\n');
fprintf('  %d passed, %d failed\n', nPass, nFail);
fprintf('===============================\n');
if nFail > 0; error('testAggregateStimGroup: %d check(s) failed', nFail); end

%% ---- nested helpers ----
    function check(cond,label)
        if cond; nPass = nPass+1; fprintf('  PASS  %s\n',label);
        else;    nFail = nFail+1; fprintf('  FAIL  %s\n',label); end
    end
end

function tf = throws(fn)
tf = false;
try, fn(); catch, tf = true; end
end

function s = firstLine(str)
c = splitlines(string(str)); s = char(c(1));
end
