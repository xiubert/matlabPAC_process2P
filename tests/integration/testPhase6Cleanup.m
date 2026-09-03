function testPhase6Cleanup()
% testPhase6Cleanup  Phase 6: manifests, retired combinetable, compileCohortData fix.
%
%   matlab -batch "addpath('tests'); runTests('Filter','testPhase6Cleanup')"

addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();
d = requireFixture(cfg.aggregateDir,'TOMT aggregate group files');
nPass = 0; nFail = 0;

%% ============ 1. manifests ============
fprintf('\n=== 1. group manifests ===\n');
for fam = {'BPN','CGC'}
    for g = {'A','B','C','D'}
        jf = fullfile(d,sprintf('%s_Group%s_manifest.json',fam{1},g{1}));
        okFile = isfile(jf);
        okContent = false;
        if okFile
            m = jsondecode(fileread(jf));
            okContent = strcmp(m.family,fam{1}) && strcmp(m.group,g{1}) && ...
                ~isempty(m.animals) && isfield(m,'cohortRoot') && isfield(m,'outDir');
        end
        check(okFile && okContent, sprintf('%s_Group%s_manifest.json valid',fam{1},g{1}));
    end
end

% manifests must match what the group files actually contain
allMatch = true;
for fam = {'BPN','CGC'}
    spec = stimGroupSpec(fam{1});
    for g = {'A','B','C','D'}
        m = jsondecode(fileread(fullfile(d,sprintf('%s_Group%s_manifest.json',fam{1},g{1}))));
        T = getfield(load(fullfile(d,sprintf(spec.groupPattern,g{1})),spec.varname),spec.varname); %#ok<GFLD>
        allMatch = allMatch && isequal(sort(string(m.animals(:))), sort(unique(string(T.animal))));
    end
end
check(allMatch,'every manifest matches its group file''s actual animals');

% the BPN/CGC Group B membership difference is now recorded, not hidden
mb = jsondecode(fileread(fullfile(d,'BPN_GroupB_manifest.json')));
mc = jsondecode(fileread(fullfile(d,'CGC_GroupB_manifest.json')));
check(isequal(setxor(string(mb.animals),string(mc.animals)),"TO0012"), ...
    'Group B family difference (TO0012, CGC only) is visible in the manifests');

% a json path is accepted by the aggregator
threw = '';
try
    aggregateStimGroup(fullfile(d,'BPN_GroupD_manifest.json'),'dryRun',true,'verbose',false);
catch ME
    threw = ME.identifier;
end
check(strcmp(threw,'aggregateStimGroup:fileNotFound'), ...
    'json manifest parses; fails only at per-animal file resolution (not present on this machine)');

%% ============ 2. combinetable retired ============
fprintf('\n=== 2. combinetable ===\n');
ct = fullfile(d,'combinetable.m');
check(isfile(ct),'combinetable.m still present (as a signpost)');
threw = '';
try
    run(ct);
catch ME
    threw = ME.identifier;
end
check(strcmp(threw,'combinetable:retired'),'combinetable.m errors with migration instructions');
src = fileread(ct);
check(contains(src,'aggregateStimGroup') && contains(src,'_manifest.json'), ...
    'retirement notice names the replacement and the manifests');
check(isfile(fullfile(d,'bkmat','retired_20260831','combinetable.m')), ...
    'original combinetable.m preserved in bkmat/');

%% ============ 3. compileCohortData signature fix ============
fprintf('\n=== 3. compileAnmlROItables ===\n');
cc = fileread(fullfile(cfg.repoRoot,'stimulusSpecific','extraCGC','compileCohortData.m'));
check(contains(cc,'compileAnmlROItables(family,params)'), ...
    'compileCohortData passes family as the first argument');
check(contains(cc,"family     = 'CGC'"),'family is an editable parameter');
check(isfile(fullfile(cfg.repoRoot,'stimulusSpecific','extraCGC','plotCohortData.m')) && ...
      ~isfile(fullfile(cfg.repoRoot,'plotCohortData.m')), ...
    'cohort-table pair moved out of the repo root into stimulusSpecific/extraCGC');
pc = fileread(fullfile(cfg.repoRoot,'stimulusSpecific','extraCGC','plotCohortData.m'));
check(contains(pc,'dFF_DRC - mean(dFF_DRC') && ~contains(pc,'dFoFcalc(F,[find((t>=tBasePT'), ...
    'plotCohortData uses the additive dFF_PT, not the divisive form');
check(~contains(pc,'params.nFramesPostPulse,params.pkPTsigSD'), ...
    'both peak calls use params.pkPTframeBin');

% the old broken call shape must still be rejected, not silently accepted
threw = '';
try
    compileAnmlROItables(struct('parentPath',d,'tableDir','.'));
catch ME
    threw = ME.identifier;
end
check(~isempty(threw), sprintf('the old one-arg call still fails loudly (%s)',threw));

% the fixed shape resolves; no animal subfolders here, so it returns empty
T = compileAnmlROItables('CGC', struct('parentPath',d,'tableDir','.'));
check(istable(T), 'compileAnmlROItables(family,params) runs and returns a table');

check(any(strcmp('BPN',stimGroupSpec())) && any(strcmp('CGC',stimGroupSpec())), ...
    'stimGroupSpec() lists the validated families');

%% ============ summary ============
fprintf('\n===============================\n');
fprintf('  %d passed, %d failed\n', nPass, nFail);
fprintf('===============================\n');
if nFail > 0; error('testPhase6Cleanup: %d check(s) failed', nFail); end

    function check(cond,label)
        if cond; nPass = nPass+1; fprintf('  PASS  %s\n',label);
        else;    nFail = nFail+1; fprintf('  FAIL  %s\n',label); end
    end
end
