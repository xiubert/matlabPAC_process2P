function testProcessAnimalStimFamilies()
% testProcessAnimalStimFamilies  processAnimal2P section 11: auto-run the
%                                per-stim process* scripts.
%
%   matlab -batch "addpath('tests'); runTests('Filter','testProcessAnimalStimFamilies')"
%
% Builds an animal folder holding _raw tables for both families (as
% stimParam2ROI would leave it), then checks that the processed tables appear,
% validate, and match what running the scripts by hand produces.

addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();
d = requireFixture(cfg.aggregateDir,'TOMT aggregate group files');
nPass = 0; nFail = 0;

[tmp,tmpCleanup] = testSandbox('processAnimalStimFamilies'); %#ok<ASGLU>
aDir = fullfile(tmp,'TO0007'); mkdir(aDir);

%% ---- lay out the animal folder as stimParam2ROI leaves it ----
for fam = {'BPN','CGC'}
    spec = stimGroupSpec(fam{1});
    G = getfield(load(fullfile(d,sprintf(spec.groupPattern,'D')), ...
        'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>
    Traw = G(string(G.animal)=="TO0007",:);
    Traw = removevars(Traw, intersect(spec.derivedVars,Traw.Properties.VariableNames));
    S = struct('anmlROIbyStim',Traw, ...
        'stimTable',unique(Traw(:,spec.stimVars)), ...
        'tifStimParamTable',table(1),'dataPath',aDir); %#ok<NASGU>
    save(fullfile(aDir,['TO0007' spec.suffixRaw]),'-struct','S','-v7.3');
end
check(numel(dir(fullfile(aDir,'*_raw.mat')))==2, 'two _raw tables in place, no processed ones');

%% ---- stage the FRA inputs too ----
% FRA has no _raw stage: it is driven from the tif inventory. processFRA reads
% the F traces out of tifFileList itself and the stim params out of the
% _Pulses.mat files beside it, so the .tif files are not needed and the
% fixture stays small (~13 MB + 200 KB rather than gigabytes of tifs).
% Without this the folder is an animal that genuinely has no FRA data, and the
% honest outcome is a skip -- which is asserted separately below.
srcAnimal = fullfile(cfg.animalsRoot,'TO0007');
haveFRAsrc = isfile(fullfile(srcAnimal,'TO0007_tifFileList.mat'));
if haveFRAsrc
    copyfile(fullfile(srcAnimal,'TO0007_tifFileList.mat'),aDir);
    pf = dir(fullfile(srcAnimal,'*_Pulses.mat'));
    for k = 1:numel(pf); copyfile(fullfile(pf(k).folder,pf(k).name),aDir); end
end

%% ---- run ----
% processAnimalStimFamilies reports one record per family in stimGroupSpec,
% including families this animal has no data for. Only the _raw-driven ones
% actually run, so split the report before asserting on it.
% plotAllROI off: this animal has enough ROIs that the per-ROI FRA tiles would
% open a stack of figures for no benefit here
res = processAnimalStimFamilies(aDir,'verbose',false, ...
    'scriptVars',struct('plotAllROI',false));
ran = res(~[res.skipped]);
expectRan = {'BPN','CGC'};
if haveFRAsrc; expectRan = {'BPN','CGC','FRA'}; end
check(isequal(sort({ran.family}),expectRan), ...
    sprintf('families that ran: %s (expected %s)', ...
        strjoin({ran.family},', '), strjoin(expectRan,', ')));
check(all([ran.ok]), 'all families processed ok');
check(all(arrayfun(@(r) isfile(r.procFile), ran)), 'processed tables written');
check(all(arrayfun(@(r) isfile(r.rawFile), ran)), 'inputs still present');

% FRA has no _raw stage (suffixRaw is ''), so its gate is the tif inventory
% rather than a _raw table. This used to be impossible to satisfy: the gate
% built <dataPath>/<animal>, which never exists, so FRA was skipped for every
% animal and the reason claimed it had no tifs of this family.
fra = res(strcmp({res.family},'FRA'));
check(isscalar(fra) && endsWith(fra.rawFile,'_tifFileList.mat'), ...
    'FRA gate checks the tif inventory, not a nonexistent _raw path');
if haveFRAsrc
    check(fra.ok && ~fra.skipped, ...
        sprintf('FRA runs when the inventory is present (reason: "%s")',fra.reason));
    check(endsWith(fra.procFile,'_anmlROI_FRAtable.mat') && isfile(fra.procFile), ...
        'FRA wrote its long-form table');
else
    check(fra.skipped && contains(fra.reason,'no tif inventory'), ...
        'FRA skips with an accurate reason when the inventory is absent');
end

% ...and it must still skip, with an ACCURATE reason, when there is nothing to
% run on. The wrong reason is what made the old behaviour costly: an animal
% with map tifs was told it had none.
sbx = testSandbox('noMapTifs');
noMap = fullfile(sbx,'AA0001'); mkdir(noMap);
sNo = processAnimalStimFamilies(noMap,'families',{'FRA'},'verbose',false);
check(sNo.skipped && contains(sNo.reason,'no tif inventory'), ...
    sprintf('missing inventory reported as such: "%s"',sNo.reason));

tifFileList = struct('stim',struct('name','x.tif'));   % inventory, but no .map
save(fullfile(noMap,'AA0001_tifFileList.mat'),'tifFileList');
sNo2 = processAnimalStimFamilies(noMap,'families',{'FRA'},'verbose',false);
check(sNo2.skipped && contains(sNo2.reason,'no map tifs'), ...
    sprintf('inventory without map tifs reported as such: "%s"',sNo2.reason));

for r = ran
    vn = stimGroupSpec(r.family).varname;   % FRA's is anmlROIbyFRA, not ...byStim
    T = getfield(load(r.procFile,vn),vn); %#ok<GFLD>
    rep = validateStimGroup(T,r.family,'verbose',false);
    check(rep.ok, sprintf('%s processed table validates',r.family));
end

% values must match the repaired group data
Tc = getfield(load(fullfile(aDir,'TO0007_anmlROI_CGCstimTable.mat'), ...
    'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>
Gc = getfield(load(fullfile(d,'CGC_GroupD.mat'),'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>
Gc = Gc(string(Gc.animal)=="TO0007",:);
dv = max(cellfun(@(a,b) max(abs(a-b),[],'all'), Tc.dFF_PT_avg, Gc.dFF_PT_avg));
check(dv < 1e-12, sprintf('CGC dFF_PT_avg matches the group data (max dev %.3g)',dv));

%% ---- the group step now succeeds without a manual step ----
outDir = fullfile(tmp,'out'); mkdir(outDir);
man = struct('group','auto','family','CGC','animals',"TO0007", ...
    'cohortRoot',tmp,'outDir',outDir);
gi = aggregateStimGroup(man,'verbose',false);
check(gi.nCells==27, ...
    sprintf('aggregateStimGroup runs straight after: %d cells',gi.nCells));

%% ---- options ----
% no figures left behind by default
nFigBefore = numel(findall(0,'Type','figure'));
processAnimalStimFamilies(aDir,'families',{'BPN'},'verbose',false);
check(numel(findall(0,'Type','figure'))==nFigBefore, ...
    'showPlots=false leaves no figures open');

% overwrite=false skips an already-processed family
res2 = processAnimalStimFamilies(aDir,'overwrite',false,'verbose',false);
done = res2(ismember({res2.family},{'BPN','CGC'}));
check(all([res2.skipped]) && all(contains({done.reason},'already exists')), ...
    'overwrite=false skips families that are already processed');

% a family with no _raw table is skipped, not an error
onlyBPN = fullfile(tmp,'TO0099'); mkdir(onlyBPN);
copyfile(fullfile(aDir,'TO0007_anmlROI_BPNstimTable_raw.mat'), ...
    fullfile(onlyBPN,'TO0099_anmlROI_BPNstimTable_raw.mat'));
res3 = processAnimalStimFamilies(onlyBPN,'verbose',false);
cgc = res3(strcmp({res3.family},'CGC'));
check(cgc.skipped && contains(cgc.reason,'no _raw table'), ...
    'a family this animal does not have is skipped with a reason');

% one failing family must not stop the other
bad = fullfile(tmp,'TO0055'); mkdir(bad);
copyfile(fullfile(aDir,'TO0007_anmlROI_BPNstimTable_raw.mat'), ...
    fullfile(bad,'TO0055_anmlROI_BPNstimTable_raw.mat'));
Sb = struct('anmlROIbyStim',table(1),'stimTable',table(1)); %#ok<NASGU>
save(fullfile(bad,'TO0055_anmlROI_CGCstimTable_raw.mat'),'-struct','Sb','-v7.3');
ws = warning('off','processAnimalStimFamilies:familyFailed');
res4 = processAnimalStimFamilies(bad,'verbose',false);
warning(ws);
okBPN = res4(strcmp({res4.family},'BPN')).ok;
badCGC = res4(strcmp({res4.family},'CGC'));
check(okBPN && ~badCGC.ok && ~isempty(badCGC.reason), ...
    'a failing family is reported while the other still processes');

% bad path
threw = false;
try, processAnimalStimFamilies(fullfile(tmp,'nope'),'verbose',false); catch, threw = true; end
check(threw,'a non-existent data folder raises');


%% ============ summary ============
fprintf('\n===============================\n');
fprintf('  %d passed, %d failed\n', nPass, nFail);
fprintf('===============================\n');
if nFail > 0; error('testProcessAnimalStimFamilies: %d check(s) failed', nFail); end

    function check(cond,label)
        if cond; nPass = nPass+1; fprintf('  PASS  %s\n',label);
        else;    nFail = nFail+1; fprintf('  FAIL  %s\n',label); end
    end
end
