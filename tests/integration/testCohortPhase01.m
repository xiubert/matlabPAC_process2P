% testCohortPhase01  Phase 0 (stimGroupSpec, validateStimGroup) and
%                    Phase 1 (cohort primitives) against real TOMT data.
%
%   matlab -batch "addpath('tests'); runTests('Filter','testCohortPhase01')"
%
% Asserts:
%   1. validator is clean on the current (repaired) group files
%   2. validator reproduces the known diagnoses on the pre-fix bkmat backups
%   3. the degenerate-n contract from etc/harmonization_plan.md holds

function testCohortPhase01()

addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();
d = requireFixture(cfg.aggregateDir,'TOMT aggregate group files');
bkT = fullfile(d,'bkmat','preTrim_20260831');
bkR = fullfile(d,'bkmat','preRecompute_20260831');

nPass = 0; nFail = 0;

%% ============ 1. current files validate clean ============
fprintf('\n=== 1. current group files ===\n');
refBPN = []; refCGC = [];
for g = {'A','B','C','D'}
    r = validateStimGroup(fullfile(d,sprintf('BPN_Group%s.mat',g{1})),'verbose',true);
    check(r.ok, sprintf('BPN_Group%s ok',g{1}));
    if isempty(refBPN); refBPN = r; end
end
for g = {'A','B','C','D'}
    r = validateStimGroup(fullfile(d,sprintf('CGC_Group%s.mat',g{1})),'verbose',true);
    check(r.ok, sprintf('CGC_Group%s ok',g{1}));
    if isempty(refCGC); refCGC = r; end
end

% cross-group: schema + time axis must agree
TA = getfield(load(fullfile(d,'CGC_GroupA.mat'),'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>
refVarsCGC = TA.Properties.VariableNames;
refAxisCGC = TA.t_dFF_DRC{1};
allOK = true;
for g = {'B','C','D'}
    r = validateStimGroup(fullfile(d,sprintf('CGC_Group%s.mat',g{1})), ...
        'refVars',refVarsCGC,'refTimeAxis',refAxisCGC,'verbose',false);
    allOK = allOK && r.ok && ~hasCode(r,"columnOrderDiffers") && ~hasCode(r,"timeAxisDiffersFromRef");
end
check(allOK,'CGC B/C/D conform to GroupA schema + time axis');

TB = getfield(load(fullfile(d,'BPN_GroupA.mat'),'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>
refAxisBPN = TB.t_dFF{1};
allOK = true;
for g = {'B','C','D'}
    r = validateStimGroup(fullfile(d,sprintf('BPN_Group%s.mat',g{1})), ...
        'refVars',TB.Properties.VariableNames,'refTimeAxis',refAxisBPN,'verbose',false);
    allOK = allOK && r.ok && ~hasCode(r,"timeAxisDiffersFromRef");
end
check(allOK,'BPN B/C/D conform to GroupA schema + time axis');

%% ============ 2. validator catches the pre-fix defects ============
fprintf('\n=== 2. pre-fix backups (known defects) ===\n');
r = validateStimGroup(fullfile(bkT,'BPN_GroupA.mat'),'verbose',false);
check(~r.ok && hasCode(r,"raggedTimeAxis"), 'pre-trim BPN_GroupA: ragged time axis detected');
r = validateStimGroup(fullfile(bkT,'BPN_GroupD.mat'),'verbose',false);
check(~r.ok && hasCode(r,"raggedTimeAxis"), 'pre-trim BPN_GroupD: ragged time axis detected');
r = validateStimGroup(fullfile(bkT,'BPN_GroupB.mat'),'verbose',false);
check(r.ok, 'pre-trim BPN_GroupB: internally uniform, no error (correct)');

% GroupC was internally uniform at 20 frames -- only the cross-group check sees it
r = validateStimGroup(fullfile(bkT,'BPN_GroupC.mat'),'verbose',false);
check(r.ok, 'pre-trim BPN_GroupC: internally uniform, no error (correct)');
r = validateStimGroup(fullfile(bkT,'BPN_GroupC.mat'),'refTimeAxis',refAxisBPN,'verbose',false);
check(hasCode(r,"timeAxisDiffersFromRef"), 'pre-trim BPN_GroupC: cross-group axis mismatch caught');

r = validateStimGroup(fullfile(bkR,'CGC_GroupC.mat'),'verbose',false);
check(~r.ok && hasCode(r,"traceTimeMismatch"), 'pre-rebuild CGC_GroupC: dFF_PT width mismatch detected');
r = validateStimGroup(fullfile(bkR,'CGC_GroupD.mat'),'verbose',false);
check(~r.ok && hasCode(r,"traceTimeMismatch"), 'pre-rebuild CGC_GroupD: dFF_PT width mismatch detected');

r = validateStimGroup(fullfile(bkR,'CGC_GroupB.mat'),'refVars',refVarsCGC,'verbose',false);
check(r.ok && hasCode(r,"columnOrderDiffers"), 'pre-reorder CGC_GroupB: column order flagged, no error');

%% ============ 3. degenerate-n contract ============
fprintf('\n=== 3. degenerate-n contract ===\n');
B = getfield(load(fullfile(d,'BPN_GroupD.mat'),'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>
C = getfield(load(fullfile(d,'CGC_GroupD.mat'),'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>

% --- groupN ---
N = groupN(B);
check(N.nCells==51 && N.nAnimals==2, sprintf('groupN BPN_GroupD: %s',N.label));
check(height(N.perAnimal)==2 && sum(N.perAnimal.nCells)==51, 'groupN per-animal breakdown sums to total');
N1 = groupN(B, string(B.animal)=="TO0006");
check(N1.singleAnimal && N1.nCells==24, sprintf('groupN single animal: %s',N1.label));
N0 = groupN(B, false(height(B),1));
check(N0.empty && N0.nCells==0, 'groupN empty selection is flagged, not an error');

% --- cohortMeanSEM: the two small-n traps ---
M = [0.1 0.2 0.3; 0.2 0.3 0.4; 0.3 0.4 0.5];
o = cohortMeanSEM(M);
check(o.n==3 && o.showBand && isequal(size(o.sem),[1 3]), 'cohortMeanSEM n=3: band shown, sem 1xN');
o1 = cohortMeanSEM(M(1,:));
check(o1.n==1 && ~o1.showBand && all(isnan(o1.sem)) && isequal(size(o1.mean),[1 3]), ...
    'cohortMeanSEM n=1: sem all NaN (not zeros), band suppressed, mean stays 1xN');
check(isequal(size(SEMcalc(M(1,:))),[1 1]) && isequal(size(o1.sem),[1 3]), ...
    'cohortMeanSEM avoids the SEMcalc row-vector transpose (1x1 vs 1x3)');
o0 = cohortMeanSEM(zeros(0,3));
check(o0.empty && o0.n==0 && ~o0.showBand, 'cohortMeanSEM n=0: empty, band suppressed');
threw = false;
try, cohortMeanSEM({1,2}); catch, threw = true; end
check(threw,'cohortMeanSEM rejects non-numeric input');

% --- gatherCellTraces ---
[M1,cells1,i1] = gatherCellTraces(B, B.BPNdBAmpl==50, 'dFF_avg', 't_dFF');
check(i1.ok && size(M1,1)==51 && size(M1,2)==15 && height(cells1)==51, ...
    sprintf('gatherCellTraces BPN dB=50: %s', mat2str(size(M1))));
check(isequal(size(i1.t),[1 15]), 'gatherCellTraces returns the common time vector');
% per-trial column must collapse to one row per cell
[M2,~,i2] = gatherCellTraces(B, B.BPNdBAmpl==50, 'dFF', 't_dFF');
check(i2.ok && isequal(size(M2),size(M1)), 'gatherCellTraces averages per-trial dFF to one row per cell');
check(max(abs(M2-M1),[],'all') < 1e-12, 'per-trial mean matches the stored dFF_avg');

% the CGC single-animal case that used to crash splitapply
keyC = strcat(string(C.animal),'_',string(C.roiID));
sigC = cellfun(@(x) logical(x(1)), C.sigPk);
ukC  = unique(keyC);
validC = arrayfun(@(i) all(sigC(keyC==ukC(i))), (1:numel(ukC))');
mask6 = ismember(keyC, ukC(validC)) & string(C.animal)=="TO0006";
[M3,~,i3] = gatherCellTraces(C, mask6, 'dFF_PT_avg', 't_dFF_DRC');
check(~i3.ok && isempty(M3) && ~isempty(i3.reason), ...
    sprintf('gatherCellTraces empty selection returns reason, not an error ("%s")', i3.reason));
o3 = cohortMeanSEM(reshape(M3,0,0));
check(o3.empty, 'cohortMeanSEM handles the empty result downstream');

% single cell end-to-end
oneCell = keyC == ukC(find(validC,1)) & C.dBdelta==10;
[M4,~,i4] = gatherCellTraces(C, oneCell, 'dFF_PT_avg', 't_dFF_DRC');
o4 = cohortMeanSEM(M4);
check(i4.ok && o4.n==1 && ~o4.showBand && all(isnan(o4.sem)), ...
    'single cell: gathered ok, band suppressed, sem NaN');

% ragged input must raise rather than silently trim
Braw = getfield(load(fullfile(bkT,'BPN_GroupA.mat'),'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>
threw = false;
try, gatherCellTraces(Braw, Braw.BPNdBAmpl==50, 'dFF_avg', 't_dFF'); catch, threw = true; end
check(threw,'gatherCellTraces raises on ragged pre-fix data');

% --- cohortStat ---
rng(0);
x = randn(20,1); y = x + 0.5 + 0.3*randn(20,1);
s = cohortStat(x,y);
check(s.ok && ~isnan(s.p) && ~isempty(s.test), sprintf('cohortStat n=20: %s p=%.3g',s.test,s.p));
sT = cohortStat(x,y,'test','ttest');
sR = cohortStat(x,y,'test','rank');
check(strcmp(sT.test,'ttest') && strcmp(sR.test,'signrank'), 'cohortStat honours an explicit test choice');
s1 = cohortStat(x(1),y(1));
check(~s1.ok && isnan(s1.p) && contains(s1.reason,'below minN'), ...
    sprintf('cohortStat n=1 refuses: "%s"',s1.reason));
s0 = cohortStat([],[]);
check(~s0.ok && isnan(s0.p), 'cohortStat n=0 refuses');
s2 = cohortStat([x;NaN],[y;NaN]);
check(s2.ok && s2.n==20, 'cohortStat drops NaN pairs');
sU = cohortStat(randn(15,1),randn(12,1),'paired',false);
check(sU.ok && numel(sU.n)==2, sprintf('cohortStat unpaired: %s n=[%d %d]',sU.test,sU.n(1),sU.n(2)));

% --- annotateN stamps the caveats ---
f = figure('Visible','off'); ax = axes(f);
h = annotateN(ax,N1);
str = char(join(string(h.String),' '));
check(contains(str,'1 animal'), 'annotateN stamps the single-animal caveat');
h2 = annotateN(axes(f), groupN(B, false(height(B),1)));
check(contains(char(join(string(h2.String),' ')),'no cells'), 'annotateN stamps the empty-group caveat');
close(f);


%% ============ 3b. sigCellCounts ============
fprintf('\n=== 3b. sigCellCounts ===\n');
TC = getfield(load(fullfile(d,'CGC_GroupD.mat'),'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>
TB = getfield(load(fullfile(d,'BPN_GroupD.mat'),'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>

cC = sigCellCounts(TC,'CGC');
check(isequal(cC.labels,["Low contrast" "High contrast"]), ...
    'CGC labels name the contrasts');
check(cC.nTotal==51 && numel(cC.nSig)==2, sprintf('CGC: %d cells, 2 levels',cC.nTotal));

% must agree with the matrix plotCGCgroup builds independently
set(0,'DefaultFigureVisible','off');
oC = plotCGCgroup(TC,'plots',{},'verbose',false);
check(isequal(cC.nSig, sum(oC.sig,1)), 'CGC per-contrast counts match plotCGCgroup''s sig matrix');
check(cC.nSigAll == sum(oC.valid), ...
    sprintf('CGC nSigAll (%d) equals plotCGCgroup''s both-contrast count (%d)', ...
    cC.nSigAll, sum(oC.valid)));
close all

cB = sigCellCounts(TB,'BPN');
r  = tableRLF(TB,'nConsec',2);
check(isequal(cB.nSig, sum(r.sigAll,1)), 'BPN per-level counts match tableRLF''s sigAll');
check(isequal(cB.labels(1),"30 dB") && cB.nTotal==51, ...
    sprintf('BPN: %d cells, labels like "%s"',cB.nTotal,cB.labels(1)));
check(max(abs(cB.pctSig - 100*cB.nSig/cB.nTotal)) < 1e-12, 'percentages are consistent');
check(cB.nSigAny >= max(cB.nSig) && cB.nSigAll <= min(cB.nSig), ...
    sprintf('nSigAll (%d) <= min per-level <= max per-level <= nSigAny (%d)', ...
    cB.nSigAll, cB.nSigAny));

% degenerate inputs
cE = sigCellCounts(TB,'BPN','rowMask',false(height(TB),1));
check(cE.nTotal==0 && isempty(cE.nSig) && cE.nSigAll==0, 'empty selection gives zero counts');
keyB = strcat(string(TB.animal),'_',string(TB.roiID));
c1 = sigCellCounts(TB,'BPN','rowMask',keyB==keyB(1));
check(c1.nTotal==1 && all(c1.nSig<=1), 'single cell: counts are 0 or 1 per level');

% a cell missing a level cannot count as significant everywhere
kC = string(TC.animal) + "_" + string(TC.roiID);
TCdrop = TC(~(kC==kC(1) & TC.dBdelta==30),:);
cD = sigCellCounts(TCdrop,'CGC');
check(cD.nSigAll <= cC.nSigAll, ...
    'dropping one contrast for a cell cannot increase nSigAll');

%% ============ summary ============
fprintf('\n===============================\n');
fprintf('  %d passed, %d failed\n', nPass, nFail);
fprintf('===============================\n');
if nFail > 0; error('testCohortPhase01: %d check(s) failed', nFail); end

%% ---- nested helpers (share nPass/nFail with the parent workspace) ----
    function check(cond,label)
        if cond; nPass = nPass+1; fprintf('  PASS  %s\n',label);
        else;    nFail = nFail+1; fprintf('  FAIL  %s\n',label); end
    end

    function tf = hasCode(rep,c)
        tf = any(rep.problems.code == c);
    end

end
