function testPlotCGCgroup()
% testPlotCGCgroup  Phase 4: CGC group plotter.
%
%   matlab -batch "addpath('tests'); runTests('Filter','testPlotCGCgroup')"
%
% The headline case is TO0006: a single-animal CGC group with ZERO cells
% significant in both contrasts, which crashed the old processCGC path with
% "Group numbers must be a vector of positive integers" from splitapply.

addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();
d = requireFixture(cfg.aggregateDir,'TOMT aggregate group files');
set(0,'DefaultFigureVisible','off'); c0 = onCleanup(@() set(0,'DefaultFigureVisible','on'));

nPass = 0; nFail = 0;

T = getfield(load(fullfile(d,'CGC_GroupD.mat'),'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>

%% ============ 1. full group ============
fprintf('\n=== 1. full group (CGC_GroupD) ===\n');
o = plotCGCgroup(fullfile(d,'CGC_GroupD.mat'),'verbose',false);
check(o.N.nCells==51 && o.N.nAnimals==2, sprintf('counts: %s',o.N.label));
check(isequal(o.levels,[10 30]), 'contrasts ascending: low (10) then high (30)');
check(sum(o.valid)==15, sprintf('both-contrast significant cells: %d (matches the repaired data)',sum(o.valid)));
check(isequal(size(o.pk),[51 2]) && isequal(size(o.sig),[51 2]), 'pk/sig are nCells x nContrast');
check(all([o.traces.n]==15), 'trace panels use the 15 significant cells');
check(o.stat.ok && o.stat.n==15, sprintf('paired test ran: %s p=%.4g n=%d',o.stat.test,o.stat.p,o.stat.n));
close all

% cross-check the peak matrix against a direct lookup
key = strcat(string(T.animal),'_',string(T.roiID));
uk  = unique(key,'stable');
pkRef = nan(numel(uk),2);
for i = 1:numel(uk)
    for k = 1:2
        s = key==uk(i) & T.dBdelta==o.levels(k);
        if any(s); pkRef(i,k) = T.pkPT{find(s,1)}(1); end
    end
end
check(isequaln(o.pk,pkRef), 'pk matrix matches a direct roiID+dBdelta lookup');

% cross-check against the regenerated _sigpk export
Sp = load(fullfile(d,'CGC_GroupD_sigpk.mat'),'Tpop');
check(height(Sp.Tpop)==sum(o.valid)*2, ...
    sprintf('valid count agrees with CGC_GroupD_sigpk.mat (%d rows = %d cells x 2)', ...
    height(Sp.Tpop),sum(o.valid)));

%% ============ 2. THE crash case: single animal, zero cells ============
fprintf('\n=== 2. single animal with zero significant cells (TO0006) ===\n');
T6 = T(string(T.animal)=="TO0006",:);
o6 = plotCGCgroup(T6,'verbose',false);
check(o6.N.singleAnimal && o6.N.nCells==24, sprintf('loaded: %s',o6.N.label));
check(sum(o6.valid)==0, 'zero cells significant in both contrasts');
check(all([o6.traces.n]==0), 'trace panels report n=0 instead of crashing');
check(all(cellfun(@(m) all(isnan(m)), {o6.traces.mean})), 'empty trace mean is NaN');
check(~o6.stat.ok && isnan(o6.stat.p) && contains(o6.stat.reason,'below minN'), ...
    sprintf('paired test refuses: "%s"',o6.stat.reason));
check(numel(o6.fig)==3, 'all three panels still rendered');
close all

% the old idiom on the same data, for contrast
validKeys6 = strcat(string(T6.animal),'_',string(T6.roiID));
sig6 = cellfun(@(x) logical(x(1)),T6.sigPk);
uk6 = unique(validKeys6);
v6 = arrayfun(@(i) all(sig6(validKeys6==uk6(i))), (1:numel(uk6))');
Tpop6 = T6(ismember(validKeys6,uk6(v6)),:);
oldCrashed = false;
try
    g = findgroups(Tpop6.dBdelta);
    splitapply(@(x) {mean(vertcat(x{:}),1,'omitnan')},Tpop6.dFF_PT_avg,g);
catch
    oldCrashed = true;
end
check(oldCrashed, 'the old splitapply idiom does crash on this same input');

%% ============ 3. single animal WITH cells ============
fprintf('\n=== 3. single animal with cells (TO0007) ===\n');
T7 = T(string(T.animal)=="TO0007",:);
o7 = plotCGCgroup(T7,'verbose',false);
check(o7.N.singleAnimal && sum(o7.valid)==15, ...
    sprintf('TO0007 contributes all 15 significant cells: %s',o7.N.label));
check(o7.stat.ok, 'paired test runs on a single-animal group');
close all

%% ============ 4. single cell ============
fprintf('\n=== 4. single cell ===\n');
firstValid = uk(find(o.valid,1));
T1 = T(key==firstValid,:);
o1 = plotCGCgroup(T1,'verbose',false);
check(all([o1.traces.n]==1), 'one cell per contrast');
check(~any([o1.traces.showBand]), 'no SEM band at n=1');
check(all(cellfun(@(s) all(isnan(s)), {o1.traces.sem})), 'trace SEM is NaN at n=1');
check(~o1.stat.ok && isnan(o1.stat.p), 'paired test refuses at n=1 (no NaN p-value on the figure)');
close all

%% ============ 5. options ============
fprintf('\n=== 5. options ===\n');
oAll = plotCGCgroup(T,'sigOnly',false,'verbose',false);
check(all([oAll.traces.n]==51), 'sigOnly=false widens the trace panel to all 51 cells');
check(oAll.stat.n==15, 'scatter/bar still use the significant set regardless of sigOnly');
close all
oSub = plotCGCgroup(T,'plots',{'scatter'},'verbose',false);
check(numel(oSub.fig)==1 && ~isfield(oSub,'traces'), 'plots subset honoured');
close all
oMin = plotCGCgroup(T,'minN',99,'verbose',false);
check(~oMin.stat.ok && contains(oMin.stat.reason,'minN'), 'minN gate refuses a large threshold');
close all

% mixed PT onset warns
Tm = T; Tm.PTsOnset(1) = 0.5;
lastwarn('');
plotCGCgroup(Tm,'plots',{'traces'},'verbose',false);
[~,wid] = lastwarn;
check(strcmp(wid,'plotCGCgroup:mixedPTonset'), 'mixed PT onset warns and omits the marker');
close all


%% ============ 5b. showCells toggle ============
fprintf('\n=== 5b. showCells ===\n');
oNo  = plotCGCgroup(T,'plots',{'traces'},'verbose',false);
oYes = plotCGCgroup(T,'plots',{'traces'},'showCells',true,'verbose',false);
check(isequal([oNo.traces.mean],[oYes.traces.mean]) && ...
      isequaln([oNo.traces.sem],[oYes.traces.sem]), ...
    'showCells changes only the drawing, not mean/SEM');
check(all(arrayfun(@(t) size(t.M,1)==t.n, oYes.traces)), ...
    'traces.M has one row per contributing cell');
check(all(arrayfun(@(t) size(t.M,2)==numel(t.t), oYes.traces)), ...
    'traces.M width matches the time vector');
check(all(arrayfun(@(t) height(t.cells)==t.n, oYes.traces)), ...
    'traces.cells is aligned with traces.M');
check(max(abs(mean(oYes.traces(1).M,1,'omitnan') - oYes.traces(1).mean)) < 1e-12, ...
    'mean of traces.M reproduces traces.mean');
% the outlier the figure is meant to reveal is reachable from the output
[mx,iMx] = max(max(oYes.traces(1).M,[],2));
check(isfinite(mx) && iMx>=1, ...
    sprintf('largest low-contrast cell peak %.2f is attributable to %s roi %s', ...
    mx, oYes.traces(1).cells.animal(iMx), oYes.traces(1).cells.roiID(iMx)));
close all

% empty and single-cell groups must still render with showCells on
o6c = plotCGCgroup(T6,'plots',{'traces'},'showCells',true,'verbose',false);
check(all([o6c.traces.n]==0), 'showCells with zero cells still reports n=0');
close all


%% ============ 5c. plotCellTrials (outlier follow-up) ============
fprintf('\n=== 5c. plotCellTrials ===\n');
oc = plotCellTrials(fullfile(d,'CGC_GroupD.mat'),"TO0007","7",'verbose',false);
check(numel(oc.trials)==2 && oc.trials(1).n==7 && oc.trials(2).n==11, ...
    sprintf('reps per contrast: %d low, %d high',oc.trials(1).n,oc.trials(2).n));
% the per-rep traces must average to the stored cell-average trace
rowLo = find(string(T.animal)=="TO0007" & string(T.roiID)=="7" & T.dBdelta==10);
check(max(abs(oc.trials(1).mean - T.dFF_PT_avg{rowLo})) < 1e-12, ...
    'mean of the plotted reps equals the stored dFF_PT_avg');
check(isequal(oc.trials(1).perRepMax, max(oc.trials(1).M,[],2)'), ...
    'perRepMax is the per-repetition maximum');
% the question the tool exists to answer: one bad rep, or many large ones?
nBig = sum(oc.trials(1).perRepMax > 1);
check(nBig >= 2, ...
    sprintf('low contrast: %d of %d reps exceed 1.0 dF/F -- not a single trace', ...
    nBig, oc.trials(1).n));
check(~isempty(oc.rawF) && size(oc.rawF(1).M,1)==oc.trials(1).n, ...
    'raw F panel has the same repetition count');
close all

% family-agnostic: same call shape works for BPN
ob = plotCellTrials(fullfile(d,'BPN_GroupD.mat'),"TO0007","7",'verbose',false);
check(numel(ob.trials)==6 && strcmp(ob.family,'BPN'), ...
    sprintf('BPN: %d levels resolved from the spec',numel(ob.trials)));
close all

threw = false;
try, plotCellTrials(fullfile(d,'CGC_GroupD.mat'),"TO9999","1",'verbose',false); catch, threw = true; end
check(threw,'unknown cell raises rather than plotting an empty panel');
close all

%% ============ 6. every group runs ============
fprintf('\n=== 6. all four groups ===\n');
for g = {'A','B','C','D'}
    og = plotCGCgroup(fullfile(d,sprintf('CGC_Group%s.mat',g{1})),'verbose',false);
    check(numel(og.fig)==3, sprintf('CGC_Group%s: %d/%d cells sig in both, test %s', ...
        g{1}, sum(og.valid), og.N.nCells, ternary(og.stat.ok,'ran','refused')));
    close all
end

%% ============ summary ============
fprintf('\n===============================\n');
fprintf('  %d passed, %d failed\n', nPass, nFail);
fprintf('===============================\n');
if nFail > 0; error('testPlotCGCgroup: %d check(s) failed', nFail); end

    function check(cond,label)
        if cond; nPass = nPass+1; fprintf('  PASS  %s\n',label);
        else;    nFail = nFail+1; fprintf('  FAIL  %s\n',label); end
    end
end

function v = ternary(c,a,b)
if c; v = a; else; v = b; end
end
