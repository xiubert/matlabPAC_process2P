% testFRAgroup
%
% Covers the FRA-by-group pipeline:
%   1. stimGroupSpec registers FRA and the spec is well formed
%   2. FRAmap2table round-trips pkDFF/sig from the struct without change,
%      and refuses a pre-fix (no .params) struct
%   3. tableFRAmetrics: BF agrees with anmlFRA2BF, bandwidth arithmetic is
%      right on a hand-built tuning curve, threshold respects minBand
%   4. the sham noise-floor control is computed and sane
%   5. real group files load, validate and plot
%
% Run: matlab -batch "addpath('tests'); runTests('Filter','testFRAgroup')"

addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();
set(0,'defaultFigureVisible','off');
nPass = 0;

%% ---------------------------------------------------------------- 1
spec = stimGroupSpec('FRA');
assert(strcmp(spec.family,'FRA'));
assert(strcmp(spec.varname,'anmlROIbyFRA'));
assert(strcmp(spec.levelVar,'dBSPL'));
assert(isequal(spec.stimVars,{'freqHz','dBSPL'}),'FRA stim is (freq,level), both');
assert(isempty(spec.suffixRaw),'FRA has no _raw stage');
assert(isempty(spec.cellTrialVar),'FRA keeps no per-trial trace column');
assert(ismember('FRA',stimGroupSpec()),'FRA missing from the family list');
fprintf('PASS 1: stimGroupSpec registers FRA\n'); nPass = nPass+1;

%% ---------------------------------------------------------------- 2
dp = fullfile(cfg.animalsRoot,'TO0003');
haveReal = isfolder(dp);
if haveReal
    R = load(fullfile(dp,'TO0003_FRAmap.mat'),'FRAmap'); R = R.FRAmap;
    T = FRAmap2table(R,'TO0003');
    nDB = numel(R.dBlist); nFreq = numel(R.freqList);
    assert(height(T) == numel(R.BFuDB)*nDB*nFreq,'row count wrong');
    ok = true;
    for dbR = 1:nDB
        for fq = 1:nFreq
            c = R.dBFreqMap{dbR,fq};
            sub = T(T.dBSPL==R.dBlist(dbR) & T.freqHz==R.freqList(fq),:);
            [~,ord] = sort(str2double(sub.roiID)); sub = sub(ord,:);
            ok = ok && isequaln(sub.pkDFF,c.pkDFF) && isequal(sub.sig,c.sigPkDFF);
        end
    end
    assert(ok,'table pkDFF/sig diverge from the FRAmap struct');

    stale = rmfield(R,'params');
    threw = false;
    try, FRAmap2table(stale,'TO0003'); catch e
        threw = strcmp(e.identifier,'FRAmap2table:staleInput'); end
    assert(threw,'a pre-fix FRAmap struct must be refused, not silently converted');
    fprintf('PASS 2: table round-trips exactly; pre-fix struct refused\n'); nPass = nPass+1;
else
    fprintf('SKIP 2: %s not mounted\n',dp);
end

%% ---------------------------------------------------------------- 3
% hand-built cell: significant band at 3 adjacent frequencies, 50 dB only
freqs = 1000*2.^(0:0.125:0.875);      % 8 freqs, 0.125 oct apart
levels = [30 50 70];
[F,D] = meshgrid(freqs,levels);
pk  = 0.01*ones(numel(F),1);
sig = false(numel(F),1);
tbl = table(repmat("X",numel(F),1),repmat("1",numel(F),1), ...
    F(:),D(:),pk,sig,'VariableNames',{'animal','roiID','freqHz','dBSPL','pkDFF','sig'});
band = ismember(tbl.freqHz,freqs(3:5)) & tbl.dBSPL==50;
tbl.sig(band) = true;  tbl.pkDFF(band) = [0.2;0.5;0.3];
% one isolated significant frequency at 30 dB -- must NOT set threshold
iso = tbl.freqHz==freqs(8) & tbl.dBSPL==30;
tbl.sig(iso) = true; tbl.pkDFF(iso) = 0.4;

% minBand pinned explicitly: these assertions depend on it, so a change to
% the DEFAULT (now 2) must not silently pass here
M = tableFRAmetrics(tbl,'nConsec',1,'minBand',3);
assert(tableFRAmetrics(tbl,'nConsec',1).minBand == 2,'minBand default should now be 2');
ci = M.cellInfo;
assert(ci.threshold == 50, ...
    'minBand=3 should ignore the isolated 30 dB frequency; got threshold %g',ci.threshold);
assert(ci.bf == freqs(4),'BF should be the strongest significant frequency');
expBW = log2(freqs(5)/freqs(3));
assert(abs(M.bwByLevel(2)-expBW) < 1e-12, ...
    'bandwidth should be %.4f oct across the 3-frequency band, got %.4f',expBW,M.bwByLevel(2));
assert(isnan(ci.bw20),'threshold 50 + 1 level = 70 dB, which has no response -> NaN');
% minBand=1 instead lets the isolated frequency win
M1 = tableFRAmetrics(tbl,'nConsec',1,'minBand',1);
assert(M1.cellInfo.threshold == 30,'minBand=1 should accept the isolated frequency');
fprintf('PASS 3: BW = %.3f oct across a 3-frequency band; minBand gates threshold\n',expBW);
nPass = nPass+1;

if haveReal
    Mr = tableFRAmetrics(T,'nConsec',1,'minBand',3);
    cir = Mr.cellInfo;
    % align by roiID: findgroups sorts roiID as STRINGS ("1","10","11",...)
    [~,ord] = sort(str2double(cir.roiID));
    assert(isequal(cir.bf(ord),R.BFuDB), ...
        'bf must agree with anmlFRA2BF once rows are aligned by numeric roiID');
    fprintf('PASS 3b: bf matches anmlFRA2BF for all %d cells\n',height(cir)); nPass = nPass+1;
end

%% ---------------------------------------------------------------- 4
if haveReal
    s = Mr.shamControl;
    assert(s.ok,'sham control should be computable for a real animal');
    assert(s.realRate>=0 && s.realRate<=1 && s.shamRate>0 && s.shamRate<=1);
    assert(isfinite(s.ratio));
    fprintf('PASS 4: sham control real %.3f / sham %.3f = %.2fx\n', ...
        s.realRate,s.shamRate,s.ratio); nPass = nPass+1;
end

%% ---------------------------------------------------------------- 5
agg = cfg.aggregateDir;
gf  = fullfile(agg,'FRA_GroupA.mat');
if isfile(gf)
    [Tg,info,rep] = loadStimGroup(gf,'family','FRA','strict',false);
    assert(istable(Tg) && height(Tg)>0);
    assert(~isempty(info) && strcmp(info.family,'FRA'),'group file should be stamped');
    assert(rep.ok,'group A failed validation: %s', ...
        strjoin(cellstr(rep.problems.message(rep.problems.severity=="error")),'; '));
    o = plotFRAgroup(gf,'verbose',false);
    assert(numel(o.fig)==4,'expected 4 figures, got %d',numel(o.fig));
    assert(o.N.nCells==info.nCells);
    assert(sum(o.threshold.counts)==o.N.nCells,'threshold counts must cover every cell');
    assert(all(o.bw.values>=0),'bandwidth cannot be negative');
    close(o.fig);
    fprintf('PASS 5: FRA_GroupA loads, validates and plots (%s)\n',o.N.label); nPass = nPass+1;
else
    fprintf('SKIP 5: %s not present\n',gf);
end

fprintf('\nAll %d test groups passed.\n',nPass);
