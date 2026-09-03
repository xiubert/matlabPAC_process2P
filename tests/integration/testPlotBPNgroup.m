function testPlotBPNgroup()
% testPlotBPNgroup  Phase 3: BPN group plotter + the RLF small-n unification.
%
%   matlab -batch "addpath('tests'); runTests('Filter','testPlotBPNgroup')"
%
% Figures are built with 'Visible','off'; the assertions are on the returned
% data struct, so the panels are checked numerically rather than by eye.

addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();
d = requireFixture(cfg.aggregateDir,'TOMT aggregate group files');
set(0,'DefaultFigureVisible','off'); c0 = onCleanup(@() set(0,'DefaultFigureVisible','on'));

nPass = 0; nFail = 0;

T = getfield(load(fullfile(d,'BPN_GroupD.mat'),'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>

%% ============ 1. full group ============
fprintf('\n=== 1. full group ===\n');
o = plotBPNgroup(fullfile(d,'BPN_GroupD.mat'),'nConsec',2,'verbose',false);
check(o.N.nCells==51 && o.N.nAnimals==2, sprintf('counts: %s',o.N.label));
check(isequal(o.levels,[30 40 50 60 70 80]), 'levels default to all present, ascending');
check(numel(o.traces)==6 && all([o.traces.n]==51), 'traces: one entry per level, 51 cells each');
check(all(cellfun(@(s) isequal(size(s),[1 15]), {o.traces.mean})), 'trace means are 1x15');
check(all([o.traces.showBand]), 'SEM bands shown at n=51');
check(numel(o.peak)==6 && all([o.peak.n]==51), 'peak: 51 cells per level');
check(o.rlf.nIncluded==25 && o.rlf.nTotal==51, ...
    sprintf('RLF matches the saved run: %d/%d included',o.rlf.nIncluded,o.rlf.nTotal));
close all

% cross-check against the researcher's saved RLF file
S = load(fullfile(d,'RLF_GroupD_nConsec2.mat'),'rlf');
check(isequal(o.rlf.nIncluded,S.rlf.nIncluded) && ...
      max(abs(o.rlf.meanRLF - S.rlf.meanRLF)) < 1e-12, ...
      'plotBPNgroup RLF reproduces RLF_GroupD_nConsec2.mat exactly');

% traces reproduce the processBPN2P Plot-3 idiom
k = 3;  % 50 dB
rows = find(T.BPNdBAmpl==o.levels(k));
M = vertcat(T.dFF_avg{rows});
check(max(abs(o.traces(k).mean - mean(M,1,'omitnan'))) < 1e-12, ...
    'trace mean matches the direct vertcat/mean idiom');
check(max(abs(o.traces(k).sem - SEMcalc(M,1))) < 1e-12, ...
    'trace SEM matches SEMcalc at n>1 (no behaviour change for real groups)');

%% ============ 2. new tableRLF fields ============
fprintf('\n=== 2. tableRLF additions ===\n');
r = tableRLF(T,'nConsec',2);
check(isfield(r,'nAnimals') && r.nAnimals==2, sprintf('nAnimals reported: %d',r.nAnimals));
check(isfield(r,'nAnimalsIncl') && r.nAnimalsIncl<=r.nAnimals, ...
    sprintf('nAnimalsIncl reported: %d',r.nAnimalsIncl));
check(isfield(r,'showBand') && r.showBand, 'showBand true at n>1');

%% ============ 3. single animal ============
fprintf('\n=== 3. single animal ===\n');
T6 = T(string(T.animal)=="TO0006",:);
o6 = plotBPNgroup(T6,'nConsec',2,'verbose',false);
check(o6.N.singleAnimal && o6.N.nCells==24, sprintf('single animal: %s',o6.N.label));
check(all([o6.traces.n]==24) && all([o6.traces.showBand]), 'traces still band at 24 cells');
check(o6.rlf.nAnimalsIncl==1, 'RLF reports 1 mouse');
close all

%% ============ 4. single cell ============
fprintf('\n=== 4. single cell ===\n');
key = strcat(string(T.animal),'_',string(T.roiID));
T1 = T(key==key(1),:);
o1 = plotBPNgroup(T1,'nConsec',2,'verbose',false);
check(all([o1.traces.n]==1), 'one cell per level');
check(~any([o1.traces.showBand]), 'no SEM band at n=1');
check(all(cellfun(@(s) all(isnan(s)), {o1.traces.sem})), 'trace SEM is NaN at n=1, not zeros');
r1 = tableRLF(T1,'nConsec',2);
check(r1.nIncluded<=1 && (r1.nIncluded==0 || all(isnan(r1.semRLF))), ...
    'tableRLF semRLF is NaN at n<2 (was zeros via SEMcalc)');
close all

%% ============ 5. zero cells pass inclusion ============
fprintf('\n=== 5. zero included cells ===\n');
o0 = plotBPNgroup(T,'nConsec',99,'traceCells','included','verbose',false);
check(o0.rlf.nIncluded==0, 'no cells pass nConsec=99');
check(o0.Nplot.empty && o0.Nplot.nCells==0, 'plotted set is empty');
check(all([o0.traces.n]==0), 'trace panels report n=0 without erroring');
check(all(cellfun(@(m) all(isnan(m)), {o0.traces.mean})), 'empty trace mean is NaN, not zeros');
check(all(isnan([o0.peak.mean])), 'empty peak mean is NaN');
close all

%% ============ 6. options ============
fprintf('\n=== 6. options ===\n');
oL = plotBPNgroup(T,'levels',[40 60],'plots',{'traces'},'verbose',false);
check(numel(oL.traces)==2 && isempty(oL.rlf), 'levels + plots subset honoured');
check(numel(oL.fig)==1, 'only the requested panel is created');
close all
oI = plotBPNgroup(T,'nConsec',2,'traceCells','included','verbose',false);
check(oI.Nplot.nCells==oI.rlf.nIncluded, ...
    sprintf('traceCells=included restricts to %d cells',oI.Nplot.nCells));
check(all([oI.traces.n]==oI.rlf.nIncluded), 'every level uses the included set');
close all

% mixed stimulus timing warns rather than mislabelling
Tm = T; Tm.BPNmsPulseLen(1) = 999;
lastwarn('');
plotBPNgroup(Tm,'plots',{'traces'},'verbose',false);
[~,wid] = lastwarn;
check(strcmp(wid,'plotBPNgroup:mixedStimTiming'), 'mixed pulse length warns and omits markers');
close all

%% ============ 7. plotRLF degenerate rendering ============
fprintf('\n=== 7. plotRLF ===\n');
f = figure; ax = axes(f);
[~,~,hErr] = plotRLF(tableRLF(T1,'nConsec',1),'ax',ax);
check(isempty(hErr), 'plotRLF draws no error bar when SEM is all NaN');
close(f);
f = figure; ax = axes(f);
[~,hLine,~] = plotRLF(tableRLF(T,'nConsec',99),'ax',ax);
check(isgraphics(hLine), 'plotRLF renders an empty RLF without erroring');
close(f);

% legacy struct without showBand/nAnimals still plots
legacy = load(fullfile(d,'RLF_GroupA_nConsec2.mat'),'rlf');
f = figure; ax = axes(f);
plotRLF(legacy.rlf,'ax',ax,'showCells',true);
check(true,'plotRLF accepts a legacy RLF struct (no showBand/nAnimalsIncl)');
close(f);

%% ============ summary ============
fprintf('\n===============================\n');
fprintf('  %d passed, %d failed\n', nPass, nFail);
fprintf('===============================\n');
if nFail > 0; error('testPlotBPNgroup: %d check(s) failed', nFail); end

    function check(cond,label)
        if cond; nPass = nPass+1; fprintf('  PASS  %s\n',label);
        else;    nFail = nFail+1; fprintf('  FAIL  %s\n',label); end
    end
end
