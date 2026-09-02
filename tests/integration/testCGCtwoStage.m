function testCGCtwoStage()
% testCGCtwoStage  Phase 5: CGC _raw two-stage split.
%
%   matlab -batch "addpath('tests'); runTests('Filter','testCGCtwoStage')"
%
% Asserts that processCGC reads the _raw artifact, writes the processed one,
% never mutates its input, and still handles a pre-split animal that has only
% the un-suffixed file (stripping its derived columns first).
%
% A per-animal CGC dataset is synthesized from CGC_GroupD by taking one
% animal and removing the derived columns, since the real per-animal tables
% are not on this machine.

addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();
d = requireFixture(cfg.aggregateDir,'TOMT aggregate group files');
set(0,'DefaultFigureVisible','off'); c0 = onCleanup(@() set(0,'DefaultFigureVisible','on'));

nPass = 0; nFail = 0;
[tmp,tmpCleanup] = testSandbox('CGCtwoStage'); %#ok<ASGLU>
animalDir = fullfile(tmp,'TO0007'); mkdir(animalDir);

%% ---- synthesize a per-animal _raw table ----
G = getfield(load(fullfile(d,'CGC_GroupD.mat'),'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>
Traw = G(string(G.animal)=="TO0007",:);
spec = stimGroupSpec('CGC');
derived = intersect(spec.derivedVars, Traw.Properties.VariableNames);
Traw = removevars(Traw, derived);
stimTable = unique(Traw(:,{'dBdelta'}));                 %#ok<NASGU>
tifStimParamTable = table((1:2)','VariableNames',{'placeholder'}); %#ok<NASGU>

rawPath  = fullfile(animalDir,'TO0007_anmlROI_CGCstimTable_raw.mat');
procPath = fullfile(animalDir,'TO0007_anmlROI_CGCstimTable.mat');
S = struct('anmlROIbyStim',Traw,'stimTable',stimTable, ...
    'tifStimParamTable',tifStimParamTable,'dataPath',animalDir);
save(rawPath,'-struct','S','-v7.3');

check(~ismember('dFF_PT_avg',Traw.Properties.VariableNames), ...
    'synthesized _raw has no derived columns');

%% ============ 1. raw -> processed ============
fprintf('\n=== 1. two-stage write ===\n');
rawHashBefore = fileHash(rawPath);
runProcessCGC(animalDir,'TO0007');

check(isfile(procPath), 'processCGC wrote the processed file');
check(isfile(rawPath) && strcmp(fileHash(rawPath),rawHashBefore), ...
    '_raw input is byte-identical after the run (never mutated)');

W = whos('-file',procPath);
names = sort({W.name});
check(isequal(names,{'anmlROIbyStim','dataPath','stimTable','tifStimParamTable'}), ...
    sprintf('processed bundle holds all four variables: %s', strjoin(names,', ')));

Tp = getfield(load(procPath,'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>
rep = validateStimGroup(Tp,'CGC','verbose',false);
check(rep.ok, 'processed per-animal table validates');
check(all(ismember({'dFF_PT','dFF_PT_avg','pkPT','sigPk'}, ...
    Tp.Properties.VariableNames)), 'derived columns present after processing');

% the processed output must match what the repaired group file holds
Gd = G(string(G.animal)=="TO0007",:);
dPT = max(cellfun(@(a,b) max(abs(a-b),[],'all'), Tp.dFF_PT_avg, Gd.dFF_PT_avg));
dPk = max(abs(cellfun(@(c) c(1),Tp.pkPT) - cellfun(@(c) c(1),Gd.pkPT)));
check(dPT < 1e-12 && dPk < 1e-12, ...
    sprintf('recomputed dFF_PT_avg and pkPT match the repaired group data (dev %.3g / %.3g)',dPT,dPk));

%% ============ 2. re-run is repeatable ============
fprintf('\n=== 2. re-run ===\n');
h1 = fileHash(procPath);
runProcessCGC(animalDir,'TO0007');
Tp2 = getfield(load(procPath,'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>
check(strcmp(fileHash(rawPath),rawHashBefore), '_raw still untouched after a second run');
check(isequaln(Tp(:,{'dFF_PT_avg','pkPT','sigPk'}),Tp2(:,{'dFF_PT_avg','pkPT','sigPk'})), ...
    'second run reproduces the same derived values');
if strcmp(h1,fileHash(procPath))
    fprintf('  (processed file is byte-identical too)\n');
end

%% ============ 3. legacy pre-split layout ============
fprintf('\n=== 3. legacy fallback ===\n');
legacyDir = fullfile(tmp,'TO0006'); mkdir(legacyDir);
Tleg = G(string(G.animal)=="TO0006",:);        % raw AND derived columns
Tleg.dFF_PT_preDRCf0 = Tleg.dFF_PT_avg;        % a stale column from another script version
Tleg.t_dFF_PT = Tleg.t_dFF_DRC;
Sl = struct('anmlROIbyStim',Tleg,'stimTable',unique(Tleg(:,{'dBdelta'})), ...
    'tifStimParamTable',tifStimParamTable,'dataPath',legacyDir);
legacyProc = fullfile(legacyDir,'TO0006_anmlROI_CGCstimTable.mat');
save(legacyProc,'-struct','Sl','-v7.3');
check(~isfile(fullfile(legacyDir,'TO0006_anmlROI_CGCstimTable_raw.mat')), ...
    'legacy animal has no _raw file');

[wid,~] = runProcessCGC(legacyDir,'TO0006');
check(strcmp(wid,'processCGC:legacyLayout'), 'legacy layout warns');

Tl = getfield(load(legacyProc,'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>
check(~any(ismember({'dFF_PT_preDRCf0','t_dFF_PT'},Tl.Properties.VariableNames)), ...
    'stale foreign columns are stripped, not carried through');
repL = validateStimGroup(Tl,'CGC','verbose',false);
check(repL.ok, 'legacy-migrated table validates');

%% ============ 4. no input at all ============
fprintf('\n=== 4. missing input ===\n');
emptyDir = fullfile(tmp,'TO0099'); mkdir(emptyDir);
threw = false;
try
    runProcessCGC(emptyDir,'TO0099');
catch ME
    threw = strcmp(ME.identifier,'processCGC:noInput');
end
check(threw,'missing input raises processCGC:noInput');

%% ============ 5. stimParam2ROI writes _raw ============
fprintf('\n=== 5. stimParam2ROI target ===\n');
src = fileread(fullfile(cfg.repoRoot,'helperFcns','dataOrg','stimParam2ROI.m'));
check(contains(src,'_anmlROI_CGCstimTable_raw.mat'), ...
    'stimParam2ROI saves the CGC table to _raw');
check(~contains(src,"'_anmlROI_CGCstimTable.mat'"), ...
    'stimParam2ROI no longer writes the un-suffixed CGC file');
check(strcmp(spec.suffixRaw,'_anmlROI_CGCstimTable_raw.mat') && ...
      strcmp(spec.suffixProcessed,'_anmlROI_CGCstimTable.mat'), ...
    'stimGroupSpec agrees with both filenames');


%% ============ summary ============
fprintf('\n===============================\n');
fprintf('  %d passed, %d failed\n', nPass, nFail);
fprintf('===============================\n');
if nFail > 0; error('testCGCtwoStage: %d check(s) failed', nFail); end

    function check(cond,label)
        if cond; nPass = nPass+1; fprintf('  PASS  %s\n',label);
        else;    nFail = nFail+1; fprintf('  FAIL  %s\n',label); end
    end
end

%% ---- helpers ----
function [wid,wmsg] = runProcessCGC(dataPath,animal) %#ok<INUSD>
% Run the processCGC script in a scope holding only dataPath and animal.
lastwarn('');
evalin('base','clearvars');
assignin('base','dataPath',dataPath);
assignin('base','animal',animal);
evalin('base','processCGC;');
evalin('base','close all');
[wmsg,wid] = lastwarn;
end

function h = fileHash(f)
d = dir(f);
fid = fopen(f,'r'); bytes = fread(fid,Inf,'*uint8'); fclose(fid);
h = sprintf('%d_%s', d.bytes, mat2str(typecast(sum(double(bytes)),'uint64')));
end
