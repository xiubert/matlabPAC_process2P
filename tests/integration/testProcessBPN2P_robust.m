% testProcessBPN2P_robust
%
% Robustness checks of the BPN pipeline under parameter changes the user
% might plausibly try:
%   - baselineSec sweep: 0.5, 1.0 (current default)
%   - baselineSec > min(BPNsOnset): expect a clean error from combineDiffOnset
%   - synthetic single-onset table (only BPNsOnset=2): no merge, dFF and
%     raw F still onset-aligned consistently
%   - sub-frame pulse length: ceil(BPNmsPulseLen/1000*fr) must produce an
%     integer >= 1 frame window
%
% These checks operate on the TO0003 raw table.

addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();

dp = requireFixture(cfg.exampleAnimalDir,'TO0003 animal folder');
animal = 'TO0003';
rawFile = fullfile(dp, [animal '_anmlROI_BPNstimTable_raw.mat']);
S = load(rawFile, 'anmlROIbyStim');
Traw = S.anmlROIbyStim;

%% Helper: replicate processBPN2P's dFF prep on a given table for a given baselineSec
function T = buildDFF(T, baselineSec)
    T.t_total = rowfun(@(F,fr) {(1:size(F,2))/fr}, ...
        T,'InputVariables',{'SCALEDfissaFroi','frameRate'}, ...
        'ExtractCellContents',true,'OutputFormat','uniform');
    T.t_dFF = rowfun(@(t,onset) ...
        {t(find((t>onset-baselineSec),1,'first'):end)}, ...
        T,'InputVariables',{'t_total','BPNsOnset'}, ...
        'ExtractCellContents',true,'OutputFormat','uniform');
    T.dFF = rowfun(@(F,t,onset) ...
        {dFoFcalc(F,[find((t>onset-baselineSec),1,'first') ...
                     find((t<=onset),1,'last')],1)}, ...
        T,'InputVariables',{'SCALEDfissaFroi','t_total','BPNsOnset'}, ...
        'ExtractCellContents',true,'OutputFormat','uniform');
end

%% Helper: structural alignment check on the merged table vs input
function assertStructuralAlignment(Tm, T, baselineSec, label)
    for r = 1:height(Tm)
        inMask = (T.roiID == Tm.roiID(r)) & (T.BPNdBAmpl == Tm.BPNdBAmpl(r));
        inIdx = find(inMask);
        expectedRows = cell(numel(inIdx),1);
        for ii = 1:numel(inIdx)
            r0 = inIdx(ii);
            F0 = T.SCALEDfissaFroi{r0};
            fr0 = T.frameRate(r0);
            t0 = (1:size(F0,2))/fr0;
            sf = find(t0 > T.BPNsOnset(r0) - baselineSec, 1, 'first');
            if isempty(sf), sf = 1; end
            expectedRows{ii} = F0(:, sf:end);
        end
        minLen = min(cellfun(@(x) size(x,2), expectedRows));
        expectedRows = cellfun(@(x) x(:, 1:minLen), expectedRows, 'UniformOutput', false);
        expected = vertcat(expectedRows{:});
        assert(isequal(Tm.SCALEDfissaFroi{r}, expected), ...
            '[%s] SCALEDfissaFroi post-merge row %d does not match expected onset-aligned vertcat', label, r);
    end
end

%% Case 1: baselineSec = 1.0 (default)
fprintf('=== Case 1: baselineSec=1.0 ===\n');
T = buildDFF(Traw, 1.0);
Tm = combineDiffOnset(T, 1.0);
assertStructuralAlignment(Tm, T, 1.0, 'baselineSec=1.0');
fprintf('  height=%d, dFF dims=%s, OK\n', height(Tm), mat2str(size(Tm.dFF{1})));

%% Case 2: baselineSec = 0.5 (shorter baseline)
fprintf('=== Case 2: baselineSec=0.5 ===\n');
T = buildDFF(Traw, 0.5);
Tm = combineDiffOnset(T, 0.5);
assertStructuralAlignment(Tm, T, 0.5, 'baselineSec=0.5');
fprintf('  height=%d, dFF dims=%s, OK\n', height(Tm), mat2str(size(Tm.dFF{1})));
% frameStart should still be derivable
for r = 1:height(Tm)
    fs = find(Tm.t_dFF{r} >= Tm.BPNsOnset(r), 1, 'first');
    assert(~isempty(fs), 'No frameStart for row %d at baselineSec=0.5', r);
end

%% Case 3: baselineSec > min(BPNsOnset) should error cleanly
fprintf('=== Case 3: baselineSec=1.5 > min(BPNsOnset)=1 ===\n');
T = buildDFF(Traw, 1.5);
errored = false;
try
    Tm = combineDiffOnset(T, 1.5);
catch ME
    errored = strcmp(ME.identifier, 'combineDiffOnset:BaselineExceedsMinOnset');
end
assert(errored, 'Expected combineDiffOnset:BaselineExceedsMinOnset, got otherwise');
fprintf('  errored as expected\n');

%% Case 4: synthetic single-onset table (drop BPNsOnset=1 rows)
fprintf('=== Case 4: single-onset (BPNsOnset=2 only) ===\n');
Tonly2 = Traw(Traw.BPNsOnset == 2, :);
% baselineSec=1 is fine here since min(BPNsOnset)=2
Tonly2 = buildDFF(Tonly2, 1.0);
Tm = combineDiffOnset(Tonly2, 1.0);
assertStructuralAlignment(Tm, Tonly2, 1.0, 'single-onset');
% No merging, so output rows == input rows
assert(height(Tm) == height(Tonly2), ...
    'Expected height %d for single-onset (no merge), got %d', height(Tonly2), height(Tm));
% frameStart derivation should still match for every row
for r = 1:height(Tm)
    fs = find(Tm.t_dFF{r} >= Tm.BPNsOnset(r), 1, 'first');
    expected = round(Tm.frameRate(r) * 1.0);  % baselineSec=1
    assert(abs(fs - expected) <= 1, ...
        'Single-onset row %d frameStart %d != expected %d', r, fs, expected);
end
fprintf('  height=%d (no merge), dFF dims=%s, OK\n', height(Tm), mat2str(size(Tm.dFF{1})));

%% Case 5: sub-frame pulse length
fprintf('=== Case 5: sub-frame pulse length ===\n');
fr = 5;
for BPNmsPulseLen = [50 100 150 200 400]
    nFramesPostPulse = 2;
    nFrameWindow = ceil(BPNmsPulseLen/1000*fr) + nFramesPostPulse;
    assert(nFrameWindow >= 3, 'nFrameWindow must be at least nFramesPostPulse+1 for any non-zero pulse');
    fprintf('  pulse=%d ms -> nFrameWindow=%d (>= 3)\n', BPNmsPulseLen, nFrameWindow);
end

fprintf('\nAll robustness assertions passed.\n');
