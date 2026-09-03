% testProcessBPN2P_e2e
%
% End-to-end test of the new processBPN2P pipeline against TO0003:
%   1. Load TO0003_anmlROI_BPNstimTable_raw.mat (produced by stimParam2ROI).
%   2. Replicate processBPN2P's compute path (without saving) up through
%      pkFcalc.
%   3. Assert key invariants that the revision plan was designed to
%      enforce:
%        - combineDiffOnset returns one row per (ROI, BPNdBAmpl) group
%        - SCALEDfissaFroi and dFF in the merged table are matrices in
%          cells, not nested cells
%        - SCALEDfissaFroi is onset-aligned post-merge: per-trial peaks
%          for sub-rows from BPNsOnset=1 and BPNsOnset=2 land at the
%          same frame neighborhood
%        - frameStart computed from t_dFF / BPNsOnset matches the dFF
%          stim-onset index for every row, independent of baselineSec
%        - resultsTable columns map to sigPeak / sig / pkResp cleanly
%
% Run after stimParam2ROI has generated the _raw.mat.

addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();

dp = requireFixture(cfg.exampleAnimalDir,'TO0003 animal folder');
animal = 'TO0003';
rawFile = fullfile(dp, [animal '_anmlROI_BPNstimTable_raw.mat']);
assert(isfile(rawFile), 'Missing %s; run stimParam2ROI(dataPath) first.', rawFile);

S = load(rawFile, 'anmlROIbyStim');
T = S.anmlROIbyStim;
fprintf('Loaded raw table: %d rows\n', height(T));
assert(ismember('BPNsOnset', T.Properties.VariableNames), 'Expected BPNsOnset column.');

%% Compute t_total / t_dFF / dFF per row (mirror processBPN2P)
baselineSec = 1;

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

%% Pre-merge sanity
nRowsPre = height(T);
nRoi = numel(unique(T.roiID));
nDb  = numel(unique(T.BPNdBAmpl));
nOnset = numel(unique(T.BPNsOnset));
fprintf('Pre-merge: %d rows = %d ROIs x %d dB x %d onsets (expected %d)\n', ...
    nRowsPre, nRoi, nDb, nOnset, nRoi*nDb*nOnset);
assert(nRowsPre == nRoi*nDb*nOnset);

%% Merge
Tm = combineDiffOnset(T, baselineSec);
fprintf('Post-merge: %d rows (expected %d = %d ROIs x %d dB)\n', ...
    height(Tm), nRoi*nDb, nRoi, nDb);
assert(height(Tm) == nRoi*nDb);

%% Cell-content shape audit
shapeChecks = {'SCALEDfissaFroi','rawFroi','moCorRawFroi','fissaFroi','dFF','t_total','t_dFF'};
for k = 1:numel(shapeChecks)
    nm = shapeChecks{k};
    if ~ismember(nm, Tm.Properties.VariableNames), continue; end
    sample = Tm.(nm){1};
    assert(isnumeric(sample), '%s should be numeric, got %s', nm, class(sample));
    assert(~iscell(sample), '%s contents must not be a nested cell', nm);
    fprintf('  %-18s -> %s, dims = %s\n', nm, class(sample), mat2str(size(sample)));
end

%% SCALEDfissaFroi onset alignment: structural check by hand-reconstruction
% For every merged row, identify the contributing input rows, apply the
% per-row trim + common minLen trim, vertcat in input order, and assert
% equality with the post-merge cell content. This verifies the
% combineDiffOnset trim algorithm exactly, independent of evoked-response
% presence in the data.
fr = unique(T.frameRate);
assert(isscalar(fr), 'Test assumes uniform frameRate across rows.');
expectStimFrame = round(fr * baselineSec);

for checkRow = 1:height(Tm)
    inMask = (T.roiID == Tm.roiID(checkRow)) & (T.BPNdBAmpl == Tm.BPNdBAmpl(checkRow));
    inIdx = find(inMask);

    expectedRows = cell(numel(inIdx), 1);
    for ii = 1:numel(inIdx)
        r0 = inIdx(ii);
        F0 = T.SCALEDfissaFroi{r0};
        fr0 = T.frameRate(r0);
        t0 = (1:size(F0,2)) / fr0;
        sf = find(t0 > T.BPNsOnset(r0) - baselineSec, 1, 'first');
        if isempty(sf), sf = 1; end
        expectedRows{ii} = F0(:, sf:end);
    end
    minLen = min(cellfun(@(x) size(x,2), expectedRows));
    expectedRows = cellfun(@(x) x(:, 1:minLen), expectedRows, 'UniformOutput', false);
    expected = vertcat(expectedRows{:});

    assert(isequal(Tm.SCALEDfissaFroi{checkRow}, expected), ...
        'SCALEDfissaFroi post-merge for row %d does not match expected onset-aligned vertcat', checkRow);
end
fprintf('Structural alignment check passed for all %d merged rows.\n', height(Tm));

%% dFF stim onset via t_dFF / BPNsOnset matches expected frame
for r = 1:height(Tm)
    t_dFF = Tm.t_dFF{r};
    fs = find(t_dFF >= Tm.BPNsOnset(r), 1, 'first');
    if isempty(fs)
        error('No frame in t_dFF{%d} reaches BPNsOnset=%g', r, Tm.BPNsOnset(r));
    end
    if r == 1
        fprintf('frameStart from t_dFF: %d (expected ~%d = frameRate*baselineSec)\n', ...
            fs, expectStimFrame);
    end
    assert(abs(fs - expectStimFrame) <= 1, ...
        'frameStart %d for row %d differs from expected %d', fs, r, expectStimFrame);
end

%% Trial-average + pkFcalc
Tm.dFF_avg = rowfun(@(F) {mean(F,1,'omitnan')}, ...
    Tm,'InputVariables',{'dFF'}, ...
    'ExtractCellContents',true,'OutputFormat','uniform');

pkBPNsigSD = 2;
nFramesPostPulse = 2;
resultsTable = rowfun(@(F, t_dFF, frameRate, BPNsOnset, BPNmsPulseLen) ...
    pkFcalc(F, ...
        find(t_dFF >= BPNsOnset, 1, 'first'), ...
        ceil(BPNmsPulseLen/1000*frameRate) + nFramesPostPulse, ...
        pkBPNsigSD), ...
    Tm, ...
    'InputVariables',{'dFF_avg','t_dFF','frameRate','BPNsOnset','BPNmsPulseLen'}, ...
    'ExtractCellContents',true,'NumOutputs',5,'OutputFormat','cell');

Tm.sigPeak = resultsTable(:,1);
Tm.sig     = resultsTable(:,2);
Tm.pkResp  = resultsTable(:,3);
assert(all(cellfun(@(x) isnumeric(x) || islogical(x), Tm.sigPeak)), 'sigPeak should be numeric/logical');
assert(all(cellfun(@(x) isnumeric(x) || islogical(x), Tm.sig)),     'sig should be numeric/logical');
assert(all(cellfun(@(x) isnumeric(x), Tm.pkResp)),                  'pkResp should be numeric');

fprintf('\nAll assertions passed.\n');
