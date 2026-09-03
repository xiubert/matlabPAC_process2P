% testFRAsignificance
%
% Covers the FRA significance rework:
%   1. pkFcalc's new baseIDX argument is backward compatible (BPN/CGC safe)
%   2. sigPkDFF is a LOGICAL mask and the sig map holds the peak, not peak^2
%   3. significance is tested once per (cell,freq,dB) on the trial average,
%      so duplicating trials does not change the result
%   4. the onset column is the first frame at or after tone onset
%   5. missing freq/dB pairs are filled without breaking the reshape
%   6. the peak window is ceil'd, not floored
%   7. the linear index ordering anmlFRA2BF / anmlFRA2dPrime rely on holds
%   8. real TO0003 values, including the ROI 1 case that the inclusive
%      baseline would reject
%
% Run: matlab -batch "addpath('tests'); runTests('Filter','testFRAsignificance')"

addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();

[scratch,scratchCleanup] = testSandbox('FRAsignificance'); %#ok<ASGLU>
nPass = 0;

%% ---------------------------------------------------------------- 1
% pkFcalc: 4-arg form must equal 5-arg form with the legacy window
rng(1);
F = cumsum(randn(6,24),2);
Fc = {F; cumsum(randn(3,24),2)};        % column cell; row cells never worked
for frameStart = [5 9]
    o = cell(1,5); n = cell(1,5);
    [o{:}] = pkFcalc(F,frameStart,4,2);
    [n{:}] = pkFcalc(F,frameStart,4,2,1:frameStart);
    assert(isequaln(o,n),'pkFcalc default baseline changed (matrix)');
    [o{1:3}] = pkFcalc(Fc,frameStart,4,2);
    [n{1:3}] = pkFcalc(Fc,frameStart,4,2,1:frameStart);
    assert(isequaln(o(1:3),n(1:3)),'pkFcalc default baseline changed (cell)');
end
% and the exclusive window must actually do something. This fixture mirrors
% the real ROI 1 geometry: 3 flat baseline frames, then a response that is
% already high at the onset frame itself.
tr = [0 0.01 -0.01 1.0 0.6 0.2];        % onset at frame 4
[~,sIncl] = pkFcalc(tr,4,3,2);           % baseline 1:4 swallows the onset frame
[~,sExcl] = pkFcalc(tr,4,3,2,1:3);       % baseline strictly pre-onset
assert(sExcl && ~sIncl, ...
    ['exclusive baseline should keep this response significant and the ' ...
     'inclusive one should reject it (got excl=%d incl=%d)'],sExcl,sIncl);
fprintf('PASS 1: pkFcalc 4-arg == 5-arg legacy window; BPN/CGC unaffected\n'); nPass = nPass+1;

%% ---------------------------------------------------------------- fixture
% Synthetic map animal: 2 freqs x 2 dB, known response injected into one
% (cell, freq, dB). Pulse names use the real format parsed by
% extractMapPulseParams.
fs = 5; ISI = 3; stimDelay = 4; nCell = 4;
framesPerPulse = ISI*fs;                  % 15
freqs = [8409 25937]; dBs = [30 70];
onsets = [0.6 1.0];                       % two onset groups, as in TO0003
msPulse = 400;
nWinExpect = ceil(msPulse/1000*fs) + 2;   % nFramesPostPulse = 2

conds = {}; k = 0;
for rep = 1:3
    for f = freqs
        for d = dBs
            k = k+1;
            conds{k} = struct('f',f,'d',d,'o',onsets(mod(k,2)+1)); %#ok<SAGROW>
        end
    end
end
nPulse = numel(conds);

TARGET.cell = 2; TARGET.f = freqs(2); TARGET.d = dBs(2); TARGET.amp = 0.5;
baselineF = 100;

nFrames = stimDelay*fs + framesPerPulse*nPulse;
% low-level noise is required, not cosmetic: with a perfectly flat baseline
% the SD is 0, pkFcalc's zero2nan turns the threshold into NaN, and nothing
% can ever be significant
rng(7);
rawF = baselineF*(1 + 0.001*randn(nCell,nFrames));
pulse = struct('pulsename',{});
for p = 1:nPulse
    c = conds{p};
    pulse(p).pulsename = sprintf( ...
        '%dHz_%ddB_TestTone_%dmsPulse_at_%gs_10msRamp_1500msTotal_Fs250kHz', ...
        c.f,c.d,msPulse,c.o);
    if c.f==TARGET.f && c.d==TARGET.d
        base = stimDelay*fs + (p-1)*framesPerPulse;
        onsetFrameLocal = find((0:framesPerPulse-1)/fs >= c.o - 1e-9,1,'first');
        hit = base + onsetFrameLocal + 1;   % 1 frame after onset, like real Ca
        rawF(TARGET.cell,hit) = baselineF*(1+TARGET.amp);
    end
end
params = struct('stimDelay',stimDelay,'totalPulses',nPulse,'ISI',ISI, ...
    'PulseHiTime',0,'PulseLoTime',0);
save(fullfile(scratch,'AA9999AAAA_00001_00001_Pulses.mat'),'pulse','params');

tifFileList = struct('map',struct( ...
    'name','AA9999AAAA_00001_00001.tif','folder',scratch, ...
    'treatment','none FRAmap','nFrames',nFrames,'frameRate',fs, ...
    'SCALEDfissaFroi',rawF));

FRA = FRAmap(tifFileList,2,2,'SCALEDfissaFroi');
nDB = numel(FRA.dBlist); nFreq = numel(FRA.freqList);

%% ---------------------------------------------------------------- 2
c11 = FRA.dBFreqMap{1,1};
assert(islogical(c11.sigPkDFF),'sigPkDFF must be logical');
assert(isequal(size(c11.sigPkDFF),[nCell 1]),'sigPkDFF must be nCell x 1');
assert(isequal(size(c11.pkDFF),[nCell 1]),'pkDFF must be nCell x 1');

tgtCol = sub2ind([nDB nFreq], find(FRA.dBlist==TARGET.d), find(FRA.freqList==TARGET.f));
got = FRA.CellSigPkLinDBfreq(TARGET.cell,tgtCol);
% tolerance covers the injected noise; the point is that the map holds the
% PEAK (~0.5) and not its square (0.25)
assert(abs(got-TARGET.amp) < 0.01, ...
    'sig map should hold the peak %.4f, got %.4f (%.4f would be its square)', ...
    TARGET.amp,got,TARGET.amp^2);
assert(abs(got-TARGET.amp^2) > 0.1,'sig map still looks squared');
sigm = ~isnan(FRA.CellSigPkLinDBfreq);
assert(isequal(FRA.CellSigPkLinDBfreq(sigm),FRA.CellPkRespLinDBfreq(sigm)), ...
    'sig map must equal the peak map where significant');
fprintf('PASS 2: sigPkDFF logical; sig map holds peak %.3f, not %.3f\n', ...
    TARGET.amp,TARGET.amp^2); nPass = nPass+1;

%% ---------------------------------------------------------------- 3
% trial-count invariance: duplicate the tif, results must not move
tfl2 = tifFileList; tfl2.map(2) = tfl2.map(1);
copyfile(fullfile(scratch,'AA9999AAAA_00001_00001_Pulses.mat'), ...
         fullfile(scratch,'AA9999AAAA_00002_00001_Pulses.mat'));
tfl2.map(2).name = 'AA9999AAAA_00002_00001.tif';
FRA2 = FRAmap(tfl2,2,2,'SCALEDfissaFroi');
% the significance PATTERN must be identical; peak VALUES only to within
% float tolerance, since averaging 6 duplicated trials rather than 3 changes
% the summation order (observed difference ~2e-19)
assert(isequal(~isnan(FRA.CellSigPkLinDBfreq),~isnan(FRA2.CellSigPkLinDBfreq)), ...
    'doubling the trial count changed which conditions are significant');
assert(max(abs(FRA.CellPkRespLinDBfreq(:)-FRA2.CellPkRespLinDBfreq(:))) < 1e-12, ...
    'doubling the trial count moved the peak map by more than float noise');
assert(size(FRA2.dBFreqMap{1,1}.pkDFFbyTrial,2) == ...
       2*size(FRA.dBFreqMap{1,1}.pkDFFbyTrial,2),'per-trial peaks not retained');
fprintf('PASS 3: peak and significance invariant to trial count\n'); nPass = nPass+1;

%% ---------------------------------------------------------------- 4
% Onset alignment. The aligned trace no longer maps onto a within-pulse frame
% index (the baseline reaches back across the pulse boundary), so check the
% alignment directly: the response was injected one frame AFTER onset, so it
% must land at column onsetCol+1 of the trial-averaged trace, for every pulse
% regardless of which onset group it came from.
onsetCol = FRA.params.onsetCol;
assert(onsetCol == FRA.params.nBaselineFrames+1,'onsetCol must be nBase+1');
cTgt = FRA.dBFreqMap{FRA.dBlist==TARGET.d, FRA.freqList==TARGET.f};
[~,hitCol] = max(cTgt.dFFavg(TARGET.cell,:));
assert(hitCol == onsetCol+1, ...
    'injected response landed at column %d, expected onsetCol+1 = %d', ...
    hitCol,onsetCol+1);
% and tPTrel must advance by exactly one frame per column
dt = diff(cTgt.tPTrel(:,1));
assert(max(abs(dt - 1/fs)) < 1e-9,'aligned time axis is not uniform at 1/fs');
fprintf('PASS 4: onset at column %d (nBase %d); injected response lands at %d\n', ...
    onsetCol,FRA.params.nBaselineFrames,hitCol); nPass = nPass+1;

%% ---------------------------------------------------------------- 5
% drop every trial of one freq/dB pair -> placeholder, no reshape break
dropF = freqs(1); dropD = dBs(1);
keep = cellfun(@(c) ~(c.f==dropF && c.d==dropD), conds);
pulse3 = pulse(keep);
params3 = params; params3.totalPulses = numel(pulse3);
raw3 = baselineF*ones(nCell, stimDelay*fs + framesPerPulse*numel(pulse3));
s3 = fullfile(scratch,'drop'); mkdir(s3);
pulse = pulse3; params = params3; %#ok<NASGU>
save(fullfile(s3,'AA9999AAAA_00003_00001_Pulses.mat'),'pulse','params');
tfl3 = struct('map',struct('name','AA9999AAAA_00003_00001.tif','folder',s3, ...
    'treatment','none FRAmap','nFrames',size(raw3,2),'frameRate',fs, ...
    'SCALEDfissaFroi',raw3));
FRA3 = FRAmap(tfl3,2,2,'SCALEDfissaFroi');
assert(isequal(size(FRA3.CellSigPkLinDBfreq), ...
    [nCell numel(FRA3.dBlist)*numel(FRA3.freqList)]), ...
    'missing-pair fill produced the wrong map size');
fprintf('PASS 5: missing freq/dB pair filled; map is %s\n', ...
    mat2str(size(FRA3.CellSigPkLinDBfreq))); nPass = nPass+1;

%% ---------------------------------------------------------------- 6
% peak window is ceil'd. 300 ms at 5 Hz = 1.5 frames -> 2, not 1.
assert(ceil(300/1000*fs)+2 == 4,'ceil arithmetic wrong');
assert(numel(6:(6+2.5-1)) == 2 && numel(6:(6+ceil(2.5)-1)) == 3, ...
    'a bare colon floors the window; ceil must be applied before it');
[~,~,pkCeil] = pkFcalc([zeros(1,5) 1 2 9 1 1],6,ceil(2.5),2);
[~,~,pkFloor] = pkFcalc([zeros(1,5) 1 2 9 1 1],6,2.5,2);
assert(pkCeil == 9 && pkFloor == 2,'ceil vs floor window did not differ as expected');
fprintf('PASS 6: window ceil'' d (peak %g) vs floored (%g)\n',pkCeil,pkFloor); nPass = nPass+1;

%% ---------------------------------------------------------------- 7
% linear index ordering assumed by anmlFRA2BF / anmlFRA2dPrime
tagged = cell(nDB,nFreq);
for r = 1:nDB, for c = 1:nFreq, tagged{r,c} = repmat(r+100*c,nCell,1); end, end
lin = reshape(cell2mat(tagged),[nCell nDB*nFreq]);
[rr,cc] = ind2sub([nDB nFreq],1:nDB*nFreq);
assert(isequal(lin(1,:),rr+100*cc), ...
    'reshape ordering no longer matches sub2ind([nDB nFreq])');
fprintf('PASS 7: linear index ordering matches sub2ind([nDB nFreq])\n'); nPass = nPass+1;

%% ---------------------------------------------------------------- 8
dp = cfg.exampleAnimalDir;
if isfolder(dp)
    S = load(fullfile(dp,'TO0003_tifFileList.mat'),'tifFileList');
    tfl = S.tifFileList; [tfl.map.folder] = deal(dp);
    R = FRAmap(tfl,2,2,'SCALEDfissaFroi');       % processFRA stock defaults

    assert(strcmp(R.params.sigMethod,'trialAveraged'));
    assert(strcmp(R.params.baseline,'preOnsetExtended'));
    % baseline now reaches back across the pulse boundary into the previous
    % pulse's silent tail, so onset sits at column nBase+1, not at 4
    assert(R.params.nBaselineFrames == 7, ...
        'TO0003 should resolve a 7-frame baseline, got %d',R.params.nBaselineFrames);
    assert(R.params.onsetCol == R.params.nBaselineFrames+1);
    assert(R.params.onsetCol == 8,'TO0003 onset column should be 8');
    assert(R.params.alignedLen == 17);

    sr = mean(cell2mat(cellfun(@(x) x.sigPkDFF(:),R.dBFreqMap(:),'uni',0)));
    assert(abs(sr-0.376) < 0.01,'TO0003 sig rate %.3f, expected ~0.376',sr);

    % the extended baseline must be genuinely silent: with F0 taken from it
    % its mean is 0 by construction, but a previous response leaking in would
    % show as structure across those columns
    G = []; for i=1:numel(R.dBFreqMap), G = cat(3,G,R.dBFreqMap{i}.dFFptRel); end
    gm = squeeze(mean(mean(G,3,'omitnan'),1,'omitnan'));
    assert(max(abs(gm(1:R.params.nBaselineFrames))) < 0.01, ...
        'baseline columns are not flat: %s',mat2str(round(gm(1:R.params.nBaselineFrames),4)));
    [~,pkCol] = max(gm);
    assert(pkCol > R.params.onsetCol, ...
        'grand-mean peak at col %d should follow onset col %d',pkCol,R.params.onsetCol);
    assert(all(any(~isnan(R.CellSigPkLinDBfreq),2)), ...
        'every TO0003 cell should have >=1 significant condition');

    % ROI 1 at 28284 Hz / 70 dB: the case the inclusive baseline rejects
    col = sub2ind([numel(R.dBlist) numel(R.freqList)], ...
        find(R.dBlist==70), find(R.freqList==28284));
    pk1 = R.CellPkRespLinDBfreq(1,col);
    assert(abs(pk1-0.3885) < 1e-3,'ROI 1 peak %.4f, expected 0.3885',pk1);
    assert(~isnan(R.CellSigPkLinDBfreq(1,col)), ...
        ['ROI 1 at 28284 Hz / 70 dB must be significant. If this fails, FRA ' ...
         'has fallen back to pkFcalc''s inclusive baseline.']);
    assert(all(isfinite(R.BFuDB)),'BF must be finite for every cell');
    fprintf('PASS 8: TO0003 sig rate %.3f; ROI 1 @28284Hz/70dB peak %.4f significant\n', ...
        sr,pk1); nPass = nPass+1;
else
    fprintf('SKIP 8: %s not mounted\n',dp);
end

%% ---------------------------------------------------------------- 9
% Extended baseline: it must reach into the PREVIOUS pulse's silent tail
% without ever absorbing that pulse's response.
%
% Inject a large response into every pulse's own response window. If the
% baseline of pulse p leaked into pulse p-1's response, F0 would be inflated
% and the measured dF/F would collapse.
fs2 = 5; ISI2 = 3; stimDelay2 = 4; nCell2 = 3;
fpp2 = ISI2*fs2; onset2 = 0.6; nPulse2 = 8;
onFrameLocal = find((0:fpp2-1)/fs2 >= onset2-1e-9,1,'first');
nFrames2 = stimDelay2*fs2 + fpp2*nPulse2;
rng(11); raw2 = 100*(1+0.001*randn(nCell2,nFrames2));
pulse = struct('pulsename',{});
freqs2 = [8409 25937]; dBs2 = [30 70];
for q = 1:nPulse2
    fq = freqs2(mod(q,2)+1); db = dBs2(mod(floor((q-1)/2),2)+1);
    pulse(q).pulsename = sprintf( ...
        '%dHz_%ddB_TestTone_400msPulse_at_%gs_10msRamp_1500msTotal_Fs250kHz',fq,db,onset2);
    hit = stimDelay2*fs2 + (q-1)*fpp2 + onFrameLocal + 1;
    raw2(:,hit) = 100*2.0;                      % big response in EVERY pulse
end
params = struct('stimDelay',stimDelay2,'totalPulses',nPulse2,'ISI',ISI2, ...
    'PulseHiTime',0,'PulseLoTime',0);
s9 = fullfile(scratch,'base'); mkdir(s9);
save(fullfile(s9,'AA9999AAAA_00009_00001_Pulses.mat'),'pulse','params');
tfl9 = struct('map',struct('name','AA9999AAAA_00009_00001.tif','folder',s9, ...
    'treatment','none FRAmap','nFrames',nFrames2,'frameRate',fs2, ...
    'SCALEDfissaFroi',raw2));
R9 = FRAmap(tfl9,2,2,'SCALEDfissaFroi');
nb = R9.params.nBaselineFrames;

% the baseline must not extend past the previous tone's response
prevGap = fpp2 - (ceil(400/1000*fs2)+2) + onFrameLocal - 1;   % frames of silence
assert(nb <= prevGap, ...
    'baseline of %d frames overruns the %d silent frames before onset',nb,prevGap);
assert(nb >= onFrameLocal-1, 'extended baseline must never be shorter than the old one');

% response amplitude must survive: a leaked previous response would inflate F0
c9 = R9.dBFreqMap{find(~cellfun(@(x) isempty(x.dFFavg),R9.dBFreqMap),1)};
assert(max(c9.pkDFF) > 0.5, ...
    'injected response collapsed to %.3f -- F0 likely absorbed a previous response', ...
    max(c9.pkDFF));

% first pulse of a tif draws its baseline from the stimDelay period
assert(stimDelay2*fs2 + onFrameLocal - nb >= 1, ...
    'first pulse baseline runs before the recording start');
fprintf('PASS 9: baseline %d frames, stays clear of the previous response (%d silent), peak %.2f\n', ...
    nb,prevGap,max(c9.pkDFF)); nPass = nPass+1;

%% ---------------------------------------------------------------- 10
% an explicit nBaselineFrames overrides the auto value
R10 = FRAmap(tfl9,2,2,'SCALEDfissaFroi','nBaselineFrames',3);
assert(R10.params.nBaselineFrames == 3);
assert(R10.params.onsetCol == 4);
assert(R10.params.alignedLen == R9.params.alignedLen - (nb-3), ...
    'aligned length must track the baseline length');
fprintf('PASS 10: nBaselineFrames override honoured (3 -> onsetCol 4)\n'); nPass = nPass+1;

fprintf('\nAll %d test groups passed.\n',nPass);
