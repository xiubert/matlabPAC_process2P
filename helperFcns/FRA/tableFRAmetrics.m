function out = tableFRAmetrics(T,varargin)
% tableFRAmetrics  Per-cell BF, threshold and bandwidth from an FRA table.
%
%   out = tableFRAmetrics(T)
%   out = tableFRAmetrics(T,'nConsec',1,'bwRefStep',1)
%
%   The FRA analogue of tableRLF: collapses a long-form (ROI x freq x dB)
%   table into one row per cell carrying the tuning descriptors the cohort
%   figures plot.
%
%   LEVEL RESOLUTION -- READ THIS BEFORE INTERPRETING BANDWIDTH.
%   This dataset samples three levels, 30/50/70 dB, 20 dB apart. The
%   conventional "BW10" (bandwidth 10 dB above threshold) is therefore NOT
%   measurable: threshold+10 dB is never a level that was presented. What is
%   computed here is bandwidth at threshold + ONE level, i.e. +20 dB with this
%   sampling, reported as `bw20` and labelled BW20 on every plot. It is the
%   nearest honest analogue, not BW10, and only exists for cells whose
%   threshold is below the top level. `bwByLevel` additionally gives the
%   bandwidth at every level so bandwidth-vs-level is visible rather than
%   collapsed into one number.
%
%   Inputs:
%     T - FRA group or per-animal table with columns animal, roiID, freqHz,
%         dBSPL, pkDFF, sig (see FRAmap2table).
%
%   Name/Value:
%     'cellIDvars' - columns identifying a unique cell. Default
%                    {'animal','roiID'}.
%     'nConsec'    - a cell's threshold is the lowest level from which it is
%                    significant on nConsec CONSECUTIVE levels. Default 1.
%                    tableRLF defaults to 3, but with only three levels here
%                    nConsec 3 would demand significance at every level; 1
%                    accepts any significant level as threshold and 2 matches
%                    the RLF group files (nConsec2).
%     'bwRefStep'  - how many levels above threshold to measure bandwidth at.
%                    Default 1 (= +20 dB with this sampling).
%     'bwDrop'     - how bandwidth is delimited, see below. Default
%                    'contiguous'.
%     'minBand'    - contiguous significant frequencies required before a
%                    level counts toward threshold. Default 3 (0.25 oct at
%                    0.125 oct spacing). With 28 frequencies tested per level
%                    a single isolated significant frequency is common by
%                    chance, so requiring a band is what keeps `threshold`
%                    from collapsing to the lowest level for every cell.
%                    A cell with no level meeting minBand falls back to
%                    any-significant-frequency so it still gets a threshold.
%
%   NOISE FLOOR -- out.shamControl. The per-condition significance test in
%   FRAmap compares a peak over 4 frames against a mean+2SD threshold built
%   from only 3 pre-onset frames, which is a weak discriminator. shamControl
%   re-runs exactly that test on a silent LATE window of the same width, where
%   no tone-evoked response can exist, and reports the rate. Read it as the
%   false-positive floor: .ratio near 1 means this group's significance mask
%   carries little stimulus information and its threshold and bandwidth
%   distributions should not be over-interpreted. Measured per animal it
%   ranges from ~1.6 (TO0001, TO0003) to ~0.8 (TO0010).
%
%   BANDWIDTH DEFINITION. At a given level, take that cell's significant
%   frequencies. 'contiguous' walks outward from the strongest significant
%   frequency at that level and stops at the first non-significant frequency
%   on each side; bandwidth is log2(fHigh/fLow) octaves across that run.
%   'total' instead spans the lowest to the highest significant frequency,
%   ignoring gaps. A cell significant at exactly one frequency has bandwidth
%   0 under both, not one frequency step.
%
%   BEST FREQUENCY. `bf` matches anmlFRA2BF: mean the significant peak across
%   levels per frequency, then take the argmax, so it is the tuning centroid
%   across the whole FRA. `bfAtThreshold` is instead the strongest significant
%   frequency at the cell's own threshold -- the "tip" of the tuning curve,
%   which is what the tuning-curve schematic in the literature marks. They
%   differ for cells whose tuning shifts with level.
%
%   Output (struct):
%     .cellInfo   table, one row per cell: cellID vars, bf, bfAtThreshold,
%                 threshold, thresholdIdx, bw20, bwAtThreshold, nSigCond,
%                 responsive, peakDFF
%     .bwByLevel  nCell x nLevel octaves, bandwidth at each level (NaN where
%                 the cell has no significant response at that level)
%     .sigByLevel nCell x nLevel count of significant frequencies per level
%     .FRA        nCell x nLevel x nFreq peak dF/F (NaN where not significant)
%     .FRAall     nCell x nLevel x nFreq peak dF/F regardless of significance
%     .levels .freqs .nCell .nResponsive .nConsec .bwRefStep .bwDrop
%
%   Cells with no significant condition anywhere are kept, with responsive
%   false and NaN descriptors, so counts of non-responding cells stay visible
%   (the ABR- cohort in the reference figure is exactly that case).
%
%   See also FRAmap2table, plotFRAgroup, tableRLF, anmlFRA2BF.

p = inputParser;
addRequired(p,'T',@istable);
addParameter(p,'cellIDvars',{'animal','roiID'},@(x) iscellstr(x)||isstring(x)); %#ok<ISCLSTR>
addParameter(p,'nConsec',1,@(x) isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'bwRefStep',1,@(x) isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'bwDrop','contiguous',@(x) any(strcmpi(x,{'contiguous','total'})));
addParameter(p,'minBand',3,@(x) isnumeric(x)&&isscalar(x)&&x>=1);
parse(p,T,varargin{:});
cellIDvars = cellstr(p.Results.cellIDvars);
nConsec    = p.Results.nConsec;
bwRefStep  = p.Results.bwRefStep;
bwDrop     = lower(p.Results.bwDrop);
minBand    = p.Results.minBand;

req = [cellIDvars,{'freqHz','dBSPL','pkDFF','sig'}];
missing = setdiff(req,T.Properties.VariableNames);
if ~isempty(missing)
    error('tableFRAmetrics:missingVars','Table is missing: %s',strjoin(missing,', '));
end

levels = reshape(unique(T.dBSPL),1,[]);          % ascending
freqs  = reshape(unique(T.freqHz),1,[]);
nL = numel(levels); nF = numel(freqs);

[gID,gT] = findgroups(T(:,cellIDvars));
nCell = height(gT);

FRAall = nan(nCell,nL,nF);
sigM   = false(nCell,nL,nF);
for r = 1:height(T)
    li = find(levels==T.dBSPL(r),1);
    fi = find(freqs ==T.freqHz(r),1);
    FRAall(gID(r),li,fi) = T.pkDFF(r);
    sigM(gID(r),li,fi)   = logical(T.sig(r));
end
FRA = FRAall; FRA(~sigM) = NaN;

[bf,bfAtThreshold,threshold,bw20,bwAtThreshold,peakDFF] = deal(nan(nCell,1));
thresholdIdx = nan(nCell,1);
nSigCond   = zeros(nCell,1);
responsive = false(nCell,1);
bwByLevel  = nan(nCell,nL);
sigByLevel = zeros(nCell,nL);

for c = 1:nCell
    sc = squeeze(sigM(c,:,:)); if nL==1, sc = reshape(sc,1,nF); end
    ac = squeeze(FRA(c,:,:));  if nL==1, ac = reshape(ac,1,nF); end
    sigByLevel(c,:) = sum(sc,2)';
    nSigCond(c) = sum(sc(:));
    responsive(c) = nSigCond(c) > 0;

    % bandwidth at every level
    for li = 1:nL
        bwByLevel(c,li) = localBW(sc(li,:),ac(li,:),freqs,bwDrop);
    end
    if ~responsive(c), continue; end

    peakDFF(c) = max(ac(:),[],'omitnan');

    % BF across levels, matching anmlFRA2BF (mean sig peak per freq, argmax)
    perFreq = mean(ac,1,'omitnan');
    [~,bi] = max(perFreq);
    bf(c) = freqs(bi);

    % threshold: lowest level starting a run of nConsec significant levels.
    % "significant at a level" requires a CONTIGUOUS band of minBand
    % frequencies, not a single isolated one: with 28 frequencies tested per
    % level, isolated significant frequencies are common by chance (see the
    % sham control in out.shamControl).
    anySig = false(1,nL);
    for li = 1:nL
        anySig(li) = localMaxRun(sc(li,:)) >= minBand;
    end
    if ~any(anySig), anySig = sigByLevel(c,:) > 0; end
    ti = NaN;
    for li = 1:nL-nConsec+1
        if all(anySig(li:li+nConsec-1)), ti = li; break; end
    end
    if isnan(ti)
        % no run long enough; fall back to the lowest significant level so a
        % responsive cell always has a threshold, and record that below
        ti = find(anySig,1,'first');
    end
    threshold(c)    = levels(ti);
    thresholdIdx(c) = ti;

    % tuning tip at threshold
    [~,fiT] = max(ac(ti,:),[],'omitnan');
    if any(sc(ti,:)), bfAtThreshold(c) = freqs(fiT); end
    bwAtThreshold(c) = bwByLevel(c,ti);

    % bandwidth bwRefStep levels above threshold (+20 dB with 3 levels)
    ri = ti + bwRefStep;
    if ri <= nL, bw20(c) = bwByLevel(c,ri); end
end

cellInfo = gT;
cellInfo.bf            = bf;
cellInfo.bfAtThreshold = bfAtThreshold;
cellInfo.threshold     = threshold;
cellInfo.thresholdIdx  = thresholdIdx;
cellInfo.bw20          = bw20;
cellInfo.bwAtThreshold = bwAtThreshold;
cellInfo.nSigCond      = nSigCond;
cellInfo.responsive    = responsive;
cellInfo.peakDFF       = peakDFF;

sham = localShamControl(T,cellIDvars);

out = struct('cellInfo',cellInfo,'bwByLevel',bwByLevel,'sigByLevel',sigByLevel, ...
    'shamControl',sham,'minBand',minBand, ...
    'FRA',FRA,'FRAall',FRAall,'levels',levels,'freqs',freqs, ...
    'nCell',nCell,'nResponsive',nnz(responsive),'nConsec',nConsec, ...
    'bwRefStep',bwRefStep,'bwDrop',bwDrop, ...
    'bwRefLabel',sprintf('BW%+d dB re threshold',bwRefStep*median(diff(levels))));
end %function

function n = localMaxRun(m)
%longest contiguous run of true in a logical row
d = diff([0 reshape(logical(m),1,[]) 0]);
n = max([0, find(d==-1)-find(d==1)]);
end

function s = localShamControl(T,cellIDvars)
%Re-run FRAmap's significance rule on a silent late window of the same width.
%Any excess of the real rate over this sham rate is the stimulus-driven part.
s = struct('realRate',NaN,'shamRate',NaN,'ratio',NaN,'nWin',NaN,'ok',false);
if ~all(ismember({'dFFavg','t_PTrel'},T.Properties.VariableNames)); return; end
w = cellfun(@numel,T.dFFavg);
if numel(unique(w))~=1; return; end            % ragged traces: skip
A = cell2mat(cellfun(@(c) reshape(c,1,[]),T.dFFavg,'uni',0));
t = reshape(T.t_PTrel{1},1,[]);
on = find(t >= -1e-9,1,'first');               % onset column (t == 0)
if isempty(on) || on < 2; return; end
nWin = min(4,size(A,2)-on+1);
shamStart = size(A,2)-nWin+1;
if shamStart <= on; return; end                % no silent window available
base = A(:,1:on-1);
thr  = mean(base,2,'omitnan') + 2*std(base,0,2,'omitnan');
pkR  = max(A(:,on:on+nWin-1),[],2);
pkS  = max(A(:,shamStart:shamStart+nWin-1),[],2);
s.realRate = mean(pkR>=thr,'omitnan');
s.shamRate = mean(pkS>=thr,'omitnan');
s.ratio    = s.realRate/s.shamRate;
s.nWin     = nWin;
s.ok       = true;
end

function bw = localBW(sigRow,ampRow,freqs,mode)
%octaves spanned by the significant frequencies at one level
bw = NaN;
if ~any(sigRow), return; end
idx = find(sigRow);
if strcmp(mode,'total')
    lo = idx(1); hi = idx(end);
else
    % contiguous run containing the strongest significant frequency
    a = ampRow; a(~sigRow) = -Inf;
    [~,pk] = max(a);
    lo = pk; while lo>1 && sigRow(lo-1), lo = lo-1; end
    hi = pk; while hi<numel(sigRow) && sigRow(hi+1), hi = hi+1; end
end
bw = log2(freqs(hi)/freqs(lo));
end
