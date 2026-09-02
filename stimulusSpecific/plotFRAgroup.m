function out = plotFRAgroup(src,varargin)
% plotFRAgroup  Cohort tuning plots for one FRA condition group.
%
%   out = plotFRAgroup(groupFile)
%   out = plotFRAgroup(T)
%   out = plotFRAgroup(...,'plots',{'bf','threshold','bw','fra'})
%
%   The FRA counterpart of plotBPNgroup / plotCGCgroup. Four panels:
%
%     'bf'        best-frequency distribution, occurrence (%)
%     'threshold' threshold distribution, occurrence (%)
%     'bw'        bandwidth distribution at threshold + one level (BW20),
%                 with bandwidth-vs-level alongside
%     'fra'       frequency response area colormaps: the group-mean FRA plus
%                 a grid of example single-cell receptive fields
%
%   Inputs:
%     src - path to an FRA_Group<g>.mat (built by aggregateStimGroup), or an
%           FRA group table.
%
%   Name/Value:
%     'plots'      - subset of the four above. Default all.
%     'nConsec'    - consecutive significant levels defining threshold.
%                    Default 1.
%     'minBand'    - contiguous significant frequencies required for a level
%                    to count. Default 2. See tableFRAmetrics.
%     'nExample'   - example single-cell FRAs to draw. Default 9.
%     'sigOnly'    - draw the colormaps from significant peaks only (true,
%                    default) or from every peak (false).
%     'bfEdges'    - octave-spaced BF bin edges. Default: the frequency axis
%                    binned into 8.
%     'bwEdges'    - bandwidth bin edges in octaves. Default 0:0.25:2.
%     'colors'     - nLevels x 3 colormap for the per-level lines.
%     'verbose'    - print the panel summary. Default true.
%
%   Output (struct) - the numbers behind every panel:
%     .N         groupN for the group
%     .Nresp     groupN restricted to responsive cells
%     .metrics   tableFRAmetrics output (per-cell BF/threshold/bandwidth)
%     .bf        .edges .counts .pct .values
%     .threshold .levels .counts .pct .values
%     .bw        .edges .counts .pct .values .byLevel
%     .fra       .meanMap .levels .freqs .exampleCells
%     .sham      the noise-floor control from tableFRAmetrics
%     .fig       figure handles created
%
%   READ THE NOISE FLOOR BEFORE INTERPRETING PANELS 2 AND 3. out.sham.ratio
%   is the real-window significance rate divided by the rate the same test
%   gives on a silent late window. A ratio near 1 means the significance mask
%   this group's thresholds and bandwidths are built from carries little
%   stimulus information. The value is printed by 'verbose' and annotated on
%   the threshold panel, because a threshold histogram computed from a mask
%   that cannot distinguish tone from silence looks perfectly reasonable.
%
%   Bandwidth is BW20, not BW10: with levels 20 dB apart, threshold+10 dB is
%   never sampled. See tableFRAmetrics.
%
%   See also loadStimGroup, tableFRAmetrics, FRAmap2table, plotBPNgroup.

p = inputParser;
addRequired(p,'src');
addParameter(p,'plots',{'bf','threshold','bw','fra'},@(x) iscellstr(x)||isstring(x)); %#ok<ISCLSTR>
addParameter(p,'nConsec',1,@(x) isnumeric(x)&&isscalar(x));
addParameter(p,'minBand',2,@(x) isnumeric(x)&&isscalar(x));
addParameter(p,'nExample',9,@(x) isnumeric(x)&&isscalar(x));
addParameter(p,'sigOnly',true,@islogical);
addParameter(p,'bfEdges',[],@(x) isempty(x)||isnumeric(x));
addParameter(p,'bwEdges',0:0.25:2,@isnumeric);
addParameter(p,'colors',[],@(x) isempty(x)||isnumeric(x));
addParameter(p,'verbose',true,@islogical);
parse(p,src,varargin{:});
plots    = cellstr(p.Results.plots);
nExample = p.Results.nExample;
sigOnly  = p.Results.sigOnly;
bwEdges  = p.Results.bwEdges;
verbose  = p.Results.verbose;

%% load
if istable(src)
    T = src; groupLabel = 'table';
else
    T = loadStimGroup(src,'family','FRA','strict',false);
    [~,groupLabel] = fileparts(char(src));
end

M   = tableFRAmetrics(T,'nConsec',p.Results.nConsec,'minBand',p.Results.minBand);
ci  = M.cellInfo;
N   = groupN(T);
resp = ci.responsive;
NrespMask = ismember(T(:,{'animal','roiID'}),ci(resp,{'animal','roiID'}));
Nresp = groupN(T,NrespMask);

levels = M.levels; freqs = M.freqs; nL = numel(levels);
colors = p.Results.colors;
if isempty(colors), colors = parula(max(nL,2)); colors = colors(1:nL,:); end

bfEdges = p.Results.bfEdges;
if isempty(bfEdges)
    bfEdges = 2.^linspace(log2(min(freqs)),log2(max(freqs)),9);
end

out = struct('N',N,'Nresp',Nresp,'metrics',M,'sham',M.shamControl, ...
    'group',groupLabel,'fig',gobjects(0));

shamTxt = 'noise floor: n/a';
if M.shamControl.ok
    shamTxt = sprintf('noise floor: real %.2f / sham %.2f = %.2f x', ...
        M.shamControl.realRate,M.shamControl.shamRate,M.shamControl.ratio);
end

%% ---- 1. best frequency distribution ----
if any(strcmpi(plots,'bf'))
    v = ci.bf(resp);
    cnt = histcounts(v,bfEdges);
    pct = 100*cnt/max(sum(cnt),1);
    out.bf = struct('edges',bfEdges,'counts',cnt,'pct',pct,'values',v);

    f = figure('Name',sprintf('%s | best frequency',groupLabel),'Color','w');
    ctr = sqrt(bfEdges(1:end-1).*bfEdges(2:end));
    bar(ctr,pct,1,'FaceColor',[.25 .45 .7],'EdgeColor','w');
    set(gca,'XScale','log','XTick',round(ctr),'XTickLabelRotation',45);
    xlabel('Best frequency (Hz)'); ylabel('Occurrence (%)');
    title(sprintf('Best frequency distribution\n%s (%d responsive)', ...
        N.label,Nresp.nCells));
    box off; out.fig(end+1) = f;
end

%% ---- 2. threshold distribution ----
if any(strcmpi(plots,'threshold'))
    v = ci.threshold(resp);
    cnt = arrayfun(@(L) nnz(v==L),levels);
    nNR = nnz(~resp);
    pct = 100*[cnt nNR]/max(numel(resp),1);
    out.threshold = struct('levels',levels,'counts',[cnt nNR],'pct',pct,'values',v);

    f = figure('Name',sprintf('%s | threshold',groupLabel),'Color','w');
    bar(1:nL+1,pct,0.7,'FaceColor',[.85 .45 .2],'EdgeColor','w');
    set(gca,'XTick',1:nL+1,'XTickLabel',[string(levels) "NR"]);
    xlabel('Threshold (dB SPL)'); ylabel('Occurrence (%)');
    title(sprintf('Threshold distribution\n%s | %s',N.label,shamTxt));
    box off; out.fig(end+1) = f;
end

%% ---- 3. bandwidth ----
if any(strcmpi(plots,'bw'))
    v = ci.bw20(~isnan(ci.bw20));
    cnt = histcounts(v,bwEdges);
    pct = 100*cnt/max(sum(cnt),1);
    out.bw = struct('edges',bwEdges,'counts',cnt,'pct',pct,'values',v, ...
        'byLevel',M.bwByLevel,'refLabel',M.bwRefLabel);

    f = figure('Name',sprintf('%s | bandwidth',groupLabel),'Color','w', ...
        'Position',[0 0 900 380]);
    subplot(1,2,1)
    ctr = bwEdges(1:end-1)+diff(bwEdges)/2;
    bar(ctr,pct,1,'FaceColor',[.35 .6 .35],'EdgeColor','w');
    xlabel('BW20 (octaves)'); ylabel('Occurrence (%)');
    title(sprintf('Bandwidth at threshold +20 dB\nn = %d cells (not BW10: 20 dB level steps)',numel(v)));
    box off
    subplot(1,2,2); hold on
    for li = 1:nL
        x = M.bwByLevel(:,li); x = x(~isnan(x));
        if isempty(x), continue; end
        plot(li+0.12*randn(numel(x),1),x,'.','Color',colors(li,:),'MarkerSize',8);
        plot(li+[-.25 .25],median(x)*[1 1],'k-','LineWidth',2);
    end
    set(gca,'XTick',1:nL,'XTickLabel',string(levels),'XLim',[.5 nL+.5]);
    xlabel('Level (dB SPL)'); ylabel('Bandwidth (octaves)');
    title('Bandwidth by level (median bar)'); box off
    out.fig(end+1) = f;
end

%% ---- 4. FRA colormaps ----
if any(strcmpi(plots,'fra'))
    src3 = M.FRA; if ~sigOnly, src3 = M.FRAall; end
    meanMap = squeeze(mean(src3,1,'omitnan'));
    if nL==1, meanMap = reshape(meanMap,1,[]); end

    [~,ord] = sort(ci.peakDFF,'descend','MissingPlacement','last');
    ex = ord(1:min(nExample,nnz(resp)))';
    out.fra = struct('meanMap',meanMap,'levels',levels,'freqs',freqs, ...
        'exampleCells',ci(ex,{'animal','roiID'}),'sigOnly',sigOnly);

    f = figure('Name',sprintf('%s | FRA maps',groupLabel),'Color','w', ...
        'Position',[0 0 1200 760]);
    subplot(3,4,[1 2 5 6])
    localFRAimage(meanMap,levels,freqs);
    c = colorbar; c.Label.String = '\DeltaF/F';
    title(sprintf('Group mean FRA (%s)\n%s',ternary(sigOnly,'significant peaks','all peaks'),N.label));

    tiles = [3 4 7 8 9 10 11 12];
    for k = 1:min(numel(ex),numel(tiles))
        subplot(3,4,tiles(k))
        localFRAimage(squeeze(src3(ex(k),:,:)),levels,freqs);
        title(sprintf('%s ROI %s | BF %.1f kHz', ...
            ci.animal(ex(k)),ci.roiID(ex(k)),ci.bf(ex(k))/1000),'FontSize',8);
    end
    out.fig(end+1) = f;
end

%% ---- report ----
if verbose
    fprintf('\n%s\n',groupLabel);
    fprintf('  %s\n',N.label);
    fprintf('  responsive: %d of %d cells\n',Nresp.nCells,N.nCells);
    fprintf('  %s\n',shamTxt);
    if M.shamControl.ok && M.shamControl.ratio < 1.2
        fprintf(['  WARNING: significance mask barely beats its own noise floor;\n' ...
                 '           threshold and bandwidth panels are largely uninformative.\n']);
    end
    if isfield(out,'threshold')
        fprintf('  threshold: %s\n', strjoin(arrayfun(@(k) ...
            sprintf('%d dB n=%d',levels(k),out.threshold.counts(k)),1:nL,'uni',0),', '));
    end
    if isfield(out,'bw')
        fprintf('  BW20 median %.3f oct (n = %d)\n',median(out.bw.values),numel(out.bw.values));
    end
    if N.singleAnimal
        fprintf('  NOTE: single-animal group -- no across-animal inference.\n');
    end
end
end %function

function localFRAimage(map,levels,freqs)
imagesc(flipud(map));
set(gca,'YTick',1:numel(levels),'YTickLabel',string(flipud(levels(:))));
xt = unique(round(linspace(1,numel(freqs),6)));
set(gca,'XTick',xt,'XTickLabel',compose('%.0f',freqs(xt)/1000));
xlabel('Frequency (kHz)'); ylabel('dB SPL');
end

function v = ternary(c,a,b)
if c, v = a; else, v = b; end
end
