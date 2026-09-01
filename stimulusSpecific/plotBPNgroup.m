function out = plotBPNgroup(src,varargin)
% plotBPNgroup  Cohort plots for one BPN condition group.
%
%   out = plotBPNgroup(groupFile)
%   out = plotBPNgroup(T)
%   out = plotBPNgroup(...,'plots',{'rlf','traces','peak'},'nConsec',3)
%
%   Replaces two ad hoc cohort paths: running processBPN2P's Plot 3 cell
%   against a group table, and processRLF's hardcoded list of per-animal
%   files. Works for any group size including a single animal, and states
%   cells/mice on every panel.
%
%   Inputs:
%     src - path to a BPN_Group<g>.mat, or a BPN group table.
%
%   Name/Value:
%     'plots'      - any of {'rlf','traces','peak'}. Default all three.
%     'nConsec'    - min consecutive significant levels for RLF inclusion.
%                    Default 3 (tableRLF's default).
%     'levels'     - dB levels to show. Default: all present, ascending.
%     'traceCells' - 'all' (default) or 'included' (RLF-included cells only)
%                    for the trace and peak panels. 'all' matches the
%                    behaviour of processBPN2P Plot 3.
%     'colors'     - nLevels x 3 colormap. Default jet(nLevels).
%     'ax'         - struct with fields rlf/traces/peak to plot into existing
%                    axes instead of creating figures.
%     'verbose'    - print the per-panel summary. Default true.
%
%   Output (struct) - the numbers behind every panel, so the figures can be
%   regenerated or checked without re-deriving them:
%     .N          groupN for the whole group
%     .Nplot      groupN for the cells actually plotted
%     .levels     dB axis
%     .rlf        tableRLF struct (empty when 'rlf' not requested)
%     .traces     per level: .t .mean .sem .n .showBand .cells
%     .peak       per level: .values .n .mean .sem, plus .cells
%     .fig        figure handles created
%
%   Degenerate groups are handled per etc/harmonization_plan.md: a level with
%   no cells is drawn as a labelled empty panel rather than erroring, a single
%   cell gets no SEM band, and a single-animal group is stamped as carrying no
%   across-animal inference.
%
%   See also loadStimGroup, tableRLF, plotRLF, gatherCellTraces, cohortMeanSEM.

p = inputParser;
addRequired(p,'src');
addParameter(p,'plots',{'rlf','traces','peak'},@(x) iscellstr(x)||isstring(x)); %#ok<ISCLSTR>
addParameter(p,'nConsec',3,@(x) isnumeric(x)&&isscalar(x));
addParameter(p,'levels',[],@(x) isempty(x)||isnumeric(x));
addParameter(p,'traceCells','all',@(x) any(strcmpi(x,{'all','included'})));
addParameter(p,'colors',[],@(x) isempty(x)||isnumeric(x));
addParameter(p,'ax',struct(),@isstruct);
addParameter(p,'verbose',true,@islogical);
parse(p,src,varargin{:});
plots      = lower(cellstr(p.Results.plots));
nConsec    = p.Results.nConsec;
levels     = p.Results.levels;
traceCells = lower(p.Results.traceCells);
colors     = p.Results.colors;
axIn       = p.Results.ax;
verbose    = p.Results.verbose;

spec = stimGroupSpec('BPN');

%% ---- load ----
if istable(src)
    T = src;
else
    T = loadStimGroup(src,'family','BPN');
end

out.N = groupN(T);
if isempty(levels); levels = sort(unique(T.(spec.levelVar)))'; else; levels = sort(levels(:))'; end
out.levels = levels;
out.fig = gobjects(0);
out.rlf = [];
if isempty(colors); colors = jet(numel(levels)); end

if verbose
    fprintf('plotBPNgroup: %s | levels %s\n', out.N.label, mat2str(levels));
end

%% ---- RLF (also decides the 'included' cell set) ----
rlf = tableRLF(T,'nConsec',nConsec,'dBlist',levels);
if any(strcmp(plots,'rlf'))
    out.rlf = rlf;
    ax = getAx(axIn,'rlf');
    if isempty(ax); f = figure('Color','w','Name','BPN RLF'); out.fig(end+1) = f; ax = axes(f); end
    plotRLF(rlf,'ax',ax,'showCells',true);
    if out.N.singleAnimal
        text(ax,0.98,0.02,'1 animal - no across-animal inference','Units','normalized', ...
            'HorizontalAlignment','right','VerticalAlignment','bottom', ...
            'FontSize',8,'Color',[0.6 0.2 0.1]);
    end
    if verbose
        fprintf('  RLF: %d of %d cells included (>=%d consecutive sig levels)\n', ...
            rlf.nIncluded, rlf.nTotal, nConsec);
    end
end

%% ---- which cells feed the trace / peak panels ----
if strcmp(traceCells,'included')
    inclKey = cellKey(rlf.cellInfo(rlf.cellInfo.included,:), spec.idVars);
    keepRow = ismember(cellKey(T,spec.idVars), inclKey);
else
    keepRow = true(height(T),1);
end
out.Nplot = groupN(T,keepRow);
if verbose && strcmp(traceCells,'included')
    fprintf('  trace/peak panels restricted to included cells: %s\n', out.Nplot.label);
end

%% ---- traces per level ----
if any(strcmp(plots,'traces'))
    ax = getAx(axIn,'traces');
    if isempty(ax)
        f = figure('Color','w','Name','BPN dF/F re sound level'); out.fig(end+1) = f; ax = axes(f);
    end
    hold(ax,'on');
    hLine = gobjects(numel(levels),1); lbl = strings(numel(levels),1);
    out.traces = struct('level',{},'t',{},'mean',{},'sem',{},'n',{},'showBand',{},'cells',{});

    for k = 1:numel(levels)
        mask = keepRow & T.(spec.levelVar) == levels(k);
        [M,cells,gi] = gatherCellTraces(T,mask,spec.cellAvgVar,spec.timeVar,'idVars',spec.idVars);
        s = cohortMeanSEM(M);
        out.traces(k) = struct('level',levels(k),'t',gi.t,'mean',s.mean,'sem',s.sem, ...
            'n',s.n,'showBand',s.showBand,'cells',{cells});
        if ~gi.ok; continue; end
        if s.showBand
            [hPatch,hp] = fillSEMplot(gi.t,s.mean,s.sem,colors(k,:),colors(k,:),ax);
            set(hPatch,'FaceAlpha',0.25); set(hp,'LineWidth',1.8);
        else
            plot(ax,gi.t,s.mean,'Color',colors(k,:),'LineWidth',1.8,'LineStyle','--');
        end
        hLine(k) = plot(ax,NaN,NaN,'Color',colors(k,:),'LineWidth',1.8);
        lbl(k) = sprintf('%g dB (n=%d)',levels(k),s.n);
    end

    drawStimMarkers(ax,T);
    xlabel(ax,'time (s)'); ylabel(ax,'\DeltaF/F');
    title(ax,'Population \DeltaF/F re sound level');
    keep = isgraphics(hLine);
    if any(keep); legend(ax,hLine(keep),lbl(keep),'Location','best'); end
    annotateN(ax,out.Nplot,'location','northeast');
    if all([out.traces.n] == 0)
        text(ax,0.5,0.5,'no cells to plot','Units','normalized', ...
            'HorizontalAlignment','center','Color',[0.6 0.2 0.1]);
    end
    grid(ax,'on'); hold(ax,'off');
end

%% ---- peak dF/F vs level ----
if any(strcmp(plots,'peak'))
    ax = getAx(axIn,'peak');
    if isempty(ax)
        f = figure('Color','w','Name','BPN peak dF/F re sound level'); out.fig(end+1) = f; ax = axes(f);
    end
    out.peak = struct('level',{},'values',{},'n',{},'mean',{},'sem',{},'cells',{});
    mu = nan(1,numel(levels)); sem = nan(1,numel(levels));

    for k = 1:numel(levels)
        mask = keepRow & T.(spec.levelVar) == levels(k);
        [v,cells,gi] = gatherCellValues(T,mask,spec.peakVar,'idVars',spec.idVars);
        s = cohortMeanSEM(v(:));          % nCells x 1 -> scalar mean/sem
        out.peak(k) = struct('level',levels(k),'values',v,'n',s.n, ...
            'mean',s.mean,'sem',s.sem,'cells',{cells});
        if ~gi.ok; continue; end
        mu(k) = s.mean; sem(k) = s.sem;
    end

    bar(ax,levels,mu,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
    hold(ax,'on');
    if any(~isnan(sem))
        errorbar(ax,levels,mu,sem,'k','LineStyle','none','LineWidth',1);
    end
    if numel(levels) > 1; dx = 0.15*min(diff(levels)); else; dx = 0.15; end
    rng(0);
    for k = 1:numel(levels)
        pts = out.peak(k).values(:);
        if isempty(pts); continue; end
        jit = (rand(numel(pts),1)*2-1)*dx;
        scatter(ax,levels(k)+jit,pts,26,'MarkerEdgeColor','none', ...
            'MarkerFaceColor',[0 0.2 0.6],'MarkerFaceAlpha',0.5);
    end
    xlabel(ax,'sound level (dB SPL)'); ylabel(ax,'peak \DeltaF/F');
    title(ax,'Peak \DeltaF/F re sound level');
    set(ax,'XTick',levels); box(ax,'on');
    annotateN(ax,out.Nplot);
    hold(ax,'off');
end
end

%% ================= helpers =================
function ax = getAx(axIn,name)
if isfield(axIn,name) && isgraphics(axIn.(name),'axes'); ax = axIn.(name); else; ax = []; end
end

function k = cellKey(T,idVars)
k = string(T.(idVars{1}));
for i = 2:numel(idVars); k = strcat(k,'_',string(T.(idVars{i}))); end
end

function drawStimMarkers(ax,T)
% BPN onset/offset. After combineDiffOnset every row shares these, but check
% rather than assume -- a group pooling animals with different pulse lengths
% would otherwise get a marker that is right for only some of the cells.
on  = unique(T.BPNsOnset);
len = unique(T.BPNmsPulseLen);
if isscalar(on) && isscalar(len)
    xline(ax,on,'--','BPN');
    xline(ax,on+len/1000,'--');
else
    warning('plotBPNgroup:mixedStimTiming', ...
        ['group mixes BPN onsets (%s s) or pulse lengths (%s ms); ' ...
         'stimulus markers omitted'], mat2str(on'), mat2str(len'));
end
end
