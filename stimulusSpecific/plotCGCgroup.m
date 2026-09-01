function out = plotCGCgroup(src,varargin)
% plotCGCgroup  Cohort plots for one CGC (contrast gain control) group.
%
%   out = plotCGCgroup(groupFile)
%   out = plotCGCgroup(T)
%   out = plotCGCgroup(...,'plots',{'traces','scatter','bar'},'sigOnly',true)
%
%   Replaces hand-running processCGC's population cells against a group
%   table. Works for any group size including a single animal, and cannot
%   crash on a group where no cell passes the significance filter -- which is
%   the case for 4 of the 13 TOMT animals taken individually.
%
%   Inputs:
%     src - path to a CGC_Group<g>.mat, or a CGC group table.
%
%   Name/Value:
%     'plots'       - any of {'traces','scatter','bar'}. Default all three.
%     'sigOnly'     - restrict panels to cells significant in EVERY contrast
%                     (processCGC's `valid = all(sigByROI,2)`). Default true.
%                     The scatter and bar always use this set; sigOnly=false
%                     relaxes only the trace panel, matching processCGC's
%                     popTraceSigOnly toggle.
%     'traceXlim'   - x-limits (s) for the trace panel. Default [1 5].
%     'scatterLim'  - axis limits for the low-vs-high scatter. Default [0 1].
%     'minN'        - minimum cells for the paired test. Default 3.
%     'showFit'     - draw a linear fit + 95% CI band on the scatter when
%                     there are enough points. Default true.
%     'fitMinN'     - minimum cells to attempt the fit. Default 3.
%     'fitColors'   - 2x3: row 1 fit line, row 2 CI fill. Default black line,
%                     grey band (regPlot's convention in plotDataTable).
%     'colors'      - struct from getContrastColors, or nLevels x 3.
%     'ax'          - struct with fields traces/scatter/bar to plot into.
%     'legendLocation' - legend placement for the trace panel. Default
%                     'northeast'. Never 'best' -- exportgraphics recomputes
%                     those and relocates them across tiles.
%     'verbose'     - print a summary. Default true.
%
%   Output (struct):
%     .N        groupN for the whole group
%     .Nplot    groupN for the cells actually plotted
%     .levels   dBdelta axis, ascending (low contrast first)
%     .cells    table of cell identities, rows aligned with .pk/.sig/.valid
%     .pk       nCells x nLevels peak dF/F (NaN where a contrast is missing)
%     .sig      nCells x nLevels logical significance
%     .valid    nCells x 1: present AND significant in every contrast
%     .traces   per level: .level .label .t .mean .sem .n .showBand
%     .stat     cohortStat result for the paired low-vs-high comparison
%     .fit      linear fit of high vs low peak: .ok .reason .slope .slopeCI
%               .n .r2 .mdl. The model is fit THROUGH THE ORIGIN
%               (y ~ x - 1), matching plotDataTable's
%               regPlot(...,'intercept',false), so the slope is directly the
%               high/low gain ratio and a CI excluding 1 means the contrasts
%               differ.
%     .fig      figure handles created
%
%   dBdelta is the dB RANGE of the DRC, so the smaller value is LOW contrast
%   (e.g. 50-60 dB = 10) and the larger is HIGH contrast (40-70 dB = 30).
%   The scatter puts LOW on x and HIGH on y, as plotDataTable.m does.
%
%   See also loadStimGroup, gatherCellTraces, gatherCellValues, cohortStat,
%   processCGC, getContrastColors.

p = inputParser;
addRequired(p,'src');
addParameter(p,'plots',{'traces','scatter','bar'},@(x) iscellstr(x)||isstring(x)); %#ok<ISCLSTR>
addParameter(p,'sigOnly',true,@islogical);
addParameter(p,'traceXlim',[1 5],@(x) isnumeric(x)&&numel(x)==2);
addParameter(p,'scatterLim',[0 1],@(x) isnumeric(x)&&numel(x)==2);
addParameter(p,'minN',3,@(x) isnumeric(x)&&isscalar(x));
addParameter(p,'showFit',true,@islogical);
addParameter(p,'fitMinN',3,@(x) isnumeric(x)&&isscalar(x));
addParameter(p,'fitColors',[0 0 0; 0.651 0.651 0.651],@(x) isequal(size(x),[2 3]));
addParameter(p,'colors',[],@(x) isempty(x)||isstruct(x)||isnumeric(x));
addParameter(p,'ax',struct(),@isstruct);
addParameter(p,'legendLocation','northeast',@(x) ischar(x)||isstring(x));
addParameter(p,'verbose',true,@islogical);
parse(p,src,varargin{:});
plots      = lower(cellstr(p.Results.plots));
sigOnly    = p.Results.sigOnly;
traceXlim  = p.Results.traceXlim;
scatterLim = p.Results.scatterLim;
minN       = p.Results.minN;
showFit    = p.Results.showFit;
fitMinN    = p.Results.fitMinN;
fitColors  = p.Results.fitColors;
colors     = p.Results.colors;
axIn       = p.Results.ax;
legendLoc  = char(p.Results.legendLocation);
verbose    = p.Results.verbose;

spec = stimGroupSpec('CGC');

%% ---- load ----
if istable(src)
    T = src;
else
    T = loadStimGroup(src,'family','CGC');
end

out.N   = groupN(T);
out.fig = gobjects(0);
levels  = sort(unique(T.(spec.levelVar)))';
out.levels = levels;
nL = numel(levels);
labels = contrastLabels(levels);

if isempty(colors); colors = getContrastColors(); end
[lineCol,bandCol] = resolveColors(colors,nL);

%% ---- per-cell peak + significance matrices ----
% Built by explicit cell-key lookup so nothing depends on table row order,
% and a cell missing a contrast stays NaN rather than shifting the matrix.
idv     = spec.idVars;
allKeys = unique(cellKey(T,idv),'stable');
nCell   = numel(allKeys);
pk  = nan(nCell,nL);
sg  = false(nCell,nL);

for k = 1:nL
    mask = T.(spec.levelVar) == levels(k);
    [v,cellsK,gi] = gatherCellValues(T,mask,spec.peakVar,'idVars',idv);
    if ~gi.ok; continue; end
    s = gatherCellValues(T,mask,spec.sigVar,'idVars',idv,'reduce','any');
    [tf,loc] = ismember(cellKey(cellsK,idv),allKeys);
    pk(loc(tf),k) = v(tf);
    sg(loc(tf),k) = s(tf) == 1;
end

out.cells = firstRowsFor(T,idv,allKeys);
out.pk    = pk;
out.sig   = sg;
out.valid = all(sg,2) & all(~isnan(pk),2);   % present AND significant everywhere

validKeys = allKeys(out.valid);
if sigOnly
    keepRow = ismember(cellKey(T,idv),validKeys);
else
    keepRow = true(height(T),1);
end
out.Nplot = groupN(T,keepRow);

if verbose
    fprintf('plotCGCgroup: %s | contrasts %s\n', out.N.label, mat2str(levels));
    fprintf('  significant in every contrast: %d of %d cells\n', sum(out.valid), nCell);
    if sum(out.valid) == 0
        fprintf('  NOTE: no cells passed; panels will be labelled empty.\n');
    end
end

%% ---- traces ----
if any(strcmp(plots,'traces'))
    ax = getAx(axIn,'traces');
    if isempty(ax)
        f = figure('Color','w','Name','CGC dF/F re contrast'); out.fig(end+1) = f; ax = axes(f);
    end
    hold(ax,'on');
    hLine = gobjects(nL,1); leg = strings(nL,1);
    out.traces = struct('level',{},'label',{},'t',{},'mean',{},'sem',{},'n',{},'showBand',{});

    for k = 1:nL
        mask = keepRow & T.(spec.levelVar) == levels(k);
        [M,~,gi] = gatherCellTraces(T,mask,spec.cellAvgVar,spec.timeVar,'idVars',idv);
        s = cohortMeanSEM(M);
        out.traces(k) = struct('level',levels(k),'label',labels(k),'t',gi.t, ...
            'mean',s.mean,'sem',s.sem,'n',s.n,'showBand',s.showBand);
        if ~gi.ok; continue; end
        if s.showBand
            fillSEMplot(gi.t,s.mean,s.sem,lineCol(k,:),bandCol(k,:),ax);
        else
            plot(ax,gi.t,s.mean,'Color',lineCol(k,:),'LineWidth',1.8,'LineStyle','--');
        end
        hLine(k) = plot(ax,NaN,NaN,'Color',lineCol(k,:),'LineWidth',1.8);
        leg(k) = sprintf('%s (n=%d)',labels(k),s.n);
    end

    drawPTmarker(ax,T);
    xlabel(ax,'time (s)'); ylabel(ax,'\DeltaF/F');
    title(ax,'Population \DeltaF/F re contrast');
    if ~isempty(traceXlim); xlim(ax,traceXlim); end
    % Explicit location, never 'best' -- see plotBPNgroup for why.
    keep = isgraphics(hLine);
    if any(keep); legend(ax,hLine(keep),leg(keep),'Location',legendLoc,'Box','off'); end
    annotateN(ax,out.Nplot,'location','southeast', ...
        'extra',ternary(sigOnly,'sig. in every contrast','all cells'));
    if all([out.traces.n] == 0)
        text(ax,0.5,0.5,'no cells passed the significance filter', ...
            'Units','normalized','HorizontalAlignment','center','Color',[0.6 0.2 0.1]);
    end
    hold(ax,'off');
end

%% ---- low vs high scatter ----
if any(strcmp(plots,'scatter'))
    ax = getAx(axIn,'scatter');
    if isempty(ax)
        f = figure('Color','w','Name','CGC peak low vs high'); out.fig(end+1) = f; ax = axes(f);
    end
    if nL ~= 2
        text(ax,0.5,0.5,sprintf('scatter needs exactly 2 contrasts (found %d)',nL), ...
            'Units','normalized','HorizontalAlignment','center','Color',[0.6 0.2 0.1]);
    else
        % LOW contrast on x, HIGH on y -- matches
        % matlabPAC_CGCplot/plotDataTable.m's regPlot(lo,hi,...) call.
        xLo = pk(out.valid,1);   % low contrast  (smaller dBdelta)
        yHi = pk(out.valid,2);   % high contrast (larger dBdelta)
        hold(ax,'on');

        % Fit first so the CI band sits behind the points.
        out.fit = fitLowHigh(ax,xLo,yHi,scatterLim,fitMinN,fitColors,showFit);

        plot(ax,scatterLim,scatterLim,'--','Color',[0.4 0.4 0.4], ...
            'LineWidth',1,'HandleVisibility','off');
        scatter(ax,xLo,yHi,45,'o','MarkerEdgeColor',[0.45 0.45 0.45], ...
            'MarkerFaceColor','none','LineWidth',1);

        xlim(ax,scatterLim); ylim(ax,scatterLim); axis(ax,'square');
        % Axis labels carry the contrast colour, as in plotDataTable.
        xlabel(ax,{char(labels(1)),'peak \DeltaF/F'},'Color',lineCol(1,:));
        ylabel(ax,{char(labels(2)),'peak \DeltaF/F'},'Color',lineCol(2,:));

        if out.fit.ok
            title(ax,sprintf('Peak \\DeltaF/F per cell (slope %.2f, CI %.2f-%.2f)', ...
                out.fit.slope, out.fit.slopeCI(1), out.fit.slopeCI(2)));
        else
            title(ax,'Peak \DeltaF/F per cell');
        end
        if isempty(xLo)
            text(ax,0.5,0.5,'no cells passed the significance filter', ...
                'Units','normalized','HorizontalAlignment','center','Color',[0.6 0.2 0.1]);
        end
        hold(ax,'off');
    end
    annotateN(ax,out.Nplot,'location','southeast','extra','sig. in every contrast');
end

%% ---- paired bar + guarded statistics ----
out.stat = struct('ok',false,'reason','not requested','p',NaN);
if ~isfield(out,'fit')
    out.fit = struct('ok',false,'reason','scatter not requested', ...
        'slope',NaN,'slopeCI',[NaN NaN],'n',0,'r2',NaN);
end
if any(strcmp(plots,'bar'))
    ax = getAx(axIn,'bar');
    if isempty(ax)
        f = figure('Color','w','Name','CGC peak per contrast'); out.fig(end+1) = f; ax = axes(f);
    end
    V = pk(out.valid,:);
    mu = nan(1,nL); sem = nan(1,nL);
    for k = 1:nL
        s = cohortMeanSEM(V(:,k));
        mu(k) = s.mean; sem(k) = s.sem;
    end

    b = bar(ax,1:nL,mu,'FaceColor','flat','EdgeColor','none'); hold(ax,'on');
    for k = 1:nL; b.CData(k,:) = lineCol(k,:); end
    if any(~isnan(sem))
        errorbar(ax,1:nL,mu,sem,'k.','LineWidth',1);
    end
    rng(0);
    for k = 1:nL
        pts = V(:,k); pts = pts(~isnan(pts));
        if isempty(pts); continue; end
        scatter(ax,k+(rand(numel(pts),1)-0.5)*0.16,pts,20,'k','filled','MarkerFaceAlpha',0.6);
    end
    set(ax,'XTick',1:nL,'XTickLabel',cellstr(labels));
    ylabel(ax,'peak \DeltaF/F'); box(ax,'on');

    % paired test, guarded: refuses below minN rather than emitting NaN
    if nL == 2
        out.stat = cohortStat(V(:,1),V(:,2),'paired',true,'minN',minN);
        if out.stat.ok
            annotateSig(ax,mu,sem,out.stat);
            title(ax,sprintf('%s vs %s: %s p = %.3g (n = %d)', ...
                labels(1),labels(2),out.stat.test,out.stat.p,out.stat.n));
        else
            title(ax,sprintf('%s vs %s: test not run',labels(1),labels(2)));
            text(ax,0.5,0.92,out.stat.reason,'Units','normalized', ...
                'HorizontalAlignment','center','Color',[0.6 0.2 0.1],'FontSize',9);
        end
    end
    annotateN(ax,out.Nplot,'extra','sig. in every contrast');
    hold(ax,'off');
end

if verbose && isfield(out.stat,'ok')
    if out.stat.ok
        fprintf('  paired %s: p = %.4g (n = %d cells)\n', out.stat.test, out.stat.p, out.stat.n);
    else
        fprintf('  paired test not run: %s\n', out.stat.reason);
    end
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

function C = firstRowsFor(T,idVars,keys)
k = cellKey(T,idVars);
[~,loc] = ismember(keys,k);
C = T(loc,idVars);
end

function lab = contrastLabels(levels)
lab = strings(numel(levels),1);
for i = 1:numel(levels); lab(i) = sprintf('%g dB range',levels(i)); end
if numel(levels) == 2
    lab(1) = "Low contrast"; lab(2) = "High contrast";
end
end

function [lineCol,bandCol] = resolveColors(colors,nL)
if isnumeric(colors) && size(colors,1) >= nL
    lineCol = colors(1:nL,:); bandCol = lineCol;
elseif isstruct(colors) && isfield(colors,'lohiPre') && nL == 2
    lineCol = colors.lohiPre; bandCol = colors.lohiTracePre;
else
    lineCol = jet(nL); bandCol = lineCol;
end
end

function drawPTmarker(ax,T)
on = unique(T.PTsOnset);
if isscalar(on)
    xline(ax,on,'--','pure tone');
else
    warning('plotCGCgroup:mixedPTonset', ...
        'group mixes PT onsets (%s s); marker omitted', mat2str(on'));
end
end

function annotateSig(ax,mu,sem,stat)
yTop = max(mu+sem);
if ~isfinite(yTop); yTop = max(mu); end
r = range([mu sem]); if r == 0 || ~isfinite(r); r = max(abs(mu))*0.2 + eps; end
yb = yTop + 0.5*r;
plot(ax,[1 2],[yb yb],'-k','LineWidth',1);
plot(ax,[1 1],[yb-0.1*r yb],'-k','LineWidth',1);
plot(ax,[2 2],[yb-0.1*r yb],'-k','LineWidth',1);
text(ax,1.5,yb+0.05*r,stat.stars,'HorizontalAlignment','center','FontSize',14);
end

function v = ternary(c,a,b)
if c; v = a; else; v = b; end
end

function fit = fitLowHigh(ax,x,y,lims,minN,col,doDraw)
% Linear fit of high-contrast peak on low-contrast peak, THROUGH THE ORIGIN,
% with a 95% CI band -- the model plotDataTable.m fits via
% regPlot(lo,hi,'intercept',false). With no intercept the slope IS the
% high/low gain ratio, so a CI excluding 1 means the two contrasts differ.
fit = struct('ok',false,'reason','','slope',NaN,'slopeCI',[NaN NaN], ...
    'n',0,'r2',NaN,'mdl',[]);
keep = ~isnan(x) & ~isnan(y);
x = x(keep); y = y(keep);
fit.n = numel(x);
if ~doDraw
    fit.reason = 'showFit is false'; return
end
if fit.n < minN
    fit.reason = sprintf('n = %d, below fitMinN = %d', fit.n, minN); return
end
try
    mdl = fitlm(table(x,y),'y ~ x - 1');
catch ME
    fit.reason = sprintf('fitlm failed: %s', ME.message); return
end
xf = linspace(min(lims),max(lims),200)';
[yf,ci] = predict(mdl,xf);
fill(ax,[xf;flipud(xf)],[ci(:,1);flipud(ci(:,2))],col(2,:), ...
    'LineStyle','none','FaceAlpha',0.4,'HandleVisibility','off');
plot(ax,xf,yf,'-','Color',col(1,:),'LineWidth',1.5,'HandleVisibility','off');
fit.ok = true;
fit.slope   = mdl.Coefficients.Estimate(1);
fit.slopeCI = mdl.coefCI;
fit.r2      = mdl.Rsquared.Ordinary;
fit.mdl     = mdl;
end
