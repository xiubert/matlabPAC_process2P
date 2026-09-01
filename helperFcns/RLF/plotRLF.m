function [hAx,hLine,hErr] = plotRLF(rlf,varargin)
% PLOTRLF  Line plot of response-level function (mean +/- SEM across cells).
%
%   [hAx,hLine,hErr] = plotRLF(rlf)
%   [...] = plotRLF(rlf,'ax',ax,'showCells',true,'color',[0 0 0])
%
%   Inputs:
%       rlf - struct returned by tableRLF.
%
%   Name/Value:
%       'ax'        - axes handle. Default: gca.
%       'showCells' - overlay individual included-cell points at each dB
%                     as gray dots (no connecting lines). Default false.
%       'color'     - line and marker color. Default [0.4 0.4 0.4].
%
%   Outputs:
%       hAx   - axes handle
%       hLine - the mean line. An errorbar object when SEM is drawn, since
%               errorbar renders line, markers and caps as one object.
%       hErr  - the same errorbar object, or empty when SEM is all NaN.
%
%   See also tableRLF, cellRLF.

p = inputParser;
addRequired(p,'rlf',@isstruct);
addParameter(p,'ax',[],@(x) isempty(x) || isgraphics(x,'axes'));
addParameter(p,'showCells',false,@islogical);
addParameter(p,'color',[0.4 0.4 0.4],@(x) isnumeric(x) && numel(x)==3);
parse(p,rlf,varargin{:});
rlf       = p.Results.rlf;
ax        = p.Results.ax;
showCells = p.Results.showCells;
col       = p.Results.color;

if isempty(ax); ax = gca; end

hold(ax,'on');

% Per-cell points first, so the mean line sits on top of them.
if showCells && ~isempty(rlf.RLFincl)
    [nC,nDB] = size(rlf.RLFincl);
    xJit = repmat(rlf.dBlist,nC,1) + (rand(nC,nDB)-0.5)*2;
    plot(ax,xJit(:),rlf.RLFincl(:),'o',...
        'MarkerSize',4,'MarkerEdgeColor','none',...
        'MarkerFaceColor',[0.7 0.7 0.7],'HandleVisibility','off');
end

% Mean as a LINE with markers. Error bars only when the SEM is meaningful:
% tableRLF returns all-NaN semRLF below 2 included cells, so a single-cell RLF
% draws a bare line rather than zero-height caps. The condition is written on
% semRLF (not the newer showBand field) so RLF structs saved before that field
% existed still render.
if any(~isnan(rlf.semRLF))
    % errorbar draws the connecting line, the markers and the caps as ONE
    % object, so hLine and hErr are the same handle in this branch.
    hLine = errorbar(ax,rlf.dBlist,rlf.meanRLF,rlf.semRLF,'-o',...
        'Color',col,'MarkerFaceColor',col,'MarkerEdgeColor',col,...
        'MarkerSize',6,'LineWidth',1.6,'CapSize',6);
    hErr = hLine;
else
    hLine = plot(ax,rlf.dBlist,rlf.meanRLF,'-o',...
        'Color',col,'MarkerFaceColor',col,'MarkerEdgeColor',col,...
        'MarkerSize',6,'LineWidth',1.6);
    hErr = gobjects(0);
end

xlabel(ax,'dB SPL');
ylabel(ax,'peak \DeltaF/F');

ttl = sprintf('RLF (n = %d / %d cells',rlf.nIncluded,rlf.nTotal);
if isfield(rlf,'nAnimalsIncl') && ~isnan(rlf.nAnimalsIncl)
    if rlf.nAnimalsIncl == 1; word = 'mouse'; else; word = 'mice'; end
    ttl = sprintf('%s, %d %s',ttl,rlf.nAnimalsIncl,word);
end
title(ax,[ttl ')']);

if rlf.nIncluded == 0
    text(ax,0.5,0.5,sprintf('no cells passed\n>= %d consecutive sig. levels',rlf.nConsec), ...
        'Units','normalized','HorizontalAlignment','center', ...
        'VerticalAlignment','middle','Color',[0.6 0.2 0.1]);
elseif rlf.nIncluded == 1
    text(ax,0.02,0.98,'n = 1 cell - no SEM','Units','normalized', ...
        'HorizontalAlignment','left','VerticalAlignment','top', ...
        'FontSize',8,'Color',[0.6 0.2 0.1]);
end

set(ax,'XTick',rlf.dBlist);
box(ax,'off');
hAx = ax;
end
