function [hAx,hBar,hErr] = plotRLF(rlf,varargin)
% PLOTRLF  Bar plot of response-level function (mean +/- SEM across cells).
%
%   [hAx,hBar,hErr] = plotRLF(rlf)
%   [...] = plotRLF(rlf,'ax',ax,'showCells',true,'color',[0 0 0])
%
%   Inputs:
%       rlf - struct returned by tableRLF.
%
%   Name/Value:
%       'ax'        - axes handle. Default: gca.
%       'showCells' - overlay individual included-cell points at each dB
%                     as gray dots (no connecting lines). Default false.
%       'color'     - bar face color. Default [0.4 0.4 0.4].
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

hBar = bar(ax,rlf.dBlist,rlf.meanRLF,...
    'FaceColor',col,'EdgeColor','none','BarWidth',0.7);
hold(ax,'on');
hErr = errorbar(ax,rlf.dBlist,rlf.meanRLF,rlf.semRLF,...
    'LineStyle','none','Color','k','LineWidth',1.2,'CapSize',8);

if showCells && ~isempty(rlf.RLFincl)
    [nC,nDB] = size(rlf.RLFincl);
    xJit = repmat(rlf.dBlist,nC,1) + (rand(nC,nDB)-0.5)*2;
    plot(ax,xJit(:),rlf.RLFincl(:),'o',...
        'MarkerSize',4,'MarkerEdgeColor','none',...
        'MarkerFaceColor',[0.6 0.6 0.6],'HandleVisibility','off');
end

xlabel(ax,'dB SPL');
ylabel(ax,'peak \DeltaF/F');
title(ax,sprintf('RLF (n = %d / %d cells)',rlf.nIncluded,rlf.nTotal));
set(ax,'XTick',rlf.dBlist);
box(ax,'off');
hAx = ax;
end
