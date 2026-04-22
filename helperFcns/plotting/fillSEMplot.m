function [hPatch,hPlot] = fillSEMplot(x,y,SEM,colorY,colorSEM,ax)
% FILLSEMPLOT  Plot a line with a shaded SEM band.
%
%   [hPatch, hPlot] = fillSEMplot(x, y, SEM, colorY, colorSEM)
%   [hPatch, hPlot] = fillSEMplot(x, y, SEM, colorY, colorSEM, ax)
%
%   Inputs:
%     x        - 1 x N vector of x values
%     y        - 1 x N vector of y values (mean)
%     SEM      - 1 x N vector of SEM values; shaded band spans y ± SEM
%     colorY   - line color (any MATLAB color spec)
%     colorSEM - fill color for the SEM patch (any MATLAB color spec)
%     ax       - axes handle to plot into (default: gca)
%
%   Outputs:
%     hPatch - handle to the filled SEM patch (FaceAlpha 0.6, excluded from legend)
%     hPlot  - handle to the mean line (LineWidth 3)
if nargin==5
    ax = gca;
end

hPatch = fill(ax,[x';flipud(x')],...
    [y'-SEM';flipud(y'+SEM')],...
    colorSEM,'FaceAlpha','0.6',...
    'linestyle','none','HandleVisibility','off');
hold(ax,'on')
hPlot = plot(ax,x,y,'LineWidth',3,...
    'Color',colorY);