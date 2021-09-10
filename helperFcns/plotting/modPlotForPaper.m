function modPlotForPaper(varargin)
p = inputParser;
addOptional(p,'ratioBool',false,@islogical)
parse(p,varargin{:})
ratioBool = p.Results.ratioBool;

if ratioBool
    hYline = yline(1);
    hYline.Color = 'k';
    hYline.LineStyle = '--';
    hYline.LineWidth = 2;
    hYline.HandleVisibility = 'off';
end

ax = gca;
ax.FontSize = 16;
ax.LineWidth = 1.5;
ax.XColor = 'k';
ax.YColor = 'k';
ax.TickDir = 'out';
set(ax,'color','none')
box off