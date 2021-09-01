function figSaveAsFigEpsPng(hFig,varargin)

if isempty(varargin)
    figLoc = 'D:\Data\tmpFigs';
else
    figLoc = varargin{1};
end

saveas(hFig,fullfile(figLoc,[hFig.Name '.eps']),'epsc')
% print(hFig,'-depsc','-r300',fullfile(figLoc,[hFig.Name '.eps'])) %same
% (hFig,fullfile(figLoc,[hFig.Name '.eps']),'epsc')

saveas(hFig,fullfile(figLoc,[hFig.Name '.png']))
saveas(hFig,fullfile(figLoc,[hFig.Name '.fig']))