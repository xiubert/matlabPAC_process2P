function figSaveAsFigEpsPng(hFig, varargin)
% FIGSAVEASFIGEPSPNG  Save a figure as .eps, .png, and .fig.
%
%   figSaveAsFigEpsPng(hFig)
%   figSaveAsFigEpsPng(hFig, figLoc)
%
%   Saves hFig in three formats using hFig.Name as the base filename.
%   hFig.Name must be set before calling this function.
%
%   Inputs:
%     hFig   - figure handle; hFig.Name is used as the output filename base
%     figLoc - (optional) directory to save into (default: current directory)
%
%   See also plotMeanImgROI, saveas

if isempty(varargin)
    figLoc = '.';
else
    figLoc = varargin{1};
end

saveas(hFig,fullfile(figLoc,[hFig.Name '.eps']),'epsc')
% print(hFig,'-depsc','-r300',fullfile(figLoc,[hFig.Name '.eps'])) %same
% (hFig,fullfile(figLoc,[hFig.Name '.eps']),'epsc')

saveas(hFig,fullfile(figLoc,[hFig.Name '.png']))
saveas(hFig,fullfile(figLoc,[hFig.Name '.fig']))