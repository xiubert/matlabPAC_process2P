function figSaveAsFigEpsPng(hFig, varargin)
% FIGSAVEASFIGEPSPNG  Save a figure as .pdf, .png, and .fig.
%
%   figSaveAsFigEpsPng(hFig)
%   figSaveAsFigEpsPng(hFig, figLoc)
%
%   Saves hFig using hFig.Name as the base filename.
%   hFig.Name must be set before calling this function.
%
%   Inputs:
%     hFig   - figure handle; hFig.Name is used as the output filename base
%     figLoc - (optional) directory to save into (default: current directory)
%
%   See also plotMeanImgROI, exportgraphics, saveas

if isempty(varargin)
    figLoc = '.';
else
    figLoc = varargin{1};
end

set(hFig, 'Renderer', 'painters')
set(findall(hFig, '-property', 'FontName'), 'FontName', 'Arial')

% Must replace cmr font in Illustrator with Arial. EPS is best for
% illustrator import.
% IGNORE: disable TeX interpreter — use Unicode characters (Δ μ etc.) in labels
            % instead of TeX sequences to avoid CMR font embedding
            % set(findall(hFig, '-property', 'Interpreter'), 'Interpreter', 'none')

set(hFig, 'Color', 'none')

% exportgraphics(hFig, fullfile(figLoc, [hFig.Name '.pdf']), ...
%     'ContentType', 'vector', 'BackgroundColor', 'none')
exportgraphics(hFig, fullfile(figLoc, [hFig.Name '.png']), ...
    'BackgroundColor', 'white', 'Resolution', 300)

saveas(hFig, fullfile(figLoc, [hFig.Name '.eps']), 'epsc')
saveas(hFig, fullfile(figLoc, [hFig.Name '.pdf']), 'pdf')
% saveas(hFig, fullfile(figLoc, [hFig.Name '.png']))
saveas(hFig, fullfile(figLoc, [hFig.Name '.fig']))
