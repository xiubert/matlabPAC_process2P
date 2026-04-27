function [hFig, img, roi] = plotMeanImgROI(varargin)
% PLOTMEANIMGROI  Plot mean projection of a tif series with ROI outlines.
%
%   plotMeanImgROI()
%   plotMeanImgROI(img, roi)
%   [hFig, img, roi] = plotMeanImgROI(...)
%
%   With no arguments, opens file browsers to select one or more tif files
%   and a saved ROI .mat file, concatenates the tifs, then plots.
%   Accepts both _moCorrROI_*.mat (moCorROI struct array) and
%   _roiOutput.mat (fluo2p.roi struct) formats.
%
%   Inputs (optional):
%     img  - H x W x nFrames numeric array
%     roi  - 1 x nROI struct array with fields: XYvertices, label, deleted
%              XYvertices is Nx2 (col1=X, col2=Y), matching roiGUI convention
%
%   Outputs (all optional):
%     hFig - figure handle; Name is set to the tif directory when called
%              with no args, otherwise set hFig.Name before calling
%              figSaveAsFigEpsPng
%     img  - H x W x nFrames array (useful when loaded via file browser)
%     roi  - ROI struct array (useful when loaded via file browser)
%
%   Non-deleted ROIs are outlined in green and labelled with their ID.

if nargin == 2
    img = varargin{1};
    roi = varargin{2};

elseif nargin == 0
    % select tif files
    [tifNames, tifDir] = uigetfile('*.tif', ...
        'Select tif file(s) to average', 'MultiSelect', 'on');
    if isequal(tifNames, 0)
        disp('No tif files selected.')
        return
    end
    if ischar(tifNames)
        tifNames = {tifNames};
    end

    % load and concatenate frames across tifs
    img = [];
    for k = 1:length(tifNames)
        frames = justLoadTif(fullfile(tifDir, tifNames{k}));
        img = cat(3, img, double(frames));
    end

    % select ROI file
    [roiFile, roiDir] = uigetfile('*.mat', 'Select ROI .mat file');
    if isequal(roiFile, 0)
        disp('No ROI file selected.')
        return
    end
    roi = loadROIfile(fullfile(roiDir, roiFile));

else
    error('plotMeanImgROI expects 0 or 2 arguments.')
end

% compute mean projection
meanImg = mean(img, 3, 'omitnan');

% derive figure name from source
if nargin == 0
    % use the tif directory as the figure name (animal/session ID)
    [~, figName] = fileparts(tifDir(1:end-1));
else
    figName = 'Mean projection with ROIs';
end

% plot
hFig = figure('Name', figName);
imagesc(meanImg)
colormap(gray)
axis image off
hold on

% overlay non-deleted ROI outlines and labels
% uses XYvertices (Nx2: col1=X, col2=Y) to match roiGUI plot_btn_callbk
for k = 1:length(roi)
    if roi(k).deleted
        continue
    end
    if isfield(roi(k), 'XYvertices') && ~isempty(roi(k).XYvertices)
        plot(roi(k).XYvertices(:,1), roi(k).XYvertices(:,2), 'g-', 'LineWidth', 1)
    elseif isfield(roi(k), 'ROIcurveOrderedXY') && ~isempty(roi(k).ROIcurveOrderedXY)
        % fallback: ROIcurveOrderedXY is 2xN (row1=X, row2=Y)
        plot(roi(k).ROIcurveOrderedXY(1,:), roi(k).ROIcurveOrderedXY(2,:), 'g-', 'LineWidth', 1)
    end
    text(roi(k).label.Position(1), roi(k).label.Position(2), ...
        roi(k).label.String, 'Color', 'g', 'FontWeight', 'bold', 'FontSize', 8)
end

end


function roi = loadROIfile(filepath)
% Handle _moCorrROI_ and _roiOutput.mat formats.
S = load(filepath);

if isfield(S, 'moCorROI')
    roi = S.moCorROI;

elseif isfield(S, 'fluo2p') && isfield(S.fluo2p, 'roi')
    roi = S.fluo2p.roi;

else
    error('Unrecognised ROI file format. Expected moCorROI or fluo2p.roi.')
end
end
