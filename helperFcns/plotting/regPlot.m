function [mdl, ciRegLine, hScatter, hRegLineCI, hRegLine, hRef, hFig] = regPlot(X,Y,varargin)
% REGPLOT  Scatter plot with a linear regression line and 95% CI band.
%
%   [mdl, ciRegLine, hScatter, hRegLineCI, hRegLine, hRef, hFig] = regPlot(X, Y)
%   [...] = regPlot(X, Y, Name, Value, ...)
%
%   Fits a linear model (Y ~ X) using fitlm, plots the scatter, regression
%   line, and shaded 95% confidence interval band.  An optional reference
%   line (e.g. unity line) can be overlaid.
%
%   Required Inputs:
%     X, Y          - numeric vectors of equal length
%
%   Name-Value Inputs:
%     'regLine'      - logical; draw regression line and CI band (default: true)
%     'intercept'    - logical; include intercept in model (default: true)
%     'colors'       - 2x3 RGB matrix; row 1 = line color, row 2 = CI fill
%                      color (default: black line, gray fill)
%     'logScale'     - logical; log10-transform X and Y before fitting
%                      (default: false)
%     'sigID'        - grouping vector for gscatter; when length > 1, points
%                      are split into two groups styled as 'o' and '.'
%                      (default: false)
%     'fName'        - figure name string (default: 'scatter plot with
%                      regression line and CI')
%     'refLine'      - logical; overlay a reference line (default: false)
%     'mRefLine'     - slope of reference line (default: 1)
%     'bRefLine'     - intercept of reference line (default: 0)
%     'squareAxis'   - logical; set equal x/y limits spanning both axes
%                      (default: false)
%
%   Outputs:
%     mdl        - LinearModel object from fitlm
%     ciRegLine  - Nx2 matrix of 95% CI bounds from predict
%     hScatter   - handle to scatter plot
%     hRegLineCI - handle to CI fill patch
%     hRegLine   - handle to regression line
%     hRef       - handle to reference line (empty if 'refLine' is false)
%     hFig       - handle to figure

p = inputParser;
addRequired(p, 'X', @isnumeric);
addRequired(p, 'Y', @isnumeric);
addParameter(p,'regLine',true,@islogical)
addParameter(p,'intercept',true,@islogical)
addParameter(p,'colors',[],@(x) isequal(size(x),[2 3]))
addParameter(p,'logScale',false,@islogical)
addParameter(p,'sigID',false,@islogical)
addParameter(p,'fName','scatter plot with regression line and CI',@ischar)
addParameter(p,'refLine',false,@islogical)
addParameter(p,'mRefLine',1,@isnumeric)
addParameter(p,'bRefLine',0,@isnumeric)
addParameter(p,'squareAxis',false,@islogical)

parse(p, X, Y, varargin{:})
X = p.Results.X;
Y = p.Results.Y;
intercept = p.Results.intercept;
colors = p.Results.colors;
logScale = p.Results.logScale;
sigID = p.Results.sigID;
fName = p.Results.fName;
refLine = p.Results.refLine;
mRefLine = p.Results.mRefLine;
bRefLine = p.Results.bRefLine;
regLine = p.Results.regLine;
squareAxis = p.Results.squareAxis;

%scatter with regression line and confidence interval
if logScale
    X = log10(X);
    Y = log10(Y);
end
tbl = table(X,Y);
if ~intercept
    mdl = fitlm(tbl,'Y ~ X-1'); % '-1' means remove intercept
else
    mdl = fitlm(tbl,'Y ~ X'); % '-1' means remove intercept
end

g = groot;
if isempty(g.Children) || ~strcmp(fName,'scatter plot with regression line and CI')
    hFig = figure('Name',fName);
end

if length(sigID)>1
    hScatter = gscatter(X,Y,sigID,'kk','o.','HandleVisibility','off');
else
    hScatter = scatter(X,Y,'k','o','HandleVisibility','off');
end
if squareAxis
    xlim(round([min([xlim ylim]) max([xlim ylim])]))
    ylim(round([min([xlim ylim]) max([xlim ylim])]))
end
hold on
if refLine
    hRef = refline(mRefLine,bRefLine,'HandleVisibility','off');
    hRef.LineStyle = '--';
    % hRef.Color = colors.ratio(1,:);
    hRef.Color = 'r';
    hRef.LineWidth = 2;
end
if min([xlim ylim])>0
    if squareAxis
        xRegLine = linspace(min([xlim ylim]),max([xlim ylim]),1000);
    else
        xRegLine = linspace(min(xlim),max(xlim),1000);
    end
else
    if squareAxis
        xRegLine = linspace(0,max([xlim ylim]),1000);
    else
        xRegLine = linspace(0,max(xlim),1000);
    end
end
% yRegLine = xRegLine.*mdl.Coefficients.Estimate;
[yRegLine,ciRegLine] = predict(mdl,xRegLine(:));

if regLine
    hRegLine = plot(xRegLine,yRegLine,...
        'LineWidth',1.5);
    if ~isempty(colors)
        hRegLine.Color = colors(1,:);
        hRegLineCI = fill([xRegLine';flipud(xRegLine')],...
            [ciRegLine(:,1);flipud(ciRegLine(:,2))],...
            colors(2,:),'linestyle','none',...
            'HandleVisibility','off');
    else
        hRegLine.Color = 'k';
        hRegLineCI = fill([xRegLine';flipud(xRegLine')],...
            [ciRegLine(:,1);flipud(ciRegLine(:,2))],...
            [0.6510    0.6510    0.6510],'linestyle','none',...
            'HandleVisibility','off');
    end
    
    chi=get(gca, 'Children');
    set(gca, 'Children',flipud(chi));
end

title('')
axis square