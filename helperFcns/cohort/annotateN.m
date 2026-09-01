function h = annotateN(ax,N,varargin)
% annotateN  Stamp cell/animal counts (and small-n caveats) onto an axes.
%
%   h = annotateN(ax,N)
%   h = annotateN(ax,N,'location','northwest','extra','sig. in both contrasts')
%
%   Every cohort figure should state its n on the figure itself, not only in
%   a caption written later. This also stamps the caveats the degenerate-n
%   contract requires, so a reader cannot mistake a single-animal or
%   single-cell panel for a population result.
%
%   Inputs:
%     ax - axes handle.
%     N  - struct from groupN.
%
%   Name/Value:
%     'location' - 'northwest' (default), 'northeast', 'southwest', 'southeast'.
%     'extra'    - extra text appended to the n label, e.g. the inclusion
%                  criterion the count refers to.
%     'fontSize' - default 8.
%
%   Output:
%     h - handle to the text object.
%
%   Stamps added automatically:
%     nAnimals == 1  ->  "1 animal - no across-animal inference"
%     nCells   == 1  ->  "n = 1 cell - no SEM"
%     nCells   == 0  ->  "no cells passed criterion"
%
%   See also groupN, cohortMeanSEM, cohortStat.

p = inputParser;
addParameter(p,'location','northwest',@(x) any(strcmpi(x, ...
    {'northwest','northeast','southwest','southeast'})));
addParameter(p,'extra','',@(x) ischar(x)||isstring(x));
addParameter(p,'fontSize',8,@isnumeric);
parse(p,varargin{:});
loc   = lower(p.Results.location);
extra = char(p.Results.extra);
fs    = p.Results.fontSize;

txt = string(N.label);
if ~isempty(extra); txt = txt + " (" + string(extra) + ")"; end

if N.nCells == 0
    txt = txt + newline + "no cells passed criterion";
else
    if N.nCells == 1
        txt = txt + newline + "n = 1 cell - no SEM";
    end
    if N.singleAnimal
        txt = txt + newline + "1 animal - no across-animal inference";
    end
end

switch loc
    case 'northwest'; xy = [0.02 0.98]; va = 'top';    ha = 'left';
    case 'northeast'; xy = [0.98 0.98]; va = 'top';    ha = 'right';
    case 'southwest'; xy = [0.02 0.02]; va = 'bottom'; ha = 'left';
    case 'southeast'; xy = [0.98 0.02]; va = 'bottom'; ha = 'right';
end

h = text(ax, xy(1), xy(2), txt, ...
    'Units','normalized', 'VerticalAlignment',va, 'HorizontalAlignment',ha, ...
    'FontSize',fs, 'Interpreter','none', 'Margin',2);

if N.nCells <= 1 || N.singleAnimal
    set(h,'Color',[0.6 0.2 0.1]);   % flag the caveat rows
end
end
