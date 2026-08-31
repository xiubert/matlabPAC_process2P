function s = cohortStat(x,y,varargin)
% cohortStat  Guarded two-sample comparison that refuses rather than returns NaN.
%
%   s = cohortStat(x,y)
%   s = cohortStat(x,y,'paired',true,'minN',3,'test','auto')
%
%   Wraps ttest / ttest2 / signrank / ranksum so that a comparison with too
%   few cells reports WHY it was skipped instead of emitting a NaN p-value
%   that reaches a figure. ttest with n <= 1 returns NaN without erroring,
%   which is exactly how "p = NaN" ends up printed on a plot.
%
%   Inputs:
%     x, y - numeric vectors. For paired tests they must be the same length
%            and are matched element-wise (pairs with a NaN in either are
%            dropped). For unpaired they may differ in length.
%
%   Name/Value:
%     'paired' - paired comparison. Default true.
%     'minN'   - minimum usable n to run the test at all. Default 3.
%     'test'   - 'auto' (default): Shapiro-free heuristic -- use the t-test
%                unless a Lilliefors test rejects normality, then fall back to
%                the rank-based equivalent. 'ttest' / 'rank' force a choice.
%     'alpha'  - significance level. Default 0.05.
%
%   Output (struct):
%     .ok      false when the test was not run
%     .reason  why it was skipped ('' when ok)
%     .test    name of the test actually used ('' when skipped)
%     .n       usable n (pairs for paired, [nx ny] for unpaired)
%     .p       p-value, NaN when skipped
%     .h       hypothesis result, NaN when skipped
%     .alpha   the alpha used
%     .stars   '***' / '**' / '*' / 'n.s.' / '' when skipped
%
%   Callers must check .ok before rendering .p. When .ok is false, print
%   .reason on the axes instead of a p-value.
%
%   See also groupN, cohortMeanSEM, sigDiffCalc.

p = inputParser;
addRequired(p,'x',@isnumeric);
addRequired(p,'y',@isnumeric);
addParameter(p,'paired',true,@islogical);
addParameter(p,'minN',3,@(v) isnumeric(v)&&isscalar(v)&&v>=1);
addParameter(p,'test','auto',@(v) any(strcmpi(v,{'auto','ttest','rank'})));
addParameter(p,'alpha',0.05,@(v) isnumeric(v)&&isscalar(v));
parse(p,x,y,varargin{:});
paired = p.Results.paired;
minN   = p.Results.minN;
which_ = lower(p.Results.test);
alpha  = p.Results.alpha;

s = struct('ok',false,'reason','','test','','n',0,'p',NaN,'h',NaN, ...
    'alpha',alpha,'stars','');

x = x(:); y = y(:);

if paired
    if numel(x) ~= numel(y)
        s.reason = sprintf('paired test needs equal lengths (got %d and %d)', numel(x), numel(y));
        return
    end
    keep = ~isnan(x) & ~isnan(y);
    x = x(keep); y = y(keep);
    s.n = numel(x);
    if s.n < minN
        s.reason = sprintf('n = %d complete pair%s, below minN = %d', ...
            s.n, plural(s.n), minN);
        return
    end
else
    x = x(~isnan(x)); y = y(~isnan(y));
    s.n = [numel(x) numel(y)];
    if any(s.n < minN)
        s.reason = sprintf('n = [%d %d], below minN = %d', s.n(1), s.n(2), minN);
        return
    end
end

% ---- pick the test ----
useRank = strcmp(which_,'rank');
if strcmp(which_,'auto')
    if paired; d = x - y; else; d = [x; y]; end
    useRank = ~isNormal(d,alpha);
end

try
    if paired && useRank
        [s.p,hh] = signrank(x,y,'alpha',alpha);  s.test = 'signrank';  s.h = double(hh);
    elseif paired
        [hh,s.p] = ttest(x,y,'Alpha',alpha);     s.test = 'ttest';     s.h = double(hh);
    elseif useRank
        [s.p,hh] = ranksum(x,y,'alpha',alpha);   s.test = 'ranksum';   s.h = double(hh);
    else
        [hh,s.p] = ttest2(x,y,'Alpha',alpha);    s.test = 'ttest2';    s.h = double(hh);
    end
catch ME
    s.reason = sprintf('%s failed: %s', s.test, ME.message);
    s.test = ''; s.p = NaN; s.h = NaN;
    return
end

if isnan(s.p)
    s.reason = sprintf('%s returned NaN (degenerate input)', s.test);
    s.test = ''; s.h = NaN;
    return
end

s.ok = true;
s.stars = starsFor(s.p);
end

% ---- helpers ----
function tf = isNormal(d,alpha)
% Lilliefors when available and usable; assume normal otherwise so the
% t-test remains the default rather than silently switching on an error.
tf = true;
if numel(d) < 4 || ~exist('lillietest','file'); return; end
try
    tf = ~lillietest(d,'Alpha',alpha);
catch
    tf = true;
end
end

function s = starsFor(p)
if     p < 0.001; s = '***';
elseif p < 0.01;  s = '**';
elseif p < 0.05;  s = '*';
else;             s = 'n.s.';
end
end

function s = plural(n)
if n == 1; s = ''; else; s = 's'; end
end
