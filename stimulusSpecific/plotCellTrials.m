function out = plotCellTrials(src,animal,roiID,varargin)
% plotCellTrials  Every repetition of ONE cell, one panel per stimulus level.
%
%   out = plotCellTrials(groupFile,animal,roiID)
%   out = plotCellTrials(T,animal,roiID,'family','CGC','showRawF',false)
%
%   Follow-up tool for outliers spotted with plotCGCgroup's / plotBPNgroup's
%   showCells view. Those panels draw each cell's TRIAL-AVERAGED trace, so a
%   cell that towers over the others could be one bad repetition or a
%   genuinely large response on most of them. This draws the repetitions
%   themselves, so the two cases are distinguishable.
%
%   Family-agnostic: the level column, the per-trial trace column and the
%   time axis all come from stimGroupSpec, so it works for BPN sound levels
%   and CGC contrasts alike.
%
%   Inputs:
%     src    - group .mat path, or a group table.
%     animal - animal ID, e.g. "TO0007".
%     roiID  - ROI ID, e.g. "7" (string or numeric; compared as a string).
%
%   Name/Value:
%     'family'   - 'BPN' | 'CGC'. Inferred from the filename when src is a path.
%     'showRawF' - add a second row with the raw fluorescence per repetition,
%                  which distinguishes a real response from a motion or
%                  bleaching artifact. Default true.
%     'colors'   - nLevels x 3 line colours. Default: contrast colours for
%                  CGC, jet for anything else.
%     'verbose'  - print the per-repetition peaks. Default true.
%
%   Output (struct):
%     .animal .roiID .family .levels .labels
%     .trials  per level: .M (nReps x nFrames) .t .mean .n .perRepMax
%     .rawF    per level: .M (nReps x nFrames) .t   (empty if showRawF false)
%     .fig     figure handle
%
%   .perRepMax is the quick answer to "is it just one trace": a single large
%   entry among small ones is one bad repetition; several large ones are not.
%
%   See also plotCGCgroup, plotBPNgroup, stimGroupSpec, loadStimGroup.

p = inputParser;
addParameter(p,'family','',@(x) ischar(x)||isstring(x));
addParameter(p,'showRawF',true,@islogical);
addParameter(p,'colors',[],@(x) isempty(x)||isnumeric(x));
addParameter(p,'verbose',true,@islogical);
parse(p,varargin{:});
family   = char(p.Results.family);
showRawF = p.Results.showRawF;
colors   = p.Results.colors;
verbose  = p.Results.verbose;

%% ---- load ----
if istable(src)
    T = src;
    if isempty(family)
        error('plotCellTrials:noFamily','Pass ''family'' when src is a table.');
    end
else
    rep = validateStimGroup(src,'verbose',false);
    if isempty(family); family = rep.family; end
    T = loadStimGroup(src,'family',family);
end
spec = stimGroupSpec(family);

animal = string(animal);
roiID  = string(roiID);
rows = find(string(T.animal)==animal & string(T.roiID)==roiID);
if isempty(rows)
    error('plotCellTrials:cellNotFound', ...
        'No rows for animal %s roiID %s.', animal, roiID);
end

levels = sort(unique(T.(spec.levelVar)(rows)))';
nL     = numel(levels);
if isempty(colors)
    if strcmp(family,'CGC') && nL == 2
        c = getContrastColors(); colors = c.lohiPre;
    else
        colors = jet(nL);
    end
end

out.animal = animal; out.roiID = roiID; out.family = family;
out.levels = levels;
out.labels = strings(1,nL);
out.trials = struct('level',{},'label',{},'M',{},'t',{},'mean',{},'n',{},'perRepMax',{});
out.rawF   = struct('level',{},'M',{},'t',{});

%% ---- figure ----
nRow = 1 + showRawF;
f = figure('Color','w','Position',[50 50 460*nL+140 420*nRow], ...
    'Name',sprintf('%s_%s_roi%s_trials',family,animal,roiID));
tl = tiledlayout(f,nRow,nL,'Padding','compact','TileSpacing','compact');
title(tl,sprintf('%s  %s  ROI %s  -  every repetition',family,animal,roiID), ...
    'FontWeight','bold','FontSize',13);
out.fig = f;

for k = 1:nL
    r = rows(T.(spec.levelVar)(rows) == levels(k));
    r = r(1);
    M = T.(spec.cellTrialVar){r};
    t = T.(spec.timeVar){r};
    lab = levelLabel(family,levels,k);
    out.labels(k) = lab;
    out.trials(k) = struct('level',levels(k),'label',lab,'M',M,'t',t, ...
        'mean',mean(M,1,'omitnan'),'n',size(M,1), ...
        'perRepMax',max(M,[],2)');

    ax = nexttile(tl,k); hold(ax,'on');
    % every repetition, then the cell mean bold on top
    plot(ax,t,M','-','Color',[0.6 0.6 0.6 0.75],'LineWidth',0.75);
    plot(ax,t,mean(M,1,'omitnan'),'-','Color',colors(k,:),'LineWidth',2.4);
    stimMarker(ax,T,r,family);
    xlabel(ax,'time (s)'); ylabel(ax,'\DeltaF/F');
    title(ax,sprintf('%s  (%d reps)',lab,size(M,1)));
    box(ax,'off'); hold(ax,'off');
end

if showRawF
    for k = 1:nL
        r = rows(T.(spec.levelVar)(rows) == levels(k)); r = r(1);
        F  = T.SCALEDfissaFroi{r};
        tF = T.t_total{r};
        out.rawF(k) = struct('level',levels(k),'M',F,'t',tF);
        ax = nexttile(tl,nL+k); hold(ax,'on');
        plot(ax,tF,F','-','Color',[0.6 0.6 0.6 0.75],'LineWidth',0.75);
        plot(ax,tF,mean(F,1,'omitnan'),'-','Color',colors(k,:),'LineWidth',2.4);
        xlabel(ax,'time (s)'); ylabel(ax,'raw F (FISSA-scaled)');
        title(ax,sprintf('%s  raw F',out.labels(k)));
        box(ax,'off'); hold(ax,'off');
    end
end

if verbose
    fprintf('%s %s ROI %s\n',family,animal,roiID);
    for k = 1:nL
        fprintf('  %-14s %2d reps | per-rep max: %s\n', out.trials(k).label, ...
            out.trials(k).n, mat2str(round(out.trials(k).perRepMax,3)));
    end
end
end

%% ---- helpers ----
function lab = levelLabel(family,levels,k)
if strcmp(family,'CGC')
    if numel(levels)==2
        opts = ["Low contrast" "High contrast"]; lab = opts(k);
    else
        lab = sprintf('%g dB range',levels(k));
    end
else
    lab = sprintf('%g dB',levels(k));
end
end

function stimMarker(ax,T,r,family)
switch family
    case 'CGC'
        if ismember('PTsOnset',T.Properties.VariableNames)
            xline(ax,T.PTsOnset(r),'--','pure tone');
        end
    case 'BPN'
        if ismember('BPNsOnset',T.Properties.VariableNames)
            xline(ax,T.BPNsOnset(r),'--','BPN');
            if ismember('BPNmsPulseLen',T.Properties.VariableNames)
                xline(ax,T.BPNsOnset(r)+T.BPNmsPulseLen(r)/1000,'--');
            end
        end
end
end
