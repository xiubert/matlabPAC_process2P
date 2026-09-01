% PROCESSCGC  CGC / pure-tone-in-contrast (contrast gain control) analysis.
%
%   Loads <animal>_anmlROI_CGCstimTable_raw.mat (written by stimParam2ROI),
%   computes dF/F referenced to the DRC baseline and then to the pre-pure-tone
%   (PT) baseline, detects peak PT responses + significance, plots per-ROI and
%   population summaries, and writes the processed bundle to
%   <animal>_anmlROI_CGCstimTable.mat.
%
%   TWO-STAGE CONVENTION (matches BPN): the _raw input is never modified, so
%   re-running this script is safe and repeatable, and an animal that has not
%   been processed has no processed file for aggregateStimGroup to pick up
%   silently. A pre-split animal with only the un-suffixed file is still
%   handled: it is re-processed with its derived columns stripped first.
%
%   METHOD (matches matlabPAC_CGCplot/plotDataTable.m and the manuscript):
%     dFF_DRC          = (F - F0_DRC)/F0_DRC, F0_DRC = mean F over [-1.2 0] s (pre DRC onset).
%     dFF_PT  = dFF_DRC - F0_PT, F0_PT = mean of that dF/F over [1 2] s
%                        (an additive 2nd baseline subtraction, NOT a divisive
%                        re-normalization off raw F). Per-trial, as plotDataTable.
%
%   SIGNIFICANCE: per-ROI PT plot, and the population trace all use the
%   CELL-AVERAGE response (dFF_PT_avg = mean over reps), not individual
%   trials. pkFcalc is fed dFF_PT_avg.
%
%   SCOPE: this script processes and plots ONE animal. Its population panels
%   are delegated to plotCGCgroup, so the single-animal and cohort cases run
%   the same code. For cohort work across animals, aggregate with
%   aggregateStimGroup and call plotCGCgroup on the group file rather than
%   running this script.
%
%   PT FREQUENCY: an animal recorded with several pure-tone frequencies has
%   more than one row per (ROI, contrast); set PTfreqSelect to choose one.
%   With PTfreqSelect empty and several present, this raises
%   processCGC:multiplePTfreq rather than picking silently.
%
%   See also plotCGCgroup, aggregateStimGroup, dFoFcalc, pkFcalc,
%   getContrastColors, fillSEMplot

%% ---- PARAMETERS ----
% dF/F baseline windows (seconds, re trial start after trigDelay correction)
tBaseDRC   = [-1.2 0];   % F0 window before DRC onset
tBasePT    = [1 2];      % F0_PT window before pure tone (valid for PTsOnset==2)
PTonsetSec = 2;          % pure-tone onset (s); used for plot markers/xlines

% PT FREQUENCY SELECTION
% An animal may have been recorded with more than one pure-tone frequency
% crossed with contrast (TO0001 has 6484 and 30844 Hz), which gives several
% rows per (ROI, contrast). The per-ROI peak matrix below assumes ONE row per
% (ROI, contrast), and the historical group tables carry a single frequency per
% animal -- so a frequency has to be chosen. [] means "use the only frequency
% present" and raises a directive error if there is more than one. Set it here,
% or define PTfreqSelect in the workspace before running (which is how
% processAnimalStimFamilies passes it in).
if ~exist('PTfreqSelect','var'); PTfreqSelect = []; end

% peak-response detection
pkPTframeBin = 4;        % peak-search window length (frames) after PT onset
pkPTsigSD    = 2;        % significance threshold, in baseline-SD multiples

% plotting
ROIperFig       = 9;     % ROI subplots per figure (3x3)
colors          = getContrastColors();  % contrast color scheme (lohiPre/lohiTracePre)
avgTraceXlim    = [1 5]; % x-limits (s) for population average-trace plot
pkScatterLim    = [0 1]; % axis limits for low-vs-high peak dF/F scatter
popTraceSigOnly = true;  % population avg trace: true = significant cells only,
                         % false = all cells. Scatter/bar/t-test always use
                         % significant-in-both-contrasts cells.

%% ---- LOAD DATA ----
% Two-stage convention (matches BPN): stimParam2ROI writes
% <animal>_anmlROI_CGCstimTable_raw.mat; this script reads it and writes the
% processed bundle to <animal>_anmlROI_CGCstimTable.mat. Re-running never
% mutates its own input.
if ~exist('dataPath','var')
    dataPath = uigetdir('','Select the animal data folder');
    if isequal(dataPath,0)
        error('processCGC:noDataPath','No data folder selected.');
    end
end
if ~exist('animal','var') || isempty(animal)
    animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');
end

rawFile  = fullfile(dataPath,[animal '_anmlROI_CGCstimTable_raw.mat']);
procFile = fullfile(dataPath,[animal '_anmlROI_CGCstimTable.mat']);

if isfile(rawFile)
    load(rawFile)
elseif isfile(procFile)
    % Legacy layout: before the split, stimParam2ROI wrote the un-suffixed
    % file and processCGC -append'ed into it, so the only file present holds
    % raw AND derived columns. Re-process from it, but strip the derived
    % columns first -- otherwise stale or foreign ones survive into the fresh
    % run, which is how Groups C/D ended up carrying both dFF_PT_avg and a
    % dFF_PT_preDRCf0 written by a different script version.
    warning('processCGC:legacyLayout',...
        ['%s not found; falling back to %s (pre-split layout). Derived '...
         'columns are stripped and recomputed. Re-run stimParam2ROI to '...
         'create the _raw artifact.'],...
        [animal '_anmlROI_CGCstimTable_raw.mat'],...
        [animal '_anmlROI_CGCstimTable.mat']);
    load(procFile)
    stale = intersect(stimGroupSpec('CGC').derivedVars, ...
        anmlROIbyStim.Properties.VariableNames);
    if ~isempty(stale)
        fprintf('  stripping stale derived column(s): %s\n', strjoin(stale,', '));
        anmlROIbyStim = removevars(anmlROIbyStim, stale);
    end
else
    error('processCGC:noInput',...
        'Neither %s nor %s found in %s. Run stimParam2ROI first.',...
        [animal '_anmlROI_CGCstimTable_raw.mat'],...
        [animal '_anmlROI_CGCstimTable.mat'], dataPath);
end

%% PT frequency selection
if ismember('PTfreq',anmlROIbyStim.Properties.VariableNames)
    fq = unique(anmlROIbyStim.PTfreq);
    if isempty(PTfreqSelect)
        if numel(fq) > 1
            error('processCGC:multiplePTfreq', ...
                ['%s was recorded with %d pure-tone frequencies (%s Hz), so each ' ...
                 'cell has %d rows per contrast and the per-ROI peak matrix is ' ...
                 'ambiguous.\nSet PTfreqSelect (top of this script, or in the ' ...
                 'workspace) to the frequency to analyse. The group tables carry ' ...
                 'one frequency per animal.'], ...
                animal, numel(fq), strjoin(compose('%g',fq'),', '), numel(fq));
        end
    else
        keep = anmlROIbyStim.PTfreq == PTfreqSelect;
        if ~any(keep)
            error('processCGC:PTfreqNotPresent', ...
                'PTfreqSelect = %g Hz is not in %s. Available: %s Hz.', ...
                PTfreqSelect, animal, strjoin(compose('%g',fq'),', '));
        end
        fprintf('processCGC: %s -- selecting PTfreq %g Hz of [%s] (%d of %d rows)\n', ...
            animal, PTfreqSelect, strjoin(compose('%g',fq'),', '), sum(keep), numel(keep));
        anmlROIbyStim = anmlROIbyStim(keep,:);
        if exist('stimTable','var') && istable(stimTable) && ...
                ismember('PTfreq',stimTable.Properties.VariableNames)
            stimTable = stimTable(stimTable.PTfreq==PTfreqSelect,:);
        end
    end
end

%% Setup

% calculate time vector for each row based on frame rate and trigDelay
% 0 s corresponds to DRC stimulus start, time vector starts at -trigDelay (delay to Ephus stimulus trigger)
anmlROIbyStim.t_total = rowfun(@(F,fr,trigDelay) ...
    {(1:size(F,2))/fr-trigDelay},...
    anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi','frameRate','trigDelay'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

anmlROIbyStim.roiID = string(strtrim(cellstr(anmlROIbyStim.roiID)));
roiList = unique(anmlROIbyStim.roiID, 'stable');
nCell = numel(roiList);   % number of cells (unique ROIs); not stored in the .mat
remROIplotNo = rem(nCell,ROIperFig);
roiFigNo = floor(nCell/ROIperFig)+(remROIplotNo>=1);
dBdeltaList = unique(stimTable.dBdelta);
ndBdelta=length(dBdeltaList);

% Guard: the PT F0 window (tBasePT) must lie entirely before pure-tone onset,
% otherwise the "baseline" would include the tone response and every dFF_PT /
% dFF_PT / peak would be mis-normalized. tBasePT is hardcoded here
% (the source indexed it by round(PTsOnset)), so verify all rows are
% compatible rather than failing silently.
if tBasePT(2) > min(anmlROIbyStim.PTsOnset)
    error('processCGC:PTbaselineWindow',...
        ['tBasePT(2)=%.3g s is after the earliest PTsOnset (%.3g s); the PT ' ...
         'baseline would overlap the tone. Set tBasePT to a window before ' ...
         'onset (or branch on PTsOnset as in plotDataTable.m).'],...
        tBasePT(2), min(anmlROIbyStim.PTsOnset));
end
% (raw-F per-ROI diagnostic plot moved to the APPENDIX section at end of file)

%% dFF DRC
% dFF re pre DRC
% fluorescence traces (F(t) = ΔF/F ) first calculated as (F − F0_DRC)/ F0_DRC, 
% where F0_DRC is the average cell fluorescence intensity before DRC sound onset across −1.2 to 0 s. 
anmlROIbyStim.t_dFF_DRC = rowfun(@(t) ...
    {t(find((t>=tBaseDRC(1)),1,'first'):end)},...
    anmlROIbyStim,'InputVariables',{'t_total'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

anmlROIbyStim.dFF_DRC = rowfun(@(F,t) ...
    {dFoFcalc(F,[find((t>=tBaseDRC(1)),1,'first')...
    find((t<=tBaseDRC(2)),1,'last')],1)},...
    anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi','t_total'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

% cell-average dF/F re DRC across stim repetitions (1 x nFrames per row);
% used only for the per-ROI overview plot below. (The PT-response average that
% feeds significance is dFF_PT_avg, computed in the next section.)
anmlROIbyStim.dFF_DRC_avg = rowfun(@(F) ...
    {mean(F,1,'omitnan')},...
    anmlROIbyStim,'InputVariables',{'dFF_DRC'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

    %% PLOT dFF DRC
%initialize subplots for multiple ROI per fig
for roiFigN = 1:roiFigNo
    figure('Name','dF/F responses for each ROI');
    for roiSubPlotN = 1:ROIperFig
        curROIno = roiFigN*ROIperFig - ROIperFig + roiSubPlotN;
        if curROIno <= nCell
            subplot(3,3,roiSubPlotN);
            roi=roiList(curROIno);
            rows=anmlROIbyStim(anmlROIbyStim.roiID==roi,:);
            hold on;
            label=strings(height(rows),1);
            for r=1:height(rows)
                x=rows.t_dFF_DRC{r};
                y=rows.dFF_DRC_avg{r};
                plot(x,y);
                if rows.dBdelta(r) == dBdeltaList(1)
                    label(r)='Low contrast';
                elseif rows.dBdelta(r) == dBdeltaList(2)
                    label(r)='High contrast';
                end
            end
            xlabel('time/s')
            ylabel('dF/F')
            xline(PTonsetSec,'--')
            hold off;
            title('ROI'+string(curROIno));
            legend(label(1), label(2),'pure tone')
        end
        clear curROIno
    end
    clear roiSubPlotN
end

%
%% dFF re PT

% To quantify responses to pure tones that were preceded by 2 s of contrast DRC, 
% we calculated F(t)-F0_PT, where F0_PT is the average of F(t) across 1–2 s for a 2 s DRC duration.
% F(t) is already dF/F re pre-DRC, so this is a second baseline subtraction (re PT baseline) on top of the first (re DRC baseline). 

% PT response = F(t) - F0_PT, where F(t) is the existing pre-DRC dF/F (dFF_DRC)
% and F0_PT is the mean of that dF/F over the PT baseline window [1 2] s. This
% is an additive second baseline subtraction on the existing dF/F (NOT a fresh
% divisive (F-F0_PT)/F0_PT off raw F), and shares the t_dFF_DRC time axis.
% Per-trial (nReps x nFrames), exactly as plotDataTable.m.
anmlROIbyStim.dFF_PT = rowfun(@(dFF_DRC,t) ...
    {dFF_DRC  - ...
    nanmean(dFF_DRC(:,t>=tBasePT(1) & t<=tBasePT(2)),2)},...
    anmlROIbyStim,'InputVariables',{'dFF_DRC','t_dFF_DRC'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

% CELL-AVERAGE PT response across reps (1 x nFrames). This is the ONLY intended
% deviation from plotDataTable.m: significance (pkFcalc), the per-ROI PT plot,
% and the population trace are all computed from each cell's averaged response,
% NOT from individual trials. Averaging commutes with the linear F0_PT
% subtraction, so this equals the per-trial dFF_PT averaged over reps.
anmlROIbyStim.dFF_PT_avg = rowfun(@(dFF) ...
    {mean(dFF,1,'omitnan')},...
    anmlROIbyStim,'InputVariables',{'dFF_PT'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

%% PLOT dFF re PT re DRCf0
%initialize subplots for multiple ROI per fig
for roiFigN = 1:roiFigNo
    figure('Name','dF/F responses re PT re DRCf0 for each ROI');
    for roiSubPlotN = 1:ROIperFig
        curROIno = roiFigN*ROIperFig - ROIperFig + roiSubPlotN;
        if curROIno <= nCell
            subplot(3,3,roiSubPlotN);
            roi=roiList(curROIno);
            rows=anmlROIbyStim(anmlROIbyStim.roiID==roi,:);
            hold on;
            label=strings(height(rows),1);
            for r=1:height(rows)
                x=rows.t_dFF_DRC{r};
                y=rows.dFF_PT_avg{r};
                plot(x,y,'LineWidth', 2);
                if rows.dBdelta(r) == dBdeltaList(1)
                    label(r)='Low contrast';
                elseif rows.dBdelta(r) == dBdeltaList(2)
                    label(r)='High contrast';
                end
            end
            xlabel('time/s')
            ylabel('dF/F')
            xline(PTonsetSec,'--')
            hold off;
            title('ROI'+string(curROIno));
            legend(label(1), label(2),'pure tone')
        end
        clear curROIno
    end
    clear roiSubPlotN
end

%% Peak dFF PT response and significance
% pkFcalc receives the CELL-AVERAGE PT response (dFF_PT_avg, 1 x nFrames
% per row), so significance is thresholded on each cell's averaged trace, NOT on
% individual trials.
%
% SIGNIFICANCE BASELINE: pkFcalc uses frames 1:frameStart of the trace it is
% given as the baseline (mean + pkPTsigSD*SD). To make that baseline the 1-2 s
% pre-PT F0_PT window (sustained DRC just before the tone) we (1) crop the
% trace+time to t>=tBasePT(1), and (2) set frameStart at PT onset (NOT
% PTonset+1/fr). This makes baseline = [tBasePT(1), PT onset] exactly; using
% PTonset+1/fr would push the baseline one frame past onset, pulling the rising
% PT response into the baseline and inflating the SD. Including the DRC-onset
% transient (the whole pre-PT period) would likewise inflate the SD. The peak
% search runs nFrameWindow frames from PT onset; the onset frame sits in both
% baseline and peak window (the usual minor pkFcalc overlap) but precedes the
% Ca rise, so it does not affect the peak value.

tmp = rowfun(@(dFF,t,PTonset) ...
    pkFcalc(dFF(:, t>=tBasePT(1)),...
    find(t(t>=tBasePT(1))>=PTonset,1,'first'),...
    pkPTframeBin,pkPTsigSD),...
    anmlROIbyStim,'InputVariables',{'dFF_PT_avg','t_dFF_DRC','PTsOnset'},...
    'ExtractCellContents',true,'OutputFormat','cell','OutputVariableNames',{'sigPk','sig','pk'});
anmlROIbyStim.pkPT_sig = tmp(:,1);
anmlROIbyStim.sigPk = tmp(:,2);
anmlROIbyStim.pkPT = tmp(:,3);

%% Per-ROI peak + significance
% Per-ROI peak dF/F and significance flag: rows in roiList order, one column
% per contrast level (dBdeltaList order). Built by explicit roiID+dBdelta
% lookup so downstream indexing does NOT depend on the row order of
% anmlROIbyStim. pkByROI holds the peak of each cell-average trace (pkPT);
% sigByROI holds the pkFcalc significance flag for that same average trace.

pkByROI  = nan(nCell,ndBdelta);
sigByROI = false(nCell,ndBdelta);
for i = 1:nCell
    rows = anmlROIbyStim(anmlROIbyStim.roiID==roiList(i),:);
    for k = 1:ndBdelta
        sel = rows.dBdelta==dBdeltaList(k);
        v = cell2mat(rows.pkPT(sel));
        if ~isempty(v)
            pkByROI(i,k)  = v;                          % errors loudly if >1 row per (ROI,contrast)
            sigByROI(i,k) = logical(cell2mat(rows.sigPk(sel)));
        end
    end
end

% SIGNIFICANCE FILTER: keep only cells whose cell-AVERAGE peak is significant
% in BOTH contrasts (all(sigByROI,2)). This is a significance criterion from
% pkFcalc, NOT an amplitude (pkPT>0) threshold, and matches the source
% matlabPAC_CGCplot/plotDataTable.m (sigCellID = all(sig,2)). The scatter,
% bar graph and paired t-test below all use this mask.
valid = all(sigByROI,2);

%% Save output
% Full write of the processed bundle (NOT -append). The _raw input is never
% touched, so re-running this script is safe and repeatable. stimTable and
% tifStimParamTable are carried through because downstream code reads them
% from this file; with -append they were only present by inheritance from
% whatever the file already contained.
save(procFile,'anmlROIbyStim','stimTable','tifStimParamTable','dataPath','-v7.3');

%% Population summaries (traces, low-vs-high scatter, paired bar)
% Delegated to plotCGCgroup so the single-animal case and the group case run
% the SAME code. The previous inline versions crashed whenever no cell was
% significant in both contrasts -- splitapply on an empty selection errors
% with "Group numbers must be a vector of positive integers" -- which is not
% a multi-animal-only problem: TO0006 alone has 0 of 24 such cells, so
% processCGC could not run to completion on that animal at all.
%
% plotCGCgroup renders a labelled empty panel instead, and its paired test
% refuses with a reason rather than emitting p = NaN.
%
% For cohort work across animals, aggregate with aggregateStimGroup and call
% plotCGCgroup on the group file directly rather than running this script.
popOut = plotCGCgroup(anmlROIbyStim, ...
    'sigOnly',    popTraceSigOnly, ...
    'traceXlim',  avgTraceXlim, ...
    'scatterLim', pkScatterLim, ...
    'verbose',    true);


%% ========================== APPENDIX ==========================
% Optional QC plots, not part of the standard analysis. Opt-in via toggle so
% a normal script run does not spawn many figures.

%% APPENDIX: raw-F per-ROI diagnostic
% Plots raw (un-normalized) average fluorescence per ROI, ROIperFig per
% figure (3x3), one trace per contrast. Useful for spotting bleaching/motion
% artifacts that dF/F normalization would hide.
runRawFdiagnostic = false;   % set true to generate the raw-F QC figures

if runRawFdiagnostic
    % raw F averaged across stim repetitions (1 x nFrames per row)
    anmlROIbyStim.F_avg = rowfun(@(F) ...
        {mean(F,1,'omitnan')},...
        anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi'},...
        'ExtractCellContents',true,'OutputFormat','uniform');

    for roiFigN = 1:roiFigNo
        figure('Name','raw F responses for each ROI');
        for roiSubPlotN = 1:ROIperFig
            curROIno = roiFigN*ROIperFig - ROIperFig + roiSubPlotN;
            if curROIno <= nCell
                subplot(3,3,roiSubPlotN);
                roi=roiList(curROIno);
                rows=anmlROIbyStim(anmlROIbyStim.roiID==roi,:);
                hold on;
                label=strings(height(rows),1);
                for r=1:height(rows)
                    x=rows.t_total{r};
                    y=rows.F_avg{r};
                    plot(x,y);
                    if rows.dBdelta(r) == dBdeltaList(1)
                        label(r)='Low contrast';
                    elseif rows.dBdelta(r) == dBdeltaList(2)
                        label(r)='High contrast';
                    end
                end
                xlabel('time/s')
                ylabel('raw F')
                xline(PTonsetSec,'--','pure tone')
                hold off;
                title('ROI'+string(curROIno));
                legend(label)
            end
            clear curROIno
        end
        clear roiSubPlotN
    end
end
