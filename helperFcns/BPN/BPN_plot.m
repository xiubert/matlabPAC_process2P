%% plot
% Inputs (set these)
targetROI = 2;      % desired roiID
targetdB  = 30;     % desired dB

% Find row(s) matching the ROI and dB
rowMask = (GroupedTbl.roiID == targetROI) & (GroupedTbl.dB == targetdB);
if ~any(rowMask)
    error('No rows found for roiID=%d and dB=%g', targetROI, targetdB)
end

r = find(rowMask,1);

% Extract the cell arrays of repeated measurements (n-by-1 cell)
timeCells = GroupedTbl.time{r};   % expected: n-by-1 cell, each cell m-by-1
dffCells  = GroupedTbl.dFF{r};    % same shape

% Convert to matrices: m-by-n
timeMat = horzcat(timeCells{:});   % each column is one repetition
dffMat  = horzcat(dffCells{:});


% use frame index to compute time
m = size(dffMat,1);
t = (1:m)'/fs;


% Compute mean across repetitions (along columns)
dffMean = mean(dffMat, 2, 'omitnan');

% Plot all repetitions and the average
figure;
hold on;
plot(t, dffMat, 'Color', [0.7 0.7 0.9]);        % thin light lines for individual repetitions
plot(t, dffMean, '-k', 'LineWidth', 2);         % thick black line for mean
xline(1,'--','BPN')
xline(1.4,'--')
hold off;
xlabel('Time (s)');
ylabel('dF/F');
title(sprintf('roiID=%d   dB=%g   %d repetitions', targetROI, targetdB, size(dffMat,2)));
legend({'Individual measurement','Average'});
grid on;
% filename=
% f = gcf; % Get the current figure handle
% exportgraphics(f,filename,'Resolution',300)

%% 
% Inputs
targetROI = 2;   % desired roiID

% Select rows for that ROI
rows = find(GroupedTbl.roiID == targetROI);
if isempty(rows)
    error('No rows found for roiID=%d', targetROI)
end

% Get unique dB values in those rows
dbVals = unique(GroupedTbl.dB(rows));

% Prepare figure
figure; hold on;
cmap = jet(numel(dbVals));
colors = cmap;

for k = 1:numel(dbVals)
    db = dbVals(k);
    % Find row(s) for this ROI and this dB
    mask = (GroupedTbl.roiID == targetROI) & (GroupedTbl.dB == db);
    rIdx = find(mask);
    if isempty(rIdx)
        continue
    end
    % If multiple rows for same dB, concatenate all repetitions from them
    allTimeCols = [];
    allDffCols  = [];
    for r = rIdx(:)'
        % Each of these is expected to be an n-by-1 cell, each cell m-by-1
        timeCells = GroupedTbl.time{r};
        dffCells  = GroupedTbl.dFF{r};
        % Convert to matrices: m-by-n
        if ~isempty(timeCells)
            allTimeCols = [allTimeCols, horzcat(timeCells{:})]; %#ok<AGROW>
            allDffCols  = [allDffCols,  horzcat(dffCells{:})];  %#ok<AGROW>
        end
    end
    if isempty(allDffCols)
        continue
    end

    m = size(allDffCols,1);
    t = (1:m)'/fs;

    % Compute mean and SEM across repetitions (columns)
    mu  = mean(allDffCols, 2, 'omitnan');
    sem = std(allDffCols, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(allDffCols),2));

    % Plot shaded SEM band
    upper = mu + sem;
    lower = mu - sem;
    x = [t; flipud(t)];
    y = [upper; flipud(lower)];
    hPatch = patch(x, y, colors(k,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');

    % Plot mean line
    plot(t, mu, 'Color', colors(k,:), 'LineWidth', 1.8);

    % optional: store handles for legend
    hLine(k) = plot(NaN, NaN, 'Color', colors(k,:), 'LineWidth', 1.8); %#ok<SAGROW>
end

% Finalize plot
xline(1,'--','BPN')
xline(1.4,'--')
xlabel('Time (s)');
ylabel('dF/F');
title(sprintf('ROI %d: average dF/F per dB', targetROI));
legend(hLine, arrayfun(@(v) sprintf('%d dB', v), dbVals, 'UniformOutput', false), 'Location', 'best');
grid on;
hold off;
%% peak dFF vs dB

% Input: GroupedTbl with variables: roiID, dB, MeansigPeak

% 1) Group and compute stats
[G, dBvals] = findgroups(GroupedTbl.dB);      % G groups rows by dB, dBvals are group keys
% data=cell2mat(GroupedTbl.sigPeak);

c = GroupedTbl.sigPeak;        % cell array, some cells empty
n = numel(c);
data = NaN(n,1);
nonEmpty = ~cellfun(@isempty, c);
if any(nonEmpty)
    data(nonEmpty) = cell2mat(c(nonEmpty));
end


mu = splitapply(@(x) mean(x,'omitnan'), data, G);     % mean per dB
n  = splitapply(@(x) sum(~isnan(x)), data, G);    % counts per dB
sd = splitapply(@(x) std(x,'omitnan'),  data, G);     % std per dB
sem = sd ./ sqrt(n);                                   % SEM per dB

% Sort dB values (optional) to ensure increasing x positions
[dBvals_sorted, idx] = sort(dBvals);
mu = mu(idx);
sem = sem(idx);

% 2) Bar plot with error bars
figure;
hb = bar(dBvals_sorted, mu, 'FaceColor',[0.7 0.7 0.9]);
hold on;
% errorbars: use errorbar with LineStyle 'none'
he = errorbar(dBvals_sorted, mu, sem, 'k', 'LineStyle', 'none', 'LineWidth',1);

% 3) Overlay individual data points (jittered horizontally)
% Extract MeansigPeak per group in sorted order
pointsByGroup = splitapply(@(x) {x}, data, G);
pointsByGroup = pointsByGroup(idx);   % reorder to match sorted dB

markerColor = [0 0.2 0.6];
for k = 1:numel(dBvals_sorted)
    x0 = dBvals_sorted(k);
    pts = pointsByGroup{k}(:);
    mcount = numel(pts);
    if mcount == 0
        continue
    end
    % jitter width: small fraction of spacing between dB values (or fixed)
    if numel(dBvals_sorted) > 1
        dx = 0.15 * min(diff(dBvals_sorted)); % scale with spacing
    else
        dx = 0.15; % fallback
    end
    jitter = (rand(mcount,1)*2-1) * dx;
    xs = x0 + jitter;
    % plot individual points
    scatter(xs, pts, 36, 'MarkerEdgeColor', markerColor, ...
        'MarkerFaceColor', markerColor, 'MarkerFaceAlpha', 0.6);
end
% % 3.5) Connect points that belong to the same roiID
% roiIDs = unique(GroupedTbl.roiID);
% for iR = 1:numel(roiIDs)
%     mask = (GroupedTbl.roiID == roiIDs(iR));
%     if nnz(mask) <= 1
%         continue
%     end
%     % use the same numeric dB values and the numeric data array 'data'
%     x_raw = GroupedTbl.dB(mask);
%     y_raw = data(mask);
% 
%     % map raw dB values to the sorted x positions (dBvals_sorted)
%     [found, loc] = ismember(x_raw, dBvals_sorted);
%     if ~all(found)
%         % if some dB values are missing from dBvals_sorted, skip those points
%         loc = loc(found);
%         y_raw = y_raw(found);
%     end
%     xplot = dBvals_sorted(loc);
% 
%     % sort by x (in case dB order differs)
%     [xplot, sidx] = sort(xplot);
%     yplot = y_raw(sidx);
% 
%     % draw thin semi-muted line connecting this ROI's points
%     plot(xplot, yplot, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8);
% end
% 4) Labels and aesthetics
xlabel('Sound level (dB)');
ylabel('peak dF/F');
% title('MeansigPeak by dB with SEM and individual data points');
box on;
hold off;