clearvars;close all;clc;
%%
dataPath = uigetdir(".","Data path...");
animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');
tifFiles = dir(fullfile(dataPath,[animal '*.tif']));
%%
%load roi
[inputROIfile,roiFilepath] = uigetfile({'*_roiOutput*.mat;*moCorrROI*.mat'},...
            'Select output file containing ROI locations', dataPath);
inputROI = load([roiFilepath inputROIfile]);
if isfield(inputROI,'fluo2p')
    loadedROI = inputROI.fluo2p.roi;
    roiSvar = 'fluo2p';
elseif isfield(inputROI,'moCorROI')
    loadedROI = inputROI.moCorROI;
    roiSvar = 'moCorROI';
end


%%
sz = [length(tifFiles), 3]; 
varTypes = ["string", "double", "double"]; 
varNames = ["Filename", "Ch1", "Ch2"];

result = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);
result(:,1) = {tifFiles.name}';

for i = 1:length(tifFiles)

    [img,~] = readSCIMtif(fullfile(tifFiles(i).folder,tifFiles(i).name));
    imgCh1 = img(1).img;
    imgCh2 = img(2).img;
    rawFroiCh1 = roiMask2rawFroi(imgCh1,loadedROI);
    rawFroiCh2 = roiMask2rawFroi(imgCh2,loadedROI);
    result(i,2) = {mean(mean(rawFroiCh1))};
    result(i,3) = {mean(mean(rawFroiCh2))};
end
%%
saveFileName = 'darknoise.mat';  % Define the name for the saved file
save(fullfile(dataPath,saveFileName),'result','-v7.3');
T= result;
%%
% Assumptions: Your original table is named 'T'
% T = table(Filename, Ch1, Ch2, PMTA, PMTB, wavelength);

% 1. Calculate the average Ch2 for each unique combination of PMTB and wavelength
% 'groupsummary' automatically names the resulting average column as 'mean_Ch2'
avgTable = groupsummary(T2filtered, {'PMTB', 'wavelength'}, 'mean', 'Ratio');

% 2. Find all unique wavelengths to determine the number of subplots
uniqueWavelengths = unique(avgTable.wavelength);
numWavelengths = length(uniqueWavelengths);

figure('Name', 'Average Ch1/Ch2 vs PMTB by Wavelength');

% 4. Dynamically determine subplot grid layout (e.g., trying to keep it square-ish)
numRows = ceil(sqrt(numWavelengths));
numCols = ceil(numWavelengths / numRows);

% 5. Loop through each wavelength and create its subplot
for i = 1:numWavelengths
    currentWavelength = uniqueWavelengths(i);
    
    % Filter the averaged table for the current wavelength
    subData = avgTable(avgTable.wavelength == currentWavelength, :);
    
    % Sort by PMTB to ensure the line plot connects dots sequentially
    subData = sortrows(subData, 'PMTB');
    subplot(numRows, numCols, i);
    plot(subData.PMTB, subData.mean_Ratio, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');

    title(['Wavelength: ' num2str(currentWavelength) ' nm']);
    xlabel('PMTB');
    ylabel('Average Ch1/Ch2');
    grid on;
end
%% same but for Ch1

avgTable = groupsummary(T2filtered, {'PMTA', 'wavelength'}, 'mean', 'Ratio');

% 2. Find all unique wavelengths to determine the number of subplots
uniqueWavelengths = unique(avgTable.wavelength);
numWavelengths = length(uniqueWavelengths);

% 3. Create a figure for the subplots
figure('Name', 'Average Ch1/Ch2 vs PMTA by Wavelength');

% 4. Dynamically determine subplot grid layout (e.g., trying to keep it square-ish)
numRows = ceil(sqrt(numWavelengths));
numCols = ceil(numWavelengths / numRows);

% 5. Loop through each wavelength and create its subplot
for i = 1:numWavelengths
    currentWavelength = uniqueWavelengths(i);
    
    % Filter the averaged table for the current wavelength
    subData = avgTable(avgTable.wavelength == currentWavelength, :);
    
    % Sort by PMTB to ensure the line plot connects dots sequentially
    subData = sortrows(subData, 'PMTA');
    
    % Create the subplot
    subplot(numRows, numCols, i);
    
    % Plot Average Ch2 (mean_Ch2) vs PMTB
    plot(subData.PMTA, subData.mean_Ratio, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    
    % Formatting the subplot
    title(['Wavelength: ' num2str(currentWavelength) ' nm']);
    xlabel('PMTA');
    ylabel('Average Ch1/Ch2');
    grid on;
end
%%
% 1. Calculate the average Ch2 for each unique combination of PMTB and wavelength
% 'groupsummary' automatically names the resulting average column as 'mean_Ch2'
avgTable = groupsummary(T, {'PMTB', 'wavelength'}, 'mean', 'Ch2');

% 2. Find all unique PMTB values to determine the number of subplots
uniquePMTB = unique(avgTable.PMTB);
numPMTB = length(uniquePMTB);

% 3. Create a figure for the subplots
figure('Name', 'Average Ch2 vs Wavelength by PMTB');

% 4. Dynamically determine subplot grid layout
numRows = ceil(sqrt(numPMTB));
numCols = ceil(numPMTB / numRows);

% 5. Loop through each PMTB and create its subplot
for i = 1:numPMTB
    currentPMTB = uniquePMTB(i);
    
    % Filter the averaged table for the current PMTB
    subData = avgTable(avgTable.PMTB == currentPMTB, :);
    
    % Sort by wavelength to ensure the line connects sequentially left-to-right
    subData = sortrows(subData, 'wavelength');
    
    % Create the subplot
    subplot(numRows, numCols, i);
    
    % Plot Average Ch2 (mean_Ch2) vs Wavelength
    plot(subData.wavelength, subData.mean_Ch2, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
    
    % Formatting the subplot
    title(['PMTB: ' num2str(currentPMTB)]);
    xlabel('Wavelength (nm)');
    ylabel('Average Ch2');
    grid on;
end

%%
% Assumptions: Your original table is named 'T'
% T = table(Filename, Ch1, Ch2, PMTA, PMTB, wavelength);

% 1. Calculate the average Ch2 for each unique combination of PMTB and wavelength
avgTable = groupsummary(T, {'PMTB', 'wavelength'}, 'mean', 'Ch2');

% 2. Find all unique PMTB values
uniquePMTB = unique(avgTable.PMTB);

% 3. Create a single figure for the overlaid lines
figure('Name', 'Ch2 Ratio vs Wavelength (Normalized to 920 nm)');
hold on; % Keep previous lines when drawing new ones

% 4. Loop through each PMTB to calculate the ratio and plot
for i = 1:length(uniquePMTB)
    currentPMTB = uniquePMTB(i);
    
    % Filter and sort data for the current PMTB by wavelength
    subData = avgTable(avgTable.PMTB == currentPMTB, :);
    subData = sortrows(subData, 'wavelength');
    
    % Find the mean_Ch2 value specifically at wavelength == 920
    baselineIdx = (subData.wavelength == 920);
    
    % Check if the 920 nm baseline exists for this PMTB
    if any(baselineIdx)
        baselineValue = subData.mean_Ch2(baselineIdx);
        
        % Calculate the ratio for all wavelengths in this PMTB group
        subData.Ratio = subData.mean_Ch2 / baselineValue;
        
        % Plot Ratio vs Wavelength as a line, adding a DisplayName for the legend
        plot(subData.wavelength, subData.Ratio, '-o', 'LineWidth', 1.5, ...
             'DisplayName', ['PMTB = ' num2str(currentPMTB)]);
    else
        warning(['PMTB ' num2str(currentPMTB) ' does not have a data point at wavelength 920 nm. Skipping.']);
    end
end

% 5. Finalize plot formatting
hold off;
xlabel('Laser Wavelength (nm)');
ylabel('F(Normalized to 920 nm)');
%title('Ch2 Ratio Relative to 920 nm by PMTB');
grid on;
legend('Location', 'best'); % Automatically creates legend using DisplayNames