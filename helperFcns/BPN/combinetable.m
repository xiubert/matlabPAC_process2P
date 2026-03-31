%%combine data from multiple animals or multiple regions
filter = {'*.mat','Data files'};
fileList = {};   % collected full paths

while true
    [files, folder] = uigetfile(filter, 'Select file(s) (Cancel to finish)', pwd, 'MultiSelect', 'on');
    if isequal(files,0)
        % User cancelled -> stop selection loop
        break
    end

    % Normalize to cell array
    if ischar(files)
        files = {files};
    end

    % Build full paths and append to fileList
    fullPaths = fullfile(folder, files);
    if ischar(fullPaths)  % fullfile may return char for single file
        fullPaths = {fullPaths};
    end
    fileList = [fileList; fullPaths(:)]; %#ok<AGROW>

end

choice = questdlg('multiple animals or multiple regions in single animal?', '?', 'multiple animals', 'multiple regions', 'multiple animals');
n=numel(fileList);
combinedTable = table();
for k = 1:n
    load(fileList{k}, 'GroupedTbl');  % Load the GroupedTbl from the selected file
    if strcmp(choice, 'multiple regions')
        GroupedTbl.roiID=GroupedTbl.roiID+100*(k-1);
    end
    combinedTable = [combinedTable; GroupedTbl]; %#ok<AGROW> % Append to the combined table 
end

GroupedTbl=combinedTable;

%% Save
savePath='';
saveFileName = 'combinedData.mat';  % Define the name for the saved file
save(fullfile(savePath,saveFileName),'GroupedTbl','-v7.3');