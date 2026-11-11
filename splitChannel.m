%  split tif containg multiple channels
clearvars;close all;clc;

dataPath = uigetdir(".","Data path...");

if ~isfolder(dataPath)
    error('Data path not found')
end

try
    animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');
catch
    animal = char(inputdlg('Enter animal ID:','Animal Input',[1 35],{'AA0000'}));
end

tifFiles = dir(fullfile(dataPath,[animal '*.tif']));
%%
[img,tmpHeader] = readSCIMtif(fullfile(tifFiles(1).folder,...
    tifFiles(1).name));

if  isstruct(img)
    splitTifs = cell(length(tifFiles),2);
    if ~(exist(fullfile(tifFiles(1).folder,'ch1'),'dir')==7)
        mkdir(fullfile(tifFiles(1).folder,'ch1'))
    end
    if ~(exist(fullfile(tifFiles(1).folder,'ch2'),'dir')==7)
        mkdir(fullfile(tifFiles(1).folder,'ch2'))
    end
    for tifM = 1:length(tifFiles)
        splitTifs(tifM,:) = splitTifChans(fullfile(tifFiles(tifM).folder,tifFiles(tifM).name));
        movefile(splitTifs{tifM,1},fullfile(tifFiles(tifM).folder,'ch1',tifFiles(tifM).name),'f');
        copyfile(fullfile(tifFiles(tifM).folder,strrep(tifFiles(tifM).name,'.tif','_PulseParams.mat')),fullfile(tifFiles(tifM).folder,'ch1'),'f');
        copyfile(fullfile(tifFiles(tifM).folder,strrep(tifFiles(tifM).name,'.tif','_Pulses.mat')),fullfile(tifFiles(tifM).folder,'ch1'),'f');
        movefile(splitTifs{tifM,2},fullfile(tifFiles(tifM).folder,'ch2',tifFiles(tifM).name),'f');
        copyfile(fullfile(tifFiles(tifM).folder,strrep(tifFiles(tifM).name,'.tif','_PulseParams.mat')),fullfile(tifFiles(tifM).folder,'ch2'),'f');
        copyfile(fullfile(tifFiles(tifM).folder,strrep(tifFiles(tifM).name,'.tif','_Pulses.mat')),fullfile(tifFiles(tifM).folder,'ch2'),'f');
    end
end