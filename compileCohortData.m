%% COMPILE FRA and STIM DATA FROM COHORT
clearvars;close all;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  EDIT HERE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cohortName = 'CaMKII';
outputTableSaveDir = 'C:\Users\JIC402\OneDrive - University of Pittsburgh\Data';
params.parentPath = 'C:\Users\JIC402\OneDrive - University of Pittsburgh\Data\CaMKII_combined';

params.tableDir = '.';
params.treatment = 'none';
params.pkPTsigSD = 2;
params.nFramesPostPulse = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DONE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compile FRA and stim/response data
[TfraAnml,TfraROI] = compileAnmlFRA(params);
Tinput = compileAnmlROItables(params);

%isolate tone repsonsive cells
TfraROI = TfraROI(TfraROI.dPrime>0,:);

Tinput = innerjoin(Tinput,TfraROI,'Keys',{'animal','roiID'},...
    'RightVariables',{'BFuDB'});

%% clean up
clear Tfra*

% animalID and roiID 2 sequence
[G,iA] = findgroups(Tinput.animal);
tmp = splitapply(@(x) {findgroups(x)},Tinput.roiID,G);
tmp = vertcat(tmp{:});

Tinput.animal = G;
Tinput.roiID = tmp;
clearvars -except Tinput params

%% save output table

if ~isfolder(outputTableSaveDir)
    mkdir(outputTableSaveDir)
end
save(fullfile(outputTableSaveDir,[cohortName '_dataTable.mat']),'Tinput','-v7.3')
save(fullfile(outputTableSaveDir,[cohortName '_params.mat']),'params')
