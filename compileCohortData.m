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

params.cohort = cohortName;
params.colors.lohiPre = [0,0.451000000000000,0.741200000000000;0.851000000000000,0.329400000000000,0.102000000000000];
params.colors.lohiPost = [0,0.302000000000000,0.490200000000000;0.588200000000000,0.231400000000000,0.078400000000000];
params.colors.lohiTracePre = [0.729400000000000,0.874500000000000,1;1,0.694100000000000,0.541200000000000];
params.colors.lohiTracePost = [0.529400000000000,0.674500000000000,0.800000000000000;0.800000000000000,0.494100000000000,0.341200000000000];
params.colors.ratio = [0.651000000000000,0.651000000000000,0.651000000000000;0.149000000000000,0.149000000000000,0.149000000000000];

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
clearvars -except Tinput params outputTableSaveDir cohortName

%% save output table

if ~isfolder(outputTableSaveDir)
    mkdir(outputTableSaveDir)
end
save(fullfile(outputTableSaveDir,[cohortName '_dataTable.mat']),'Tinput','-v7.3')
save(fullfile(outputTableSaveDir,[cohortName '_params.mat']),'params')
