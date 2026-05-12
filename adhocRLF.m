% adhocRLF: build response-level functions from an anmlROIbyStim table.
%
% Loads a stim table, computes per-cell RLFs and dB thresholds, and plots
% the cohort-mean RLF across cells that have >= nConsec consecutive
% significant dB levels.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'helperFcns')));

%% params
matFile  = '/media/DATA/Ophys/Jinbo/AA0067_anmlROI_BPNstimTable.mat';
nConsec  = 3;      % min consecutive sig==1 dB levels for inclusion
dBlist   = [];     % [] -> use sort(unique(T.dBampl))

%% load
S = load(matFile,'anmlROIbyStim');
T = S.anmlROIbyStim;

%% compute
rlf = tableRLF(T,'nConsec',nConsec,'dBlist',dBlist);

fprintf('Included %d of %d cells (>=%d consecutive sig dB levels)\n',...
    rlf.nIncluded, rlf.nTotal, nConsec);
disp(rlf.cellInfo);

%% plot
figure('Color','w');
plotRLF(rlf,'showCells',true);
