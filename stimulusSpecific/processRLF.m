% processRLF: response-level function across cells pooled from
% multiple animals.
%
% Loads each animal's anmlROIbyStim table, vertically concatenates them,
% computes per-cell RLFs and dB thresholds (cells identified by
% animal+roiID so collisions across animals are kept distinct), and plots
% the cohort-mean RLF across included cells.

addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'helperFcns')));

%% params
matFiles = {
    '/media/DATA/Ophys/Jinbo/AA0067_anmlROI_BPNstimTable.mat'
    % '/media/DATA/Ophys/Jinbo/AA####_anmlROI_BPNstimTable.mat'
    % '/media/DATA/Ophys/Jinbo/AA####_anmlROI_BPNstimTable.mat'
    };
nConsec = 3;      % min consecutive sig==1 dB levels for inclusion
dBlist  = [];     % [] -> use sort(unique(T.dBampl))

%% load & concatenate
tables = cell(numel(matFiles),1);
for i = 1:numel(matFiles)
    S = load(matFiles{i},'anmlROIbyStim');
    tables{i} = S.anmlROIbyStim;
end
T = vertcat(tables{:});

% sanity: report cells per animal
[~,perAnml] = findgroups(T(:,'animal'));
nCellPerAnml = splitapply(@(a,r) numel(unique(strcat(a,'_',r))),...
    T.animal, T.roiID, findgroups(T.animal));
fprintf('Loaded %d animals, %d total rows\n', numel(matFiles), height(T));
for i = 1:height(perAnml)
    fprintf('  %s: %d cells\n', perAnml.animal(i), nCellPerAnml(i));
end

%% compute
rlf = tableRLF(T,'nConsec',nConsec,'dBlist',dBlist);

fprintf('Included %d of %d cells (>=%d consecutive sig dB levels)\n',...
    rlf.nIncluded, rlf.nTotal, nConsec);
disp(rlf.cellInfo);

%% plot
figure('Color','w');
plotRLF(rlf,'showCells',true);
