%compile stimParamROI tables from animals in a cohort directory
function T = compileAnmlROItables(varargin)
defParams.tableDir = '.';

switch nargin
    case 0
        pDir = uigetdir('D:\Data\','Select cohort directory');
        params = defParams;
        params.parentPath = pDir;
    case 1
        if ischar(varargin{1}) && isfolder(varargin{1})
            warning('no params entered, using default')
            params = defParams;
            params.parentPath = varargin{1};
        elseif isstruct(varargin{1})
            params = varargin{1};
        end
end

aDir = dir(params.parentPath);
aDir = aDir(cell2mat({aDir.isdir}));
aDir = aDir(cellfun(@any,regexp({aDir.name}','[A-Z]{2}\d{4}')));
T = struct('anmlROIbyStim',[]);

for aNum = 1:length(aDir)
    disp(['Compiling stim/response data for: ' aDir(aNum).name]);
    if isfolder(fullfile(params.parentPath,aDir(aNum).name,params.tableDir))

        aData = load(fullfile(params.parentPath,aDir(aNum).name,params.tableDir,...
            [aDir(aNum).name ...
            '_anmlROI_stimTable.mat']),'anmlROIbyStim');
        aData.anmlROIbyStim.F = aData.anmlROIbyStim.SCALEDfissaFroi;
        aData.anmlROIbyStim = removevars(aData.anmlROIbyStim,...
            aData.anmlROIbyStim.Properties.VariableNames(...
            contains(aData.anmlROIbyStim.Properties.VariableNames,'Froi')));
        aData.anmlROIbyStim = removevars(aData.anmlROIbyStim,...
            {'Pulse','dBrange','PTampl','stimID'});
        T = cat(1,T,aData);
     
        clear aData  
       
    else
        disp([params.tableDir ' directory does not exist for ' ...
            aDir(aNum).name])
        disp([aDir(aNum).name ' not loaded into table'])
    end

end
clear aNum
T = T(2:end);
T = vertcat(T.anmlROIbyStim);
T.animal = string(T.animal);
T.roiID = str2num(T.roiID);