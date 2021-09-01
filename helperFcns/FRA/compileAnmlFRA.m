function [TfraAnml,TfraROI] = compileAnmlFRA(varargin)
defParams.tableDir = '.';
defParams.treatment = 'pre';
defParams.pkPTsigSD = 2;
defParams.nFramesPostPulse = 3;

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

%get animal list
aDir = dir(params.parentPath);
aDir = aDir(cell2mat({aDir.isdir}));
aDir = aDir(cellfun(@any,regexp({aDir.name}','[A-Z]{2}\d{4}')));

%initialize tables and cells
TfraAnml = struct('FRAoutput',[]);
TfraROI = struct('animal',[],'roiID',[],'treatment',[],...
    'BFuDB',[],'dPrime',[]);
animals = cell(0);

for aNum = 1:length(aDir)
    disp(['Compiling FRA for: ' aDir(aNum).name]);
    
    if isfolder(fullfile(params.parentPath,aDir(aNum).name,params.tableDir))
        
        tmpDirs = dir(fullfile(params.parentPath,aDir(aNum).name,params.tableDir));
        tmp = load(fullfile(params.parentPath,aDir(aNum).name,params.tableDir,...
            [aDir(aNum).name ...
            '_anmlROI_stimTable.mat']),'tifFileList','fissaScaleFactor');
        tifFileList = tmp.tifFileList;
        clear tmp
        
        animals = cat(1,animals,regexp(aDir(aNum).name,'[A-Z]{2}\d{4}','match','once'));
        anmlTreatment = regexprep(tifFileList.stim(1).treatment,...
            '(pre|post)','');
        
        %get roiIDs
        if numel({tmpDirs(contains({tmpDirs.name},'_moCorrROI')).name})>1
            roiFileName = tmpDirs(contains({tmpDirs.name},...
                ['_moCorrROI_' strrep(params.treatment,anmlTreatment,'')])).name;
        else
            roiFileName = tmpDirs(contains({tmpDirs.name},'_moCorrROI')).name;
        end       
        sROI = load(fullfile(params.parentPath,aDir(aNum).name,params.tableDir,...
            roiFileName),'moCorROI');
        roiID = {sROI.moCorROI.ID}';
        clear sROI
        
        tifFileList.map = tifFileList.map(...
            contains({tifFileList.map.treatment},params.treatment));
        
        %get FRA for all cells for animal
        [sFRA.FRAoutput,~,~] = ...
            FRAmap(tifFileList,...
            params.pkPTsigSD,...
            params.nFramesPostPulse);
        
        %Load FRA for animal into table
        TfraAnml = cat(1,TfraAnml,sFRA);
        
        %Combine FRA for animal into structure ordered by ROI
        TfraROI = cat(1,TfraROI,cell2struct([repmat({aDir(aNum).name},...
            size(sFRA.FRAoutput.BFuDB))...
            roiID,...
            repmat({params.treatment},...
            size(sFRA.FRAoutput.BFuDB)),...
            num2cell(sFRA.FRAoutput.BFuDB) ...
            num2cell(sFRA.FRAoutput.dPrime)],...
            fieldnames(TfraROI)',2));
        clear sFRA roiID
             
    else
        disp([params.tableDir ' directory does not exist for ' ...
            aDir(aNum).name])
        disp([aDir(aNum).name ' not loaded into table'])
    end
end

TfraAnml = TfraAnml(2:end);
TfraAnml = vertcat(TfraAnml.FRAoutput);
[TfraAnml.animal] = animals{:};
TfraAnml = orderfields(TfraAnml, [12,1:11]);
TfraAnml = denseAmplTo3ampl(TfraAnml);

TfraROI = struct2table(TfraROI(2:end));
TfraROI.animal = string(TfraROI.animal);
TfraROI.treatment = string(TfraROI.treatment);
TfraROI.roiID = str2num(char(TfraROI.roiID));