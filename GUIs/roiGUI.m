function [] = roiGUI(varargin)
% roiGUI  Load ScanImage .tif for roi drawing.
%   [] = readSCIMtif(varargin)
%   loads .tif into GUI for drawing and saving ROI
%
%   Additional input arguments: 
%       'inputFileName' --> loads specified .tif file
%       '1' --> loads specified channel, replace with desired channel number
%               channel id stored in imHeader
%       'merge' --> merges image data from 2 channels, assumes red,green
%
%   PAC_20200213
%
%   See also meanfluoROIvt.m


%PAC_20190522 overhaul for new 2P pipeline, removed epifluorescence
%also passing data outside guidata
%changed name to roiGUI

%%%%%%%%%%%%%%%%%%%%%%%%%% DEAL WITH INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(varargin) || ~any(cellfun(@ischar,varargin)) ...
        || ~any(cellfun(@(c) contains(c,'.tif'),varargin))
    [inputFileName, path] = uigetfile({'*.tif';'*roiOutput.mat'},'File Selector');
    if isnumeric(inputFileName) && inputFileName==0
        disp('No file chosen...')
        return
    end
    inputFileName = fullfile(path,inputFileName);
end

if any(cellfun(@ischar,varargin)) && any(cellfun(@(c) contains(c,'.tif'),varargin))
    inputFileName = varargin{cellfun(@(c) contains(c,'.tif'),varargin)};
    if exist(inputFileName,'file')~=2
        error('...File doesn''t exist...')
    end
end
    
if any(cellfun(@isnumeric,varargin)) && ...
        numel(varargin{cellfun(@isnumeric,varargin)})==1
    chanID = varargin{cellfun(@isnumeric,varargin)};
else
    %always load channel 2 unless otherwise specified
    chanID = 2;
end

if any(cellfun(@ischar,varargin)) && any(cellfun(@(c) contains(c,'merge'),varargin))
    mergeView = 1;
end

sROI.file = inputFileName;

%%%%%%%%%%%%%%%%%%%%%%%%%% LOOK FOR STIMULUS INFO %%%%%%%%%%%%%%%%%%%%%%%%%

%look for saved pulse and trigger info and if exsists put in base
%workspace
if exist(strrep(sROI.file,'.tif','_Pulses.mat'), 'file') == 2
    %trigger parameters are saved as [tif file name]_PulseParams.mat
    %via Scanimage user function 'digtrig_stimPulse_train'
    %ephus calls 'savePulseDetails' (located in ephus userFcns
    %directory) this will look for [tif file name]_PulseParams.mat and
    %included it in [tif file name]_Pulses.mat along with the pulse
    %info
    temp = load(strrep(sROI.file,'.tif','_Pulses.mat'));
    
    %show output of trigger parameters and pulse
    temp.params
    temp.pulse
    
    %store pulse and params in ROI gui data
    sROI.TriggerParams = temp.params;
    sROI.Stim = temp.pulse;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD TIF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%THIS SECTION IMPORTS THE TIF IMAGE DATA INTO A MATRIX

if exist('mergeView','var')==1
    [img,imHeader] = ...
        readSCIMtif(sROI.file);
    sROI.imgMultChan = img;
    img = img(chanID).img;
    sROI.rawFroiChan = chanID;
    
else
    [img,imHeader] = ...
        readSCIMtif(sROI.file,chanID);
end
  
sROI.nFrames = imHeader.nFrames;
sROI.imgWidth = imHeader.imWidth;
sROI.imgHeight = imHeader.imHeight;

%state is variable containing ScanImage state information, eg version,
%acquisition settings etc
sROI.state = imHeader;

%establish settings from header/state information
sROI.frameRate = round(imHeader.frameRate);
sROI.t = 0:1/sROI.frameRate:(sROI.nFrames/sROI.frameRate)-1/sROI.frameRate;
sROI.img = double(img);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   SETUP FIGURE  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize figure
%POSITION:  [X Y Width Height]

%Plot raw fluorescence traces with slider and frame ID
ui.roiGUI.fh = figure('Name',['Choose frame and draw ROI']);
if exist('mergeView','var')==1
    ui.roiGUI.plt = image(redGreenMerge(sROI.imgMultChan(1).img(:,:,1),...
        sROI.imgMultChan(2).img(:,:,1)));
else
    ui.roiGUI.plt = imagesc(sROI.img(:,:,1));
end
ui.roiGUI.ax = findall(ui.roiGUI.fh,'type','axes');

%initialize GUI settings
ui.roiGUI.curROItype = 'Ellipse'; %default ROI on start

%time and frame stamp in ROI gui
ui.roiGUI.frameNo = text(ui.roiGUI.ax,0.5,sROI.imgHeight-5,'Frame: 1','Color','g');
if isfield(sROI,'state')
    ui.roiGUI.timeStamp = text(ui.roiGUI.ax,0.5,5,['Time (s): ' num2str(sROI.t(1))],'Color','g');
end
colormap(gray)

%ROI gui features
ui.roiGUI.slider = uicontrol('Parent',ui.roiGUI.fh,'Style','slider',...
    'Units','Normalized',...
    'Position',[ui.roiGUI.ax.Position(1),...
    ui.roiGUI.ax.Position(2)-0.10,...
    ui.roiGUI.ax.Position(3),0.05],...
    'value',1, 'min',1, 'max',sROI.nFrames,...
    'SliderStep',[1/sROI.nFrames 1/sROI.nFrames],'callback',@slider_callback);

ui.roiGUI.exitBtn = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', 'Exit / Close All','Tag','exitbutton1',...
    'Units','Normalized',...
    'Position', [ui.roiGUI.ax.Position(1)+...
    ui.roiGUI.ax.Position(3)-0.045, ...
    ui.roiGUI.ax.Position(2),...
    0.14, 0.04],...
    'Callback',@exit_btn_callbk);

%resize image and figure
ui.roiGUI.ax.Position(3) = ui.roiGUI.ax.Position(3)-0.1;
ui.roiGUI.ax.Position(1) = ui.roiGUI.ax.Position(1)+0.05;
ui.roiGUI.ax.Position(4) = ui.roiGUI.ax.Position(4)-0.1;
ui.roiGUI.fh.Units = 'Normalized';
ui.roiGUI.fh.Position = [0.2099    0.0639    0.5521    0.8519];

%ROI Selection
ui.roiGUI.roiType = uicontrol('Parent',ui.roiGUI.fh,'Style', 'text',...
    'String','Choose ROI Type:','Units','Normalized',...
    'Position',[0.0302,0.2033,0.0915,0.0228]);

ui.roiGUI.bg = uibuttongroup('Parent',ui.roiGUI.fh,...
    'Visible','off',...
    'Position',[0.022641509433962,0.106521739130435,0.109433962264151,0.101086956521739],...
    'SelectionChangedFcn',@roiTypeSelection);

% Create two radio buttons in the button group.
ui.roiGUI.ellipse1 = uicontrol(ui.roiGUI.bg,'Style',...
    'radiobutton',...
    'String','Ellipse',...
    'Units','Normalized',...
    'Position',[0.205701791148595,0.580599555676005,0.870699999999998,0.1997]);
ui.roiGUI.freehand1 = uicontrol(ui.roiGUI.bg,'Style','radiobutton',...
    'String','Free Hand',...
    'Units','Normalized',...
    'Position',[0.211428571428571,0.231709432451245,0.8707,0.1997]);

%Channel Selection if 'Merge'
if isfield(sROI,'imgMultChan')
    ui.roiGUI.chanSel = uicontrol('Parent',ui.roiGUI.fh,'Style', 'text',...
        'String','Channel:','Units','Normalized',...
        'Position',[0.0481245,0.58,0.0915,0.0228]);
    
    ui.roiGUI.chanSelBG = uibuttongroup('Parent',ui.roiGUI.fh,...
        'Visible','on',...
        'Position',[0.036792452830188,0.491304347826087,0.109433962264151,0.0914],...
        'SelectionChangedFcn',@channelSelection);


    % Create two radio buttons in the button group.
    ui.roiGUI.merge = uicontrol(ui.roiGUI.chanSelBG,'Style',...
        'radiobutton',...
        'String','merge',...
        'Units','Normalized',...
        'Position',[0.1164,0.7,0.870699999999998,0.1997]);
    ui.roiGUI.chan1 = uicontrol(ui.roiGUI.chanSelBG,'Style','radiobutton',...
        'String','ch. 1 (red)',...
        'Units','Normalized',...
        'Position',[0.1164,0.45,0.870699999999998,0.1997]);
    ui.roiGUI.chan2 = uicontrol(ui.roiGUI.chanSelBG,'Style','radiobutton',...
        'String','ch. 2 (green)',...
        'Units','Normalized',...
        'Position',[0.1164,0.20,0.870699999999998,0.1997]);
end
%tif file in title
ui.roiGUI.title = uicontrol('Parent',ui.roiGUI.fh,'Style','text',...
    'String', sROI.file,...
    'Units','Normalized',...
    'Position', [ui.roiGUI.ax.Position(1), ...
    ui.roiGUI.ax.Position(2)+...
    ui.roiGUI.ax.Position(4)+0.12,...
    0.7, 0.04]);

%Show motor position information in the ROI gui if scanimage tif
if isfield(sROI,'state')
    ui.roiGUI.relMotorPosX = uicontrol('Parent',ui.roiGUI.fh,'Style','text',...
        'String', ['Relative X Pos: ' num2str(sROI.state.motorPos.relative(1))],...
        'Units','Normalized',...
        'Position', [ui.roiGUI.ax.Position(1)-0.1, ...
        ui.roiGUI.ax.Position(2)+...
        ui.roiGUI.ax.Position(4)+0.06,...
        0.4, 0.04]);
    ui.roiGUI.relMotorPosY = uicontrol('Parent',ui.roiGUI.fh,'Style','text',...
        'String', ['Relative Y Pos: ' num2str(sROI.state.motorPos.relative(2))],...
        'Units','Normalized',...
        'Position', [ui.roiGUI.ax.Position(1)-0.1, ...
        ui.roiGUI.ax.Position(2)+...
        ui.roiGUI.ax.Position(4)+0.03,...
        0.4, 0.04]);
    ui.roiGUI.relMotorPosZ = uicontrol('Parent',ui.roiGUI.fh,'Style','text',...
        'String', ['Relative Z Pos: ' num2str(sROI.state.motorPos.relative(3))],...
        'Units','Normalized',...
        'Position', [ui.roiGUI.ax.Position(1)-0.1, ...
        ui.roiGUI.ax.Position(2)+...
        ui.roiGUI.ax.Position(4),...
        0.4, 0.04]);
    ui.roiGUI.absMotorPosX = uicontrol('Parent',ui.roiGUI.fh,'Style','text',...
        'String', ['Absolute X Pos: ' num2str(sROI.state.motorPos.absolute(1))],...
        'Units','Normalized',...
        'Position', [ui.roiGUI.ax.Position(3)-0.05, ...
        ui.roiGUI.ax.Position(2)+...
        ui.roiGUI.ax.Position(4)+0.06,...
        0.25, 0.04]);
    ui.roiGUI.absMotorPosY = uicontrol('Parent',ui.roiGUI.fh,'Style','text',...
        'String', ['Absolute Y Pos: ' num2str(sROI.state.motorPos.absolute(2))],...
        'Units','Normalized',...
        'Position', [ui.roiGUI.ax.Position(3)-0.05, ...
        ui.roiGUI.ax.Position(2)+...
        ui.roiGUI.ax.Position(4)+0.03,...
        0.25, 0.04]);
    ui.roiGUI.absMotorPosZ = uicontrol('Parent',ui.roiGUI.fh,'Style','text',...
        'String', ['Absolute Z Pos: ' num2str(sROI.state.motorPos.absolute(3))],...
        'Units','Normalized',...
        'Position', [ui.roiGUI.ax.Position(3)-0.05, ...
        ui.roiGUI.ax.Position(2)+...
        ui.roiGUI.ax.Position(4),...
        0.25, 0.04]);  
end %if it's a ScanImage tif file

%show option buttons
ui.roiGUI.drawRoiBtn = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', 'Draw ROI','Tag','roibutton1',...
    'Units','Normalized',...
    'Position', [0.02811320754717,0.244565217391305,0.1,0.04],...
    'Callback',@roi_btn_callbk);
ui.roiGUI.loadRoiBtn = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', 'Load ROI Set','Tag','loadroibutton1',...
    'Units','Normalized',...
    'Position', [ui.roiGUI.ax.Position(1)-0.17, ...
    ui.roiGUI.ax.Position(4)+...
    ui.roiGUI.ax.Position(2)-0.04,...
    0.17, 0.04],...
    'Callback',@loadroi_btn_callbk);
ui.roiGUI.loadPrevRoiSet = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', 'Load Prev. ROI Set','Tag','loadPrevROISetButton1',...
    'Units','Normalized',...
    'Position', [ui.roiGUI.ax.Position(1)-0.17, ...
    ui.roiGUI.ax.Position(4)+...
    ui.roiGUI.ax.Position(2)-0.08,...
    0.17, 0.04],...
    'Callback',@loadPrevROI_btn_callbk);
ui.roiGUI.loadDelRoiBtn = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', 'Load Deleted ROI','Tag','roiHistButton1',...
    'Units','Normalized',...
    'Position', [ui.roiGUI.ax.Position(1)-0.17, ...
    ui.roiGUI.ax.Position(4)+...
    ui.roiGUI.ax.Position(2)-0.12,...
    0.17, 0.04],...
    'Callback',@loadDelROI_btn_callbk);
ui.roiGUI.delRoiBtn = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', 'Delete ROI(s)','Tag','delroibutton1',...
    'Units','Normalized',...
    'Position', [ui.roiGUI.ax.Position(1)+...
    ui.roiGUI.ax.Position(3)+0.005, ...
    ui.roiGUI.ax.Position(4)+...
    ui.roiGUI.ax.Position(2)-0.04,...
    0.14, 0.04],...
    'Callback',@delRoi_btn_callbk);
ui.roiGUI.clearRoiBtn = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', 'Clear ROIs','Tag','clrbutton1',...
    'Units','Normalized',...
    'Position', [ui.roiGUI.ax.Position(1)+...
    ui.roiGUI.ax.Position(3)+0.005, ...
    ui.roiGUI.ax.Position(4)+...
    ui.roiGUI.ax.Position(2)-0.08,...
    0.14, 0.04],...
    'Callback',@clr_btn_callbk);
ui.roiGUI.meanToggle = uicontrol('Parent',ui.roiGUI.fh,'Style','togglebutton',...
    'String','Show mean image','Tag','meanImgToggle','Units','Normalized',...
    'Position', [ui.roiGUI.ax.Position(1)-0.17, ...
    ui.roiGUI.ax.Position(4)+...
    ui.roiGUI.ax.Position(2)-0.2,...
    0.17, 0.04],...
    'Callback',@meanImgToggle_callbk);
ui.roiGUI.closeFigBtn = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', 'Close ROI Plots','Tag','closebutton1',...
    'Units','Normalized',...
    'Position', [ui.roiGUI.ax.Position(1)+...
    ui.roiGUI.ax.Position(3)+0.005, ...
    ui.roiGUI.ax.Position(4)+...
    ui.roiGUI.ax.Position(2)-0.16,...
    0.14, 0.04],...
    'Callback',@close_fig_callbk);
ui.roiGUI.plotBtn = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', 'Plot','Tag','donebutton1',...
    'Units','Normalized',...
    'Position', [ui.roiGUI.ax.Position(1)+...
    ui.roiGUI.ax.Position(3)+0.005, ...
    ui.roiGUI.ax.Position(4)+...
    ui.roiGUI.ax.Position(2)-0.32,...
    0.14, 0.04],...
    'Callback',@plot_btn_callbk);
ui.roiGUI.saveBtn = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', 'Save Output','Tag','savebutton1',...
    'Units','Normalized',...
    'Position', [ui.roiGUI.ax.Position(1)+...
    ui.roiGUI.ax.Position(3)+0.005, ...
    ui.roiGUI.ax.Position(4)+...
    ui.roiGUI.ax.Position(2)-0.40,...
    0.14, 0.04],...
    'Callback',@saveOutput_btn_callbk);
ui.roiGUI.outputBtn = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', 'Output to Workspace','Tag','savebutton1',...
    'Units','Normalized',...
    'Position', [0.86,0.465217391304348,0.14,0.04],...
    'Callback',@OutputWorkspace_btn_callbk);
ui.roiGUI.loadNextBtn = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', 'Load Next File','Tag','loadNextButton1',...
    'Units','Normalized',...
    'Position', [ui.roiGUI.ax.Position(1)+...
    ui.roiGUI.ax.Position(3)+0.005, ...
    ui.roiGUI.ax.Position(2)+0.05,...
    0.14, 0.04],...
    'Callback',@loadNextCallback);
ui.roiGUI.loadPrevBtn = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', 'Load Previous File','Tag','loadPrevButton1',...
    'Units','Normalized',...
    'Position', [ui.roiGUI.ax.Position(1)+...
    ui.roiGUI.ax.Position(3)+0.005, ...
    ui.roiGUI.ax.Position(2)+0.1,...
    0.14, 0.04],...
    'Callback',@loadPrevCallback);

%make figure visible
set(ui.roiGUI.fh,{'Name','Visible'},{['Choose frame and draw ROI'],'on'})
ui.roiGUI.bg.Visible = 'on';
%% CALLBACK FUNCTIONS for GUI

%closes everything without saving; leaves variables in workspace
    function [] = exit_btn_callbk(varargin)
        closeROIplots()
        close(ui.roiGUI.fh)
    end

%show mean of frames across time
    function [] = meanImgToggle_callbk(hObject,~)
        if get(hObject,'value')==1
            ui.roiGUI.slider.Visible = 0;
            if exist('mergeView','var')==1
                if ~isfield(ui.roiGUI,'curChannel') || contains(ui.roiGUI.curChannel,'merge')
                    ui.roiGUI.plt.CData = redGreenMerge(nanmean(sROI.imgMultChan(1).img,3),...
                        nanmean(sROI.imgMultChan(2).img,3));
                elseif strcmp(ui.roiGUI.curChannel,'ch. 1 (red)')
                    ui.roiGUI.plt.CData = nanmean(sROI.imgMultChan(1).img,3);
                elseif strcmp(ui.roiGUI.curChannel,'ch. 2 (green)')
                    ui.roiGUI.plt.CData = nanmean(sROI.imgMultChan(2).img,3);
                end
            else
                ui.roiGUI.plt.CData = nanmean(sROI.img,3);
            end
        elseif get(hObject,'value')==0
            ui.roiGUI.slider.Visible = 1;
            if ~isfield(ui.roiGUI,'frame')
                ui.roiGUI.frame = 1;
            end
            if exist('mergeView','var')==1
                
                if ~isfield(ui.roiGUI,'curChannel') || contains(ui.roiGUI.curChannel,'merge')
                    ui.roiGUI.plt.CData = redGreenMerge(sROI.imgMultChan(1).img(:,:,ui.roiGUI.frame),...
                        sROI.imgMultChan(2).img(:,:,ui.roiGUI.frame));
                elseif strcmp(ui.roiGUI.curChannel,'ch. 1 (red)')
                    ui.roiGUI.plt.CData = sROI.imgMultChan(1).img(:,:,ui.roiGUI.frame);
                elseif strcmp(ui.roiGUI.curChannel,'ch. 2 (green)')
                    ui.roiGUI.plt.CData = sROI.imgMultChan(2).img(:,:,ui.roiGUI.frame);
                end
                
            else
                ui.roiGUI.plt.CData = sROI.img(:,:,ui.roiGUI.frame);
            end
        end
    end %meanImgToggle_callbk

%slider for going through image frames
    function [] = slider_callback(hObject,~)
        %get frame number from slider position
        frame = round(get(hObject,'value'));
        
        %set frame number in GUI data so it is known which frame ROI was
        %drawn on
        ui.roiGUI.frame = frame;
        
        %update plot image data to respective frame number
        if exist('mergeView','var')==1
            if ~isfield(ui.roiGUI,'curChannel') || contains(ui.roiGUI.curChannel,'merge')
                ui.roiGUI.plt.CData = redGreenMerge(sROI.imgMultChan(1).img(:,:,frame),...
                    sROI.imgMultChan(2).img(:,:,frame));
                %could rewrite to regexp for channel and just use one line
            elseif strcmp(ui.roiGUI.curChannel,'ch. 1 (red)')
                ui.roiGUI.plt.CData = sROI.imgMultChan(1).img(:,:,frame);
            elseif strcmp(ui.roiGUI.curChannel,'ch. 2 (green)')
                ui.roiGUI.plt.CData = sROI.imgMultChan(2).img(:,:,frame);
            end
        else
            ui.roiGUI.plt.CData = sROI.img(:,:,frame);
        end
        
        %change frame number on plot
        ui.roiGUI.frameNo.String = ['Frame: ' num2str(frame)];
        
        %Show timestamp for tif files
        if isfield(sROI,'state')
            ui.roiGUI.timeStamp.String = ['Time (s): ' num2str(sROI.t(frame))];
        end
    end %slider_callback


%function to load ROIs from _roiOutput.mat file
%keeps deleted ROIs, but they are still deleted unless load deleted ROI
%button is chosen
    function [] = loadroi_btn_callbk(~,~)
        [inputROIfile,roiFilepath] = uigetfile({'*_roiOutput*.mat;*moCorrROI*.mat'},...
            'Select output file containing ROI locations');
        
        try
            %load fluo2p structure from _roiOutput.mat file
            inputROI = load([roiFilepath inputROIfile]);
            if isfield(inputROI,'fluo2p')
                loadedROI = inputROI.fluo2p.roi;
                roiSvar = 'fluo2p';
            elseif isfield(inputROI,'moCorROI')
                loadedROI = inputROI.moCorROI;
                roiSvar = 'moCorROI';
            end
            
            
            %for each ROI in roiOutput.mat file
            for roiN = 1:length(loadedROI)
                
                %get draw shape function from ROI shape type
                f = str2func(['@' loadedROI(roiN).type]);
                
                %call shape function on current ROI gui
                hROI = f(ui.roiGUI.ax,loadedROI(roiN).pos);
                
                %get ROI info
                ID = loadedROI(roiN).ID;
                pos = loadedROI(roiN).pos;
                frame = loadedROI(roiN).frame;
                deleted = loadedROI(roiN).deleted;
                
                %save ROI info to sROI structure
                sROI.roi(roiN) = struct('object',hROI,...
                    'ID',ID,'pos',pos,...
                    'label',text(ui.roiGUI.ax,...
                    'Position',loadedROI(roiN).label.Position,...
                    'String',loadedROI(roiN).label.String,'Color','g','FontWeight','bold'),...
                    'frame',frame,'deleted',deleted);
                
                %if the ROI was deleted don't show it on the ROI figure
                if loadedROI(roiN).deleted == 1
                    sROI.roi(roiN).object.delete;
                end
            end
            clear roiN
            
            %add position callback functions to each ROI so that they can
            %be repositioned and those positions will be saved
            for roiN = 1:length(sROI.roi)
                if sROI.roi(roiN).deleted == 0
                    addNewPositionCallback(sROI.roi(roiN).object,...
                        @(pos) getpos(pos,sROI.roi(roiN).object,roiN));
                end
            end
            clear roiN
            if isfield(inputROI.(roiSvar),'baselineIDX')
                ui.roiGUI.baselineIDX = sROI.t(inputROI.(roiSvar).baselineIDX);
            end
            if ~isfield(sROI,'frame')
                ui.roiGUI.frame = 1;
            end
        catch
            disp('No File Selected...')
        end
    end %load_ROI_callback

%looks for _roiOutput.mat file that is one trace number lower than the
%current trace number and loads ROIs from that file
    function [] = loadPrevROI_btn_callbk(~,~)
        
        %determine previous file number
        fileNo = str2double(string(regexprep(regexp(sROI.file,'[A-Z]{4}\d{4}','match'),'[A-Z]{4}','')));
        prevFileNo = fileNo-1;
        prevfileNoStr = sprintf('%04d', prevFileNo);
        
        tIDX = regexp(sROI.file,'[A-Z]{4}\d{4}');
        prevFileName = [sROI.file(1:tIDX+3) prevfileNoStr '.tif'];
        prevFileROI = strrep(prevFileName,'.tif','_roiOutput.mat');
        
        disp(['Previous ROI Set: ' prevFileROI])
        
        try
            inputROI = load(prevFileROI);
            
            %loads all ROI regardless of whether they've been deleted
            for roiN = 1:length(inputROI.fluo2p.roi)
                f = str2func(['@' inputROI.fluo2p.roi(roiN).type]);
                hROI = f(ui.roiGUI.ax,inputROI.fluo2p.roi(roiN).pos);
                ID = inputROI.fluo2p.roi(roiN).ID;
                pos = inputROI.fluo2p.roi(roiN).pos;
                try
                    ROIvertices = inputROI.fluo2p.roi(roiN).XYvertices;
                catch
                    ROIvertices = hROI.getVertices;
                end
                frame = inputROI.fluo2p.roi(roiN).frame;
                deleted = 0;
                sROI.roi(roiN) = struct('object',hROI,...
                    'ID',ID,'pos',pos,'XYvertices',ROIvertices,...
                    'label',text(ui.roiGUI.ax,...
                    'Position',inputROI.fluo2p.roi(roiN).label.Position,...
                    'String',inputROI.fluo2p.roi(roiN).label.String,'Color','g','FontWeight','bold'),...
                    'frame',frame,'deleted',deleted);
            end
            clear roiN
            
            %add position callback functions to each ROI so that they can
            %be repositioned and those positions will be saved
            for roiN = 1:length(sROI.roi)
                addNewPositionCallback(sROI.roi(roiN).object,...
                    @(pos) getpos(pos,sROI.roi(roiN).object,roiN));
            end
            clear roiN
            
            if isfield(inputROI.fluo2p,'baselineIDX')
                ui.roiGUI.baselineIDX = sROI.t(inputROI.fluo2p.baselineIDX);
            end
            if ~isfield(sROI,'frame')
                ui.roiGUI.frame = 1;
            end
        catch
            disp('Previous ROI set does not exist...')
        end
    end %loadPrevROIset

%This updates the ROI shape type when changed in settings
    function [] = roiTypeSelection(~,event)
        disp(['Current ROI Type: ' event.NewValue.String]);
        ui.roiGUI.curROItype = event.NewValue.String;
    end %roi radiobutton selection end

%This updates the plot with the respective channel
    function [] = channelSelection(~,event)
        disp(['Current channel: ' event.NewValue.String]);
        ui.roiGUI.curChannel = event.NewValue.String;
        if ~isfield(ui.roiGUI,'frame')
            ui.roiGUI.frame = 1;
        end
        
        %update plot image data to respective frame number
        if contains(ui.roiGUI.curChannel,'merge')
            ui.roiGUI.plt.CData = redGreenMerge(sROI.imgMultChan(1).img(:,:,ui.roiGUI.frame),...
                        sROI.imgMultChan(2).img(:,:,ui.roiGUI.frame));
        elseif strcmp(ui.roiGUI.curChannel,'ch. 1 (red)')
            ui.roiGUI.plt.CData = sROI.imgMultChan(1).img(:,:,ui.roiGUI.frame);
        elseif strcmp(ui.roiGUI.curChannel,'ch. 2 (green)')
            ui.roiGUI.plt.CData = sROI.imgMultChan(2).img(:,:,ui.roiGUI.frame);
        end
                
    end %roi radiobutton selection end

%function for drawing ROIs
    function [] = roi_btn_callbk(~,~)
        
        %draw ROI of shape given by settings
        disp('Draw ROI...')
        switch ui.roiGUI.curROItype
            case 'Ellipse'
                hROI = imellipse(ui.roiGUI.ax);
                ROIvertices = hROI.getVertices;
            case 'Free Hand'
                hROI = imfreehand(ui.roiGUI.ax);
                ROIvertices = hROI.getPosition;
        end
        
        hROI.Deletable = 'false';
        ROIpos = hROI.getPosition;
        %if there are existing ROIs append ROI structure data
        if isfield(sROI,'roi')
            roiIDX = length(sROI.roi)+1;
            ID = num2str(str2double(sROI.roi(end).ID)+1);
            sROI.roi(roiIDX).object = hROI;
            sROI.roi(roiIDX).ID = ID;
            sROI.roi(roiIDX).pos = ROIpos;
            sROI.roi(roiIDX).XYvertices = ROIvertices;
            sROI.roi(roiIDX).label = text(ROIpos(1),ROIpos(2)-5,ID,'Color','g','FontWeight','bold');
            sROI.roi(roiIDX).frame = ui.roiGUI.frame;
            sROI.roi(roiIDX).deleted = 0;
        else
            %create ROI structure if it's the first ROI
            if ~isfield(sROI,'frame')
                ui.roiGUI.frame = 1;
            end
            roiIDX = 1;
            sROI.roi = struct('object',hROI,...
                'ID','1','pos',ROIpos,'XYvertices',ROIvertices,...
                'label',text(ROIpos(1),ROIpos(2)-5,'1','Color','g','FontWeight','bold'),...
                'frame',ui.roiGUI.frame,'deleted',0);
        end
        
        if strmatch(class(sROI.roi(roiIDX).object),'imfreehand')
            sROI.roi(roiIDX).label.Position = [ROIpos(1,1) ROIpos(1,2)+5 0];
        end
        
        %create nested function for keeping the label in the right spot
        addNewPositionCallback(sROI.roi(roiIDX).object,@(pos) getpos(pos,sROI.roi(roiIDX).object,roiIDX));
    end %roi_btn_callbk

%This is a nested function called when ROI shapes are repositioned
%It will also move the label associated with the ROI
    function getpos(~, handle, roiIDX)
        pos = getPosition(handle);
        if strmatch(class(sROI.roi(roiIDX).object),'imfreehand')
            sROI.roi(roiIDX).label.Position = [pos(1,1) pos(1,2)+5 0];
        else
            sROI.roi(roiIDX).label.Position = [pos(1) pos(2)-5 0];
        end
        sROI.roi(roiIDX).pos = pos;
    end


%Allows the deletion of ROIs from list of current ROIs drawn on ROI gui
    function [] = delRoi_btn_callbk(~,~)
        if isfield(sROI,'roi')
            if ~isfield(sROI.roi,'deleted')
                sROI.roi = markDeletedROI(sROI.roi);
            end
            deletedROI = cell2mat({sROI.roi.deleted});
            allROI = {sROI.roi.ID};
            existingROIids = allROI(~deletedROI);
            
            %shows ROIs that haven't yet been deleted
            [delROI,~] = listdlg('PromptString','Select ROI(s)',...
                'ListString',existingROIids);
            delROIids = existingROIids(delROI);
            
            if ~isempty(delROIids)
                for roiN = 1:length(delROIids)
                    ID = delROIids(roiN);
                    delIDX = find(strcmp(ID,allROI));
                    sROI.roi(delIDX).object.delete;
                    sROI.roi(delIDX).deleted = 1;
                end
                clear roiN
                sROI.roi = markDeletedROI(sROI.roi);
            end
            
            %closes all figures other than settings and ROI figure
            h = findobj('-not','Name',['Choose frame and draw ROI'],'Type','Figure');
            if ~isempty(h)
                close(h)
            end
            
        else %if there are no ROIs available to delete
            disp('No ROIs are available to delete...')
        end
    end %delRoi_btn_callbk

%This function allows the re-loading of ROIs that have been deleted
%They will be shown again on the ROI gui at their last position
    function [] = loadDelROI_btn_callbk(~,~)
        if isfield(sROI,'roi')
            
            deletedROI = logical(cell2mat({sROI.roi.deleted}));
            if ~any(deletedROI)
                disp('No ROIs have been deleted yet...')
                return
            end
            
            allROI = {sROI.roi.ID};
            delROIids = allROI(logical(deletedROI));
            
            %Show list of ROIs that have been deleted
            [loadROI,~] = listdlg('PromptString','Select ROI(s)',...
                'ListString',delROIids);
            loadROIids = delROIids(loadROI);
            
            if ~isempty(loadROIids)
                %for each ROI chosen for re-loading initialize shape
                %function of that ROI on the ROI gui and create nested
                %repositioning function
                for roiN = 1:length(loadROIids)
                    ID = loadROIids(roiN);
                    loadIDX = find(strcmp(ID,allROI));
                    
                    type = class(sROI.roi(loadIDX).object);
                    
                    f = str2func(['@' type]);
                    pos = sROI.roi(loadIDX).pos;
                    
                    sROI.roi(loadIDX).object = f(ui.roiGUI.ax,pos);
                    
                    addNewPositionCallback(sROI.roi(loadIDX).object,@(pos)...
                        getpos(pos,sROI.roi(loadIDX).object,loadIDX));
                    
                    sROI.roi(loadIDX).deleted = 0;
                    labelpos = sROI.roi(loadIDX).label.Position;
                    labelstr = sROI.roi(loadIDX).label.String;
                    sROI.roi(loadIDX).label.delete;
                    sROI.roi(loadIDX).label = text(ui.roiGUI.ax,...
                        'Position',labelpos,...
                        'String',labelstr,'Color','g','FontWeight','bold');
                end
                clear roiN
                
                %Close all other figures besides the settings and ROI gui
                h = findobj('-not','Name',['Choose frame and draw ROI'],'Type','Figure');
                if ~isempty(h)
                    close(h)
                end
                
            else
                disp('No ROIs were chosen for loading...')
            end %if ROIs were chosen for re-loading
        else
            disp('No ROIs have been drawn or deleted...')
        end %if there are any ROIs
    end %load deleted ROI function callback


%loads the next tif file in the file numbering sequence
    function [] = loadNextCallback(~,~)
        %figure handle of current ROI gui figure / any other figures
        hPre = findall(0,'type','figure');
        
        %determine next file name
        % Capture the 5-digit group you want to increment
        tok = regexp(sROI.file, '[A-Z]{2}\d{4}[A-Z]{4}_(\d{5}).*\.tif$', 'tokens', 'once');
        if isempty(tok)
            error('Pattern not found in filename.');
        end
        
        nextFileNo = str2double(tok{1}) + 1;
        nextGroup = sprintf('%05d', nextFileNo);  % keep 5 digits; use %04d if 4 digits needed
        % Replace only that captured occurrence using regexprep with a function
        nextFileName = regexprep(sROI.file, '([A-Z]{2}\d{4}[A-Z]{4}_)(\d{5})(.*\.tif$)', ...
            ['$1' nextGroup '$3']);
        while nextFileNo<9999 && exist(nextFileName,'file')~=2
            nextFileNo = nextFileNo-1;
        end  
        nextGroup = sprintf('%05d', nextFileNo);  % keep 5 digits; use %04d if 4 digits needed
        % Replace only that captured occurrence using regexprep with a function
        nextFileName = regexprep(sROI.file, '([A-Z]{2}\d{4}[A-Z]{4}_)(\d{5})(.*\.tif$)', ...
            ['$1' nextGroup '$3']);
        
        %close current ROI gui figure and other figures
        close(hPre)
        
        %call this function on the next file
        roiGUI(nextFileName)
        
    end %loadnextCallbk


%loads the previous tif file in the file numbering sequence
    function [] = loadPrevCallback(~,~)
        %figure handle of current ROI gui figure / any other figures
        hPre = findall(0,'type','figure');
        
        %determine previous file name
        % Capture the 5-digit group you want to increment
        tok = regexp(sROI.file, '[A-Z]{2}\d{4}[A-Z]{4}_(\d{5}).*\.tif$', 'tokens', 'once');
        if isempty(tok)
            error('Pattern not found in filename.');
        end

        prevFileNo = str2double(tok{1}) - 1;
        prevGroup = sprintf('%05d', prevFileNo);  % keep 5 digits; use %04d if 4 digits needed
        % Replace only that captured occurrence using regexprep with a function
        prevFileName = regexprep(sROI.file, '([A-Z]{2}\d{4}[A-Z]{4}_)(\d{5})(.*\.tif$)', ...
            ['$1' prevGroup '$3']);
               
        while prevFileNo>0 && exist(prevFileName,'file')~=2
            prevFileNo = prevFileNo-1;
        end
        prevGroup = sprintf('%05d', prevFileNo);  % keep 5 digits; use %04d if 4 digits needed
        % Replace only that captured occurrence using regexprep with a function
        prevFileName = regexprep(sROI.file, '([A-Z]{2}\d{4}[A-Z]{4}_)(\d{5})(.*\.tif$)', ...
            ['$1' prevGroup '$3']);
        
        %close current ROI gui figure and other figures
        close(hPre)
        
        %call this function on the previous file
        roiGUI(prevFileName)
        
    end %loadPrevCallbk


%function to clear all ROIs and remove history for reloading
    function [] = clr_btn_callbk(~,~)
        if isfield(sROI,'roi')
            for roiN = 1:length(sROI.roi)
                if isvalid(sROI.roi(roiN).object)
                    sROI.roi(roiN).object.delete;    
                end
                sROI.roi(roiN).label.delete;
            end
            clear roiN
            %removes ROI history
            sROI = rmfield(sROI,'roi');
        end
        
        %close any other figures than settings and ROI gui
        closeROIplots()
    end %clr_btn_callbk


%function to close all figures other than settings and ROI gui
    function [] = close_fig_callbk(~,~)
        closeROIplots()
    end

%function to plot ROI output
    function [] = plot_btn_callbk(~,~)
        
        %closes all other figures than settings and ROI gui
        h = findobj('-not','Name',['Choose frame and draw ROI'],'Type','Figure');
        if ~isempty(h)
            close(h)
        end
        
        if isfield(sROI,'roi')
            sROI.roi = markDeletedROI(sROI.roi);
            sROI.roi = createMaskFromROIobject(sROI.roi);
            sROI.roi = getROIstructVertices(sROI.roi);
            sROI.rawFroi = roiMask2rawFroi(sROI.img,sROI.roi);

            %plot only those ROIs that haven't been deleted and find their
            %respective IDs
            validROI = ~cell2mat({sROI.roi.deleted});
            plotROIidx = find(validROI);
            roiIDs = {sROI.roi.ID};
            savedROIs = roiIDs(plotROIidx);
            
            if ~isfield(sROI,'state')
                answer = inputdlg({'Frame Rate (Hz)','dFF baseline start time (s)','dFF baseline end time (s)'},...
                    'Plot Param Input',[1 35],{'5','2','3'});
                sROI.frameRate = str2double(answer{1});
                sROI.t = 0:1/sROI.frameRate:(sROI.nFrames/sROI.frameRate)-1/sROI.frameRate;
                ui.roiGUI.dFFtBaseline(1) = str2double(answer{2});
                ui.roiGUI.dFFtBaseline(2) = str2double(answer{3});
            elseif ~isfield(ui.roiGUI,'dFFtBaseline')
                answer = inputdlg({'dFF baseline start time (s)','dFF baseline end time (s)'},...
                    'Plot Param Input',[1 35],{'2','3'});
                ui.roiGUI.dFFtBaseline(1) = str2double(answer{1});
                ui.roiGUI.dFFtBaseline(2) = str2double(answer{2});
            end
            
            %determine indices of whole tif trace baseline
            ui.roiGUI.baselineIDX = [find(sROI.t>round(ui.roiGUI.dFFtBaseline(1),1)-0.01 & sROI.t<round(ui.roiGUI.dFFtBaseline(1),1)+0.01)...
                find(sROI.t>round(ui.roiGUI.dFFtBaseline(2),1)-0.01 & sROI.t<round(ui.roiGUI.dFFtBaseline(2),1)+0.01)];
            dFF = dFoFcalc(sROI.rawFroi,ui.roiGUI.baselineIDX,1);
            
            %make plots for ROI traces across entire time duration
            %initialize average roi plot
            ui.outputFigs.meanROIoutput.fh = figure('Name',sROI.file);
            ui.outputFigs.meanROIoutput.ax1 = subplot(2,2,3);
            ui.outputFigs.meanROIoutput.ax3 = subplot(2,2,2);
            ui.outputFigs.meanROIoutput.ax4 = subplot(2,2,4);
            set(ui.outputFigs.meanROIoutput.fh,'Position',[100 100 755 540])
            
            %plot of current tif frame with all ROIs drawn
            imagesc(ui.outputFigs.meanROIoutput.ax1,mean(sROI.img,3))
            colormap(ui.outputFigs.meanROIoutput.ax1,gray)
            title(ui.outputFigs.meanROIoutput.ax1,sROI.file,'FontSize',10,'FontWeight','normal','Interpreter','none')
            
            %for each present ROI, draw on mean tif image
            for plotROIn = 1:length(plotROIidx)
                hold(ui.outputFigs.meanROIoutput.ax1,'on')
                plot(ui.outputFigs.meanROIoutput.ax1,sROI.roi(plotROIidx(plotROIn)).XYvertices(:,1),...
                    sROI.roi(plotROIidx(plotROIn)).XYvertices(:,2));
                %show label of ROI
                text(ui.outputFigs.meanROIoutput.ax1,...
                    'Position',sROI.roi(plotROIidx(plotROIn)).label.Position,...
                    'String',sROI.roi(plotROIidx(plotROIn)).label.String,'Color','g','FontWeight','bold')
            end
            ui.outputFigs.meanROIoutput.ax1.XLim = [0 size(sROI.img,2)];
            ui.outputFigs.meanROIoutput.ax1.YLim = [0 size(sROI.img,1)];
            axis(ui.outputFigs.meanROIoutput.ax1,'square')
            
            hLines = plot(ui.outputFigs.meanROIoutput.ax3, sROI.t, sROI.rawFroi(plotROIidx,:));
            xlabel(ui.outputFigs.meanROIoutput.ax3, 'Time (s)')
            ylabel(ui.outputFigs.meanROIoutput.ax3, 'Raw Average Pixel Intensity')
            % clickableLegend(hLines, savedROIs); 
            legend(hLines, savedROIs);
            lgd = findobj(gcf,'Type','Legend');
            lgd.Position(1:2) = [0.91 0.1];
            
            plot(ui.outputFigs.meanROIoutput.ax4,sROI.t(ui.roiGUI.baselineIDX(1):end),dFF(plotROIidx,:))
            hold on
            plot(ui.outputFigs.meanROIoutput.ax4,sROI.t(ui.roiGUI.baselineIDX(1):end),mean(dFF(plotROIidx,:)),'k-','LineWidth',3)
            xlabel(ui.outputFigs.meanROIoutput.ax4,'Time (s)')
            ylabel(ui.outputFigs.meanROIoutput.ax4,'\DeltaF/F')
            set(ui.outputFigs.meanROIoutput.ax1,'Position',[0.039735099337748,0.205555555555556,0.451655629139073,0.631481481481482]);
            
            %show vertical red bar at stim onset time
            if isfield(sROI,'TriggerParams') && isfield(sROI,'t')
                title(ui.outputFigs.meanROIoutput.ax4,...
                    ['Stim onset: ' num2str(sROI.TriggerParams.stimDelay) 's'],'FontWeight','normal')
                yl = ylim;
                plot(ui.outputFigs.meanROIoutput.ax4,...
                    [sROI.t(sROI.t == sROI.TriggerParams.stimDelay) sROI.t(sROI.t == sROI.TriggerParams.stimDelay)],...
                    [yl(1) yl(2)],'-r');
            end
            
            %add single pulse on top of plot
            if isfield(sROI,'Stim') && numel({sROI.Stim.pulsename})==1
                uicontrol('Style', 'text', 'String', sROI.Stim.pulsename, ...
                    'HorizontalAlignment', 'center', 'Units', 'normalized', ...
                    'Position', [0 .97 1 .03], 'BackgroundColor', [.88 .88 .88]);
            end
            
            %deal with mapping tif file multiple pulses
            %scrape frequency and amplitude of pulses, calculate dFF all pulses, and
            %output multiple pulse plots
            if size(sROI.Stim,2)>1
                [sROI.dFFptRel,sROI.rawFptRel,sROI.tPTrel,relPTonsetIDX] = ...
                    multPulseDFFcalcLocal(sROI.rawFroi,sROI.Stim,sROI.TriggerParams,sROI.frameRate);
                [PTfreq,PTdBampl,~,~,~,~,~] = extractMapPulseParams(sROI.Stim);
                nPulses = size(sROI.Stim,2);
                sROI.pulses.pulseNames = {sROI.Stim.pulsename}';
                sROI.pulses.pulseFreq = PTfreq;
                sROI.pulses.pulseDBampl = PTdBampl;
                
                %initialize subplots
                %determine number of plots needed for 9 pulses/plot
                pulsesPerFig = 9;
                remPulsePlotNo = rem(nPulses,pulsesPerFig);
                pulseFigNo = floor(nPulses/pulsesPerFig)+(remPulsePlotNo>=1);
                
                for pulseFigN = 1:pulseFigNo
                    ui.outputFigs.MultPulseOutput.delFoF(pulseFigN).fh = figure('Name',...
                        ['Multiple Pulses | delFoF For: ' sROI.file],...
                        'Position',[777 103 864 706]);
                    
                    %initialize subplots on multiple pulse figure for
                    %given pulse
                    for pulseSubPlotN = 1:pulsesPerFig
                        curPulseNo = pulseFigN*pulsesPerFig - pulsesPerFig + pulseSubPlotN;
                        if curPulseNo <= nPulses
                            figure(ui.outputFigs.MultPulseOutput.delFoF(pulseFigN).fh)
                            %NOTE SUBPLOT SIZE HARDCODED TO 3X3
                            ui.outputFigs.MultPulseOutput.delFoF(pulseFigN).(['ax' num2str(pulseSubPlotN)])...
                                = subplot(3,3,pulseSubPlotN);
                            xlabel('Time (s)')
                            ylabel('\DeltaF/F')
                            title([num2str(sROI.pulses.pulseFreq(curPulseNo)) ' Hz at '...
                                num2str(sROI.pulses.pulseDBampl(curPulseNo)) ' dB'])
                            hold(ui.outputFigs.MultPulseOutput.delFoF(pulseFigN).(['ax' num2str(pulseSubPlotN)]),'on')
                        end
                        clear curPulseNo
                    end %pulseSubPlotN
                    clear pulseSubPlotN
                end %pulseFigN 
                clear remPulsePlotNo pulseFigNo
                
                %roi traces for each pulse
                for plotROIn = 1:length(plotROIidx)
                    for pulseN = 1:length(sROI.pulses.pulseNames)
                        %PLOT
                        pulseFigN = ceil((pulseN/pulsesPerFig));
                        pulseSubPlotN = rem(pulseN,pulsesPerFig) +...
                            (rem(pulseN,pulsesPerFig)<1)*pulsesPerFig;
                        plot(ui.outputFigs.MultPulseOutput.delFoF(pulseFigN).(['ax' num2str(pulseSubPlotN)]),...
                            sROI.tPTrel(:,pulseN)',...
                            sROI.dFFptRel(plotROIidx(plotROIn),:,pulseN)');
                    end %for each pulse
                end %for each ROI
                
                %show red stim bars for multiple pulse plots
                for pulseN = 1:size(sROI.Stim,2)
                    pulseFigN = ceil((pulseN/pulsesPerFig));
                    pulseSubPlotN = rem(pulseN,pulsesPerFig) +...
                        (rem(pulseN,pulsesPerFig)<1)*pulsesPerFig;
                    
                    yl = ylim(ui.outputFigs.MultPulseOutput.delFoF(pulseFigN).(['ax' num2str(pulseSubPlotN)]));
                    plot(ui.outputFigs.MultPulseOutput.delFoF(pulseFigN).(['ax' num2str(pulseSubPlotN)]),...
                        [sROI.tPTrel(relPTonsetIDX,pulseN)...
                        sROI.tPTrel(relPTonsetIDX,pulseN)],...
                        [yl(1) yl(2)],'-r');
                    clear yl
                end
            end %multiple pulses
            
%             assignin('base','ui',ui)
            
        else
            disp('No ROIs available to plot...') %if no ROIs exist
        end %if there are ROIs in the data structure/drawn ROIs
        
    end %roiplot2P / plot button


    function [] = OutputWorkspace_btn_callbk(~,~)
        if isfield(sROI,'roi')
            %decided to re-scrape ROI info incase updates aren't tracked
            sROI.roi = markDeletedROI(sROI.roi);
            sROI.roi = createMaskFromROIobject(sROI.roi);
            sROI.roi = getROIstructVertices(sROI.roi);
            sROI.rawFroi = roiMask2rawFroi(sROI.img,sROI.roi);
        end
        assignin('base','sROI',sROI)
    end


    function [] = saveOutput_btn_callbk(~,~)
        %decided to re-scrape ROI info incase updates aren't tracked
        sROI.roi = markDeletedROI(sROI.roi);
        sROI.roi = createMaskFromROIobject(sROI.roi);
        sROI.roi = getROIstructVertices(sROI.roi);
        sROI.rawFroi = roiMask2rawFroi(sROI.img,sROI.roi);

        fluo2p = struct('file',sROI.file,...
            'nFrames',sROI.nFrames,...
            'imgWidth',sROI.imgWidth,...
            'imgHeight',sROI.imgHeight,...
            'roi',rmfield(sROI.roi,'label'),...
            'rawFroi',sROI.rawFroi);
        
        if isfield(sROI,'state')
            fluo2p.frameRate = sROI.frameRate;
            fluo2p.t = sROI.t;
            fluo2p.state = sROI.state;
        else
            disp('Not a ScanImage Tif, state and frameRate not saved.')
        end
        
        %for each roi
        for roiN = 1:length(sROI.roi)
            fluo2p.roi(roiN).type = class(sROI.roi(roiN).object);
            fluo2p.roi(roiN).label.String = sROI.roi(roiN).label.String;
            fluo2p.roi(roiN).label.Position = sROI.roi(roiN).label.Position;
        end
        fluo2p.roi = rmfield(fluo2p.roi,'object');
        
        %look for saved pulse and trigger info and if exsists, save in output structure
        %this should already have been scraped from directory, can just
        %grab from sROI STRUCT
        if isfield(sROI,'TriggerParams') && isfield(sROI,'Stim')
            fluo2p.Stim = sROI.Stim;
            fluo2p.TriggerParams = sROI.TriggerParams;
        end
        
        %write appended output file if one already exists
        [savepath,trace] = fileparts(sROI.file);
        if exist(fullfile(savepath,[trace '_roiOutput.mat']),'file')==2
            
            save(fullfile(savepath,[trace...
                '_roiOutput_' num2str(length(dir([savepath filesep trace '_roiOutput*']))+1)...
                '.mat']),'fluo2p')
        else
            save(fullfile(savepath,[trace '_roiOutput.mat']),'fluo2p')
        end
        
        disp('Done Saving.')
    end %save output callbk

%%%%%%%%%%%%%%% ADDITIONAL HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mark deleted ROI
    function roiStruct = markDeletedROI(roiStruct)
        validROI = arrayfun(@(r) isvalid(r.object),roiStruct);
        tmp = num2cell(~validROI)';
        [roiStruct.deleted] = tmp{:};
    end

%function to create masks
    function roiStruct = createMaskFromROIobject(roiStruct)
        validROI = ~cell2mat({roiStruct.deleted});
        masks = arrayfun(@(r) r.object.createMask,roiStruct(validROI),'uni',0);
        [roiStruct(validROI).mask] = masks{:};
    end

%function to get vertices for roiStruct
    function roiStruct = getROIstructVertices(roiStruct)
        for roiN = 1:length(roiStruct)
            if roiStruct(roiN).deleted==0
                if any(strmatch(class(roiStruct(roiN).object),'imellipse'))
                    roiStruct(roiN).XYvertices = roiStruct(roiN).object.getVertices;
                else
                    roiStruct(roiN).XYvertices = roiStruct(roiN).object.getPosition;
                end
            end
        end
    end


%function for calculating rawF from img and mask | decided whether to
%include deleted?
    function rawFroi = roiMask2rawFroi(img,roiStruct)
        nROI = length(roiStruct);
        nFrames = size(img,3);
        rawFroi = NaN(nROI,nFrames);
        
        for roiN = 1:length(roiStruct)
            if roiStruct(roiN).deleted==0
                for nFrame = 1:size(img,3)
                    f = img(:,:,nFrame);
                    rawFroi(roiN,nFrame) = mean(f(roiStruct(roiN).mask));
                    clear f
                end
            end
        end
        
    end

%function for outputting pulse synced dFF for multiple pulses in one
%rawFroi trace
    function [dFFptRel,rawFptRel,tPTrel,relPTonsetIDX] = multPulseDFFcalcLocal(rawFroi,stimPulseStruct,triggerParams,fs)
        [PTfreq,~,PTonsetInPulse,~,~,~,~] = extractMapPulseParams(stimPulseStruct);
        nCell = size(rawFroi,1);
        nPulses = length(PTfreq);
        framesPerPulse = triggerParams.ISI*fs;
        PTonsetIDX = PTonsetInPulse*fs;
        F = rawFroi(:,triggerParams.stimDelay*fs+1:triggerParams.stimDelay*fs+framesPerPulse*nPulses,:);
        F = reshape(F,nCell,framesPerPulse,nPulses);
        tAbs = 0:1/fs:(size(rawFroi,2)/fs)-1/fs;
        tAbs = tAbs(triggerParams.stimDelay*fs+1:triggerParams.stimDelay*fs+framesPerPulse*nPulses);
        tAbsTracePulse = reshape(tAbs,framesPerPulse,nPulses);
        
        %for PT relative traces:
        maxFramesAfterOnset = framesPerPulse-max(PTonsetIDX);
        maxFramesBeforeOnset = min(PTonsetIDX);
        relPTonsetIDX = maxFramesBeforeOnset+1;
        
        tPTrel = deal(zeros(maxFramesAfterOnset+maxFramesBeforeOnset,nPulses));
        [rawFptRel,dFFptRel] = deal(zeros(nCell,maxFramesAfterOnset+maxFramesBeforeOnset,nPulses));
        for pulseN = 1:nPulses
            tPTrel(:,pulseN) = tAbsTracePulse(PTonsetIDX(pulseN)-(maxFramesBeforeOnset-1):PTonsetIDX(pulseN)+maxFramesAfterOnset,pulseN);
            dFF = dFoFcalc(F(:,:,pulseN),[1 PTonsetIDX(pulseN)],1);
            rawFptRel(:,:,pulseN) = F(:,PTonsetIDX(pulseN)-(maxFramesBeforeOnset-1):PTonsetIDX(pulseN)+maxFramesAfterOnset,pulseN);
            dFFptRel(:,:,pulseN) = dFF(:,PTonsetIDX(pulseN)-(maxFramesBeforeOnset-1):PTonsetIDX(pulseN)+maxFramesAfterOnset);
            clear dFF
        end
    end

    %function to close plot figs
    function closeROIplots()
        if isfield(ui,'outputFigs')
            if isvalid(ui.outputFigs.meanROIoutput.fh)
            close(ui.outputFigs.meanROIoutput.fh)
            end
            if isfield(ui.outputFigs,'MultPulseOutput') && isvalid(ui.outputFigs.MultPulseOutput.delFoF(1).fh)
                for fN = 1:length(ui.outputFigs.MultPulseOutput.delFoF)
                    close(ui.outputFigs.MultPulseOutput.delFoF(fN).fh)
                end
            end
            
        end
    end




end %ENTIRE FUNCTION

