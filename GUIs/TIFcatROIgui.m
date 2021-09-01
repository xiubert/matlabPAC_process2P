function [] = TIFcatROIgui(TIFcat)

%draw ROIs on concatenated motion corrected image data
sROI.img = TIFcat;
sROI.nFrames = size(TIFcat,3);
sROI.imgHeight = size(TIFcat,1);

%this will be the main figure for drawing ROI
ui.roiGUI.fh = figure('Name','Choose frame and draw ROI','Visible','off');
ui.roiGUI.plt = imagesc(sROI.img(:,:,1));
ui.roiGUI.ax = findall(ui.roiGUI.fh,'type','axes');
ui.roiGUI.frameNoTxt = text(ui.roiGUI.ax,0.5,sROI.imgHeight-5,...
    'Frame: 1','Color','g');
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
    'Style', 'pushbutton', 'String', 'Abort and close all','Tag','exitbutton1',...
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

%show option buttons
ui.roiGUI.drawRoiBtn = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', 'Draw ROI','Tag','roibutton1',...
    'Units','Normalized',...
    'Position', [0.5+ui.roiGUI.ax.Position(1)-0.2, ...
    ui.roiGUI.ax.Position(2)+...
    ui.roiGUI.ax.Position(4),...
    0.1, 0.04],...
    'Callback',@roi_btn_callbk);
ui.roiGUI.loadRoiBtn = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', 'Load ROI Set','Tag','loadroibutton1',...
    'Units','Normalized',...
    'Position', [ui.roiGUI.ax.Position(1)-0.17, ...
    ui.roiGUI.ax.Position(4)+...
    ui.roiGUI.ax.Position(2)-0.04,...
    0.17, 0.04],...
    'Callback',@loadroi_btn_callbk);
ui.roiGUI.loadDelRoiBtn = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', 'Load Deleted ROI','Tag','roiHistButton1',...
    'Units','Normalized',...
    'Position', [ui.roiGUI.ax.Position(1)-0.17, ...
    ui.roiGUI.ax.Position(4)+...
    ui.roiGUI.ax.Position(2)-0.08,...
    0.17, 0.04],...
    'Callback',@loadDelROI_btn_callbk);
ui.roiGUI.meanToggle = uicontrol('Parent',ui.roiGUI.fh,'Style','togglebutton',...
    'String','Show mean image','Tag','meanImgToggle','Units','Normalized',...
    'Position',[ui.roiGUI.ax.Position(1)-0.17, ...
    ui.roiGUI.ax.Position(4)+...
    ui.roiGUI.ax.Position(2)-0.12,...
    0.17, 0.04],...
    'Callback',@meanImgToggle_callbk);
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
ui.roiGUI.saveRoiBtn = uicontrol('Parent',ui.roiGUI.fh,...
    'Style', 'pushbutton', 'String', ...
    '<html><b>Save ROI to workspace','Tag','donebutton1',...
    'Units','Normalized',...
    'Position', [ui.roiGUI.ax.Position(1)+...
    ui.roiGUI.ax.Position(3)+0.005, ...
    ui.roiGUI.ax.Position(4)+...
    ui.roiGUI.ax.Position(2)-0.32,...
    0.14, 0.04],...
    'Callback',@saveROI_btn_callbk);

%ROI Selection
%initialize GUI settings
ui.roiGUI.curROItype = 'Ellipse'; %default ROI on start
ui.roiGUI.roiType = uicontrol('Parent',ui.roiGUI.fh,'Style', 'text',...
    'String','Choose ROI Type:','Units','Normalized',...
    'Position',[0.029245283018868,0.203260869565217,0.091509433962264,0.022826086956522]);
ui.roiGUI.bg = uibuttongroup('Parent',ui.roiGUI.fh,...
    'Visible','off',...
    'Position',[0.022641509433962,0.13804347826087,0.109433962264151,0.069565217391304],...
    'SelectionChangedFcn',@roiTypeSelection);

% Create three radio buttons in the button group.
ui.roiGUI.ellipse1 = uicontrol(ui.roiGUI.bg,'Style',...
    'radiobutton',...
    'String','Ellipse',...
    'Units','Normalized',...
    'Position',[0.205701791148595,0.587139024032356,0.870699999999998,0.1997]);
ui.roiGUI.freehand1 = uicontrol(ui.roiGUI.bg,'Style','radiobutton',...
    'String','Free Hand',...
    'Units','Normalized',...
    'Position',[0.211428571428571,0.158122830089119,0.870700000000001,0.1997]);
% % Make the uibuttongroup visible after creating child objects.
set(ui.roiGUI.fh,'Visible','on')
ui.roiGUI.bg.Visible = 'on';

%% CALLBACKS

%CALLBACKS
%closes everything without saving; leaves variables in workspace
    function [] = exit_btn_callbk(varargin)
        close all
    end

%slider for going through image frames
    function [] = slider_callback(hObject,~)
        %get frame number from slider position
        frame = round(get(hObject,'value'));
        
        %set frame number in GUI data so it is known which frame ROI was
        %drawn on
        ui.roiGUI.frame = frame;
        
        %update plot image data to respective frame number
        ui.roiGUI.plt.CData = sROI.img(:,:,frame);
        
        %change frame number on plot
        ui.roiGUI.frameNoTxt.String = ['Frame: ' num2str(frame)];
        
    end %slider_callback

%This updates the ROI shape type when changed in settings
    function [] = roiTypeSelection(~,event)
        disp(['Current ROI Type: ' event.NewValue.String]);
        ui.roiGUI.curROItype = event.NewValue.String;
    end %roi radiobutton selection end

%function for drawing ROIs
    function [] = roi_btn_callbk(hObject,~)
        %draw ROI of shape given by settings
        disp('Draw ROI...')
        switch ui.roiGUI.curROItype
            case 'Ellipse'
                hROI = imellipse(ui.roiGUI.ax);
                ROIpos = hROI.getPosition;
                ROIvertices = hROI.getVertices;
            case 'Free Hand'
                hROI = imfreehand(ui.roiGUI.ax);
                ROIpos = hROI.getPosition;
                ROIvertices = hROI.getPosition;
        end
        hROI.Deletable = 'false';

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
            if ~isfield(ui.roiGUI,'frame')
                ui.roiGUI.frame = 1;
            end
            roiIDX = 1;
            sROI.roi = struct('object',hROI,....
                'ID','1','pos',ROIpos,'XYvertices',ROIvertices,...
                'label',text(ROIpos(1),ROIpos(2)-5,'1','Color','g','FontWeight','bold'),....
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

%function to load ROIs from _roiOutput.mat file
%keeps deleted ROIs, but they are still deleted unless load deleted ROI
%button is chosen
    function [] = loadroi_btn_callbk(hObject,~)
        [inputROIfile,roiFilepath] = uigetfile('*.mat',...
            'Select output file containing ROI locations');
        if isempty(inputROIfile)
            disp('No File Selected...')
            return
        end
        
        %handle different roi file types
        if contains(inputROIfile,'roiOutput')
            %load fluo2p structure from _roiOutput.mat file
            inputROI = load([roiFilepath inputROIfile]);
            
            %for each ROI in roiOutput.mat file
            for roiN = 1:length(inputROI.fluo2p.roi)
                
                %get draw shape function from ROI shape type
                f = str2func(['@' inputROI.fluo2p.roi(roiN).type]);
                
                %call shape function on current ROI gui
                hROI = f(ui.roiGUI.ax,inputROI.fluo2p.roi(roiN).pos);
                
                %get ROI info
                ID = inputROI.fluo2p.roi(roiN).ID;
                pos = inputROI.fluo2p.roi(roiN).pos;
                frame = inputROI.fluo2p.roi(roiN).frame;
                deleted = inputROI.fluo2p.roi(roiN).deleted;
                
                %save ROI info to sROI structure
                sROI.roi(roiN) = struct('object',hROI,...
                    'ID',ID,'pos',pos,...
                    'label',text(ui.roiGUI.ax,...
                    'Position',inputROI.fluo2p.roi(roiN).label.Position,...
                    'String',inputROI.fluo2p.roi(roiN).label.String,'Color','g','FontWeight','bold'),...
                    'frame',frame,'deleted',deleted);
                
                %if the ROI was deleted don't show it on the ROI figure
                if inputROI.fluo2p.roi(roiN).deleted == 1
                    sROI.roi(roiN).object.delete;
                end
            end
            clear roiN
            
        elseif contains(inputROIfile,'moCorrROI')
            %load moCorrROI
            S = load([roiFilepath inputROIfile],'moCorROI');
            inputROI = S.moCorROI;
            clear S
            %for each ROI in file
            for roiN = 1:length(inputROI)
                %get draw shape function from ROI shape type
                f = str2func(['@' inputROI(roiN).type]);
                
                %call shape function on current ROI gui
                hROI = f(ui.roiGUI.ax,inputROI(roiN).pos);
                
                %get ROI info
                ID = inputROI(roiN).ID;
                pos = inputROI(roiN).pos;
                frame = inputROI(roiN).frame;
                deleted = inputROI(roiN).deleted;
                
                %save ROI info to sROI structure
                sROI.roi(roiN) = struct('object',hROI,...
                    'ID',ID,'pos',pos,...
                    'label',text(ui.roiGUI.ax,...
                    'Position',inputROI(roiN).label.Position,...
                    'String',inputROI(roiN).label.String,'Color','g','FontWeight','bold'),...
                    'frame',frame,'deleted',deleted);
                
                %if the ROI was deleted don't show it on the ROI figure
                if inputROI(roiN).deleted == 1
                    sROI.roi(roiN).object.delete;
                end
            end
        else
            disp('No File Selected...')
        end
        
        %add position callback functions to each ROI so that they can
        %be repositioned and those positions will be saved
        for roiN = 1:length(sROI.roi)
            if sROI.roi(roiN).deleted == 0
                addNewPositionCallback(sROI.roi(roiN).object,...
                    @(pos) getpos(pos,sROI.roi(roiN).object,roiN));
            end
        end
        clear roiN
        
        
        if ~isfield(ui.roiGUI,'frame')
            ui.roiGUI.frame = 1;
        end
    end %load_ROI_callback

%This function allows the re-loading of ROIs that have been deleted
%They will be shown again on the ROI gui at their last position
    function [] = loadDelROI_btn_callbk(hObject,~)
        if isfield(sROI,'roi')
            deletedROI = extractfield(sROI.roi,'deleted');
            assignin('base','deletedROI',deletedROI)
            if iscell(deletedROI)
                deletedROI = [deletedROI{:}];
            end
%             deletedROI
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
                
                %Close all other figures besides the ROI gui
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

%Allows the deletion of ROIs from list of current ROIs drawn on ROI gui
    function [] = delRoi_btn_callbk(hObject,~)
        if isfield(sROI,'roi')
            deletedROI = extractfield(sROI.roi,'deleted');
%             deletedROI = cellfun(@any,deletedROI);
            if iscell(deletedROI)
                deletedROI = [deletedROI{:}];
            end
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
            end
            
            %closes all figures other than ROI figure
            h = findobj('-not','Name',['Choose frame and draw ROI'],'Type','Figure');
            if ~isempty(h)
                close(h)
            end
            
        else %if there are no ROIs available to delete
            disp('No ROIs are available to delete...')
        end
    end %delRoi_btn_callbk

%function to clear all ROIs and remove history for reloading
    function [] = clr_btn_callbk(hObject,~)
        
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
        
        %close any other figures than ROI gui
        h = findobj('-not','Name',['Choose frame and draw ROI'],'Type','Figure');
        if ~isempty(h)
            close(h)
        end
        
    end %clr_btn_callbk

%show mean of frames across time
    function [] = meanImgToggle_callbk(hObject,~)
        if get(hObject,'value')==1
            ui.roiGUI.slider.Visible = 0;
            ui.roiGUI.plt.CData = nanmean(sROI.img,3);
        elseif get(hObject,'value')==0
            ui.roiGUI.slider.Visible = 1;
            try
                ui.roiGUI.plt.CData = sROI.img(:,:,ui.roiGUI.frame);
            catch
                ui.roiGUI.plt.CData = sROI.img(:,:,1);
            end
        end
    end %meanImgToggle_callbk

%save ROIs
    function [] = saveROI_btn_callbk(hObject,~)
        
        %closes all other figures than ROI gui
        h = findobj('-not','Name',['Choose frame and draw ROI'],'Type','Figure');
        if ~isempty(h)
            close(h)
        end
        
        if isfield(sROI,'roi')
            
            %determine valid/still present ROIs
            validROI = arrayfun(@(r) isvalid(r.object),sROI.roi);
            tmp = num2cell(~validROI)';
            [sROI.roi.deleted] = tmp{:};
            
            %save only those ROIs that haven't been deleted and find their
            %respective IDs
            savROIidx = find(validROI);
            moCorROI = rmfield(sROI.roi(savROIidx),{'label','object'});
            
            %for each ROI that has not been deleted:
            for Nroi = 1:length(savROIidx)
                %get image mask of ROI
                moCorROI(Nroi).type = class(sROI.roi(savROIidx(Nroi)).object);
                moCorROI(Nroi).mask = sROI.roi(savROIidx(Nroi)).object.createMask;
                if strmatch(moCorROI(Nroi).type,'imellipse')
                    moCorROI(Nroi).XYvertices = sROI.roi(savROIidx(Nroi)).object.getVertices;
                else
                    moCorROI(Nroi).XYvertices = sROI.roi(savROIidx(Nroi)).object.getPosition;
                end
                moCorROI(Nroi).ROIxyCoord = mask2polyCoord(moCorROI(Nroi).mask);
                moCorROI(Nroi).ROIcurveOrderedXY = orderEllipsePtOnCurve(moCorROI(Nroi).ROIxyCoord);
                moCorROI(Nroi).label.String = sROI.roi(savROIidx(Nroi)).label.String;
                moCorROI(Nroi).label.Position = sROI.roi(savROIidx(Nroi)).label.Position;
            end
            assignin('base','moCorROI',moCorROI);
            disp('"moCorROI" saved to workspace...')
            disp('Done.')
            close all
        end
    end %saveROI btn callbk

end