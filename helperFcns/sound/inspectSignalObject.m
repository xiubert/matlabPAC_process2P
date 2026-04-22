function varargout = inspectSignalObject(varargin)
% INSPECTSIGNALOBJECT  Load and plot signal amplitude over time from an Ephus .signal file.
%
%   inspectSignalObject()
%   inspectSignalObject(Name, Value, ...)
%   y = inspectSignalObject(...)
%   [y, signalSampleRate] = inspectSignalObject(...)
%   [y, signalSampleRate, sigFilePath] = inspectSignalObject(...)
%
%   With no arguments, opens a file browser to select one or more .signal
%   files.  For DRC stimuli, colors the plot by contrast level using
%   getContrastDRCvars.
%
%   NOTE: the Ephus @signalObject class library must be on the MATLAB path
%   (e.g. ephus_library).
%
%   Name-Value Inputs:
%     'signalPath'  - path to a .signal file (default: false, opens browser)
%     'plot'        - logical; plot signal (default: true)
%     'plotNearPT'  - logical; for DRC stim files, zoom plot to ±0.5/1.5 s
%                     around the pure-tone onset (default: false)
%
%   Outputs (all optional):
%     y                - signal vector (amplitude, in raw units / 1000 = V)
%     signalSampleRate - sample rate in Hz
%     sigFilePath      - full path to the loaded .signal file
%
%   See also getContrastDRCvars

p = inputParser;
validationFcnPath = @(x) isfile(x) || islogical(x);
validationFcnPlot = @(x) islogical(x);

addParameter(p,'signalPath',false,validationFcnPath);
addParameter(p,'plot',true,validationFcnPlot);
addParameter(p,'plotNearPT',false,validationFcnPlot);

parse(p,varargin{:});

if ~p.Results.signalPath
    [sigfile, sigpath] = uigetfile('C:\Data\Rig Software\250kHzPulses\*.signal','Choose signal to analyze...','MultiSelect','on');
else
    [path, file, ext] = fileparts(p.Results.signalPath);
    sigpath = path;
    sigfile = [file ext];
end

if ~iscell(sigfile)
    tmp = sigfile;
    clear sigfile
    sigfile{1} = tmp;
end

for fileNo = 1:length(sigfile)
    sigS = load(fullfile(sigpath,sigfile{fileNo}),'-mat');
    signalduration = get(sigS.signal,'length'); %in seconds
    signalSampleRate = get(sigS.signal,'sampleRate'); %in Hz
    time = 0:1/signalSampleRate:signalduration-1/signalSampleRate;
    disp(sigfile{fileNo})
    %get time series data
    y = getdata(sigS.signal,signalduration); %or get(signal,'signal');
    
    if p.Results.plot
        figure('Name',strrep(sigfile{fileNo},'.signal',''))
        
        if contains(sigfile{fileNo},'DRC')
            %determine contrast level
            [colors,~] = getContrastDRCvars();
            dBspan = regexp(sigfile{fileNo},'(?<a>\d{2}-\d{2})dB','tokens','once');
            dBdelta = abs(eval(string(dBspan{1})));
            pColor = colors.lohiPre((dBdelta>15)+1,:);
            
            
            if contains(sigfile{fileNo},'stim') && ~p.Results.plotNearPT
                plot(time,y/1000,'Color',pColor);
                xlim([-2 12])
                hold on
                plot([-2 12],[0 0],'-','Color',pColor)
                set(gcf,'Position',[301 321 1054 200])
                
            elseif contains(sigfile{fileNo},'stim') && p.Results.plotNearPT
                ptOnset = str2double(strrep(string(regexp(sigfile{fileNo},...
                    '_at_(?<a>(\dpt\d|\d))s','tokens')),'pt','.'));
                
                startT = find(time==ptOnset-0.5);
                endT = find(time==ptOnset+1.5);
                
                plot(time(startT:endT),y(startT:endT)/1000,'Color',pColor);
            else
                startT = find(time==1.5);
                endT = find(time==3.5);
                plot(time(startT:endT),y(startT:endT)/1000,'Color',pColor);
                %                         startT = find(time==1.5);
                %             endT = find(time==3.5);
                plot(time,y/1000,'Color',pColor);
            end
            
        else
            plot(time,y/1000);
            yl = ylim;
            if yl(2)==1
                ylim([-1.2 1.2])
            end
        end
        
        ylabel('Voltage (V)')
        xlabel('Time (s)')
        modPlotForPaper(false);
        %     title(strrep(sigfile{fileNo},'.signal',''),'Interpreter','none','FontWeight','normal','FontSize',8)
    end
end
% clear signalduration signalSampleRate time y sigS

%handle output
switch nargout
    case 0
        varargout{1} = '';
    case 1
        varargout{1} = y;
    case 2
        varargout{1} = y;
        varargout{2} = signalSampleRate;
    case 3
        varargout{1} = y;
        varargout{2} = signalSampleRate;
        varargout{3} = fullfile(sigpath,sigfile);
end


end
