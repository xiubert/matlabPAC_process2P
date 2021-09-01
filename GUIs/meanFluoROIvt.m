function [rawintensityROI, movie, header] = meanFluoROIvt(varargin)

p = inputParser;
addParameter(p,'file','',@(x) ischar(x))
addParameter(p,'baseline',[2 3],@(x) isnumeric(x) && isvector(x) && length(x)==2);
addParameter(p,'filtOrder',4,@(x) isnumeric(x));
addParameter(p,'filtCutoffFreq',5,@(x) isnumeric(x));
addParameter(p,'stimlen',0.4,@(x) isnumeric(x));
addParameter(p,'fr',20,@(x) isnumeric(x));
addParameter(p,'temporalAvgFrameWindow',10,@(x) isnumeric(x));

parse(p,varargin{:});

file = p.Results.file;

if isempty(file)
    [fName,fDir] = uigetfile({'*.qcamraw';'*.tif'},'Load fluorescence movie');
    file = fullfile(fDir,fName);
end
if exist(file,'file')~=2
    error('...File doesn''t exist...')
end

%determine filetype
[~,~,ext] = fileparts(file);

if strmatch('.tif',ext) %or strncmp('.tif',ext,length(ext))
    imginfo = imfinfo(file);
    nFrames = length(imginfo);
    imgWidth = imginfo(1).Width;
    imgHeight = imginfo(1).Height;
    
    warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
    
    curtiff = Tiff(file,'r');
    
    header = strsplit(imginfo(1).ImageDescription,'\r');
    if length(header)>1 %scanimage v~3
        
        for k = 1:length(header)
            if ~isempty(header{k})
                evalc(header{k});
            end
        end
        clear header
        
        assignin('base','state',state);
        fr = round(state.acq.frameRate);
        
    else % %scanimage v~5
        tmp = regexp(imginfo(1).Software,'SI.hRoiManager.scanFrameRate = (?<fr>\d{1,4})','names');
        fr = str2double(tmp.fr);
    end
    
    for k = 1:nFrames
        curtiff.setDirectory(k);
        img(:,:,k) = curtiff.read();
    end
    curtiff.close();
    
    warning('on','MATLAB:imagesci:tiffmexutils:libtiffWarning')
    
    img = double(img);
    movie = img;
    t = [1:nFrames]*(1/fr);
    
    %ROI
    h1 = figure('Name',['Spatial Fluorescence:  DRAW ROI for ' file]);
    imagesc(img(:,:,1));
    colormap(gray)
    
    disp('Draw ROI')
    rect = round(getrect); %output: [X1 Y1 width height] (top left, W H)
    roi = [rect(2) rect(4)+rect(2) rect(1) rect(3)+rect(1)]; %as row1 row2 col1 col2
    disp('Done')
    close(h1)
    
elseif strmatch('.qcamraw',ext)
    fr = p.Results.fr;
    fid = fopen(file,'r');
    %get line length of qcamraw file header (goes until it finds an empty line)
    headcount = 0;
    while ~isempty(fgetl(fid))
        headcount = headcount+1;
    end
    frewind(fid)
    
    %fill header fieldnames with field data
    for k = 1:headcount
        line = fgetl(fid);
        colonLoc = regexp(line, ':');
        spaces = regexp(line, ' ');
        
        if strfind(line,'-') %fix for incompatible fieldnames
            line(1:colonLoc-1) = strrep(line(1:colonLoc-1),'-','_');
        end
        
        %fill fieldname field with field data
        if isequal(colonLoc(1),spaces(end)-1)
            header.(line(1:colonLoc-1)) = line(colonLoc+2:end);
        elseif numel(spaces)>2 %for ROI field
            header.(line(1:colonLoc-1)) = line(colonLoc+2:end);
        else %for fields with unit description in brackets
            header.(line(1:colonLoc-1)) = line(colonLoc+2:spaces(end)-1);
        end
    end
    
    %so we know what bytes to actually use for the image
    headerSize = str2num(header.Fixed_Header_Size);
    
    % set file pointer to EOF
    fseek(fid,0,1);
    % get the number of bytes from start of file
    numBytes = ftell(fid);
    % set file pointer back at start of file
    frewind(fid)
    
    imgvec=fread(fid,numBytes,'uint16','ieee-le'); %results in double presicion column vector of image data
    
    fclose(fid);
    fclose('all');
    
    %establish pixel width and height from file header ROI field
    totROI = str2num(header.ROI);
    imgWidth = totROI(3);
    imgHeight = totROI(4);
    
    %get frames from read file and ensure they were read as they are in header
    nFrames = (length(imgvec)-headerSize/2)/imgWidth/imgHeight;
    %should equal:  (numBytes-headerSize)/str2num(header.Frame_Size)
    if nFrames ~= (numBytes-headerSize)/str2num(header.Frame_Size)
        disp('Something went wrong w/r/t file size and pixel depth')
        return
    end
    
    img = reshape(imgvec((headerSize/2)+1:end), [imgWidth, imgHeight, nFrames]);
    img = permute(img,[2 1 3]); %ensure image matrix is ordered as [height, width, frame]
    movie = img;
    t = [1:nFrames]*(1/p.Results.fr);
    
    %spatial fluorescence
    normCutoff = p.Results.filtCutoffFreq/(fr/2);
    framespan = p.Results.temporalAvgFrameWindow;

    baselineIDX = (t>=p.Results.baseline(1) & t<=p.Results.baseline(2));

    spatialbase = mean(img(:,:,baselineIDX),3);
    dFpixel = bsxfun(@minus,img,spatialbase);
    dFFpixel = bsxfun(@rdivide,dFpixel,spatialbase);
    
    %butterworth filter
    [a,b] = butter(p.Results.filtOrder, normCutoff, 'low'); %low pass butterworth filter
    filteredDelFoF=filtfilt(a,b,dFFpixel);
    
    %temporal average of framespan consecutive frames immediately following sound
    %stimulus
    %assumes stim immediately proceeds baseline(2)
    spatialDelFoF = mean(filteredDelFoF(:,:,...
        (t>=(p.Results.baseline(2)+p.Results.stimlen) & t<=(p.Results.baseline(2)+p.Results.stimlen+framespan*(1/fr)))...
        ),3);
    
    h1 = figure('Name',file);
    imagesc(img(:,:,1));
    colormap(gray)
%     movegui(h1,'southwest')

    
    h2 = figure('Name',['Spatial Fluorescence:  DRAW ROI for ' file]);
%     movegui(h2,'northwest')
    imagesc(spatialDelFoF)
    colormap(jet)
    
    disp('Draw ROI')
    rect = round(getrect); %output: [X1 Y1 width height] (top left, W H)
    %THIS WAS WRONG:  roi = [rect(2) rect(3)+rect(2) rect(1) rect(4)+rect(1)]; %as row1 row2 col1 col2
    roi = [rect(2) rect(4)+rect(2) rect(1) rect(3)+rect(1)]; %as row1 row2 col1 col2
    disp('Done')
    close(h1)
    close(h2)
else
    error('...This function does not handle the provided file type...')
end

%spatial mean rawF within ROI
imgroi = img(roi(1):roi(2),roi(3):roi(4),:);
imgroi2 = reshape(imgroi,(roi(2)+1-roi(1))*(roi(4)+1-roi(3)),nFrames);
rawintensityROI = mean(imgroi2,1);

h1 = figure('Name',file);
set(h1,'Position',[572 292 755 540])
ax1 = subplot(2,2,1);
imagesc(img(:,:,1));
colormap(ax1,gray)

if strmatch('.tif',ext)
    rectangle('Position',rect,'EdgeColor','r')
end

if strmatch('.qcamraw',ext)
    ax2 = subplot(2,2,2);
    imagesc(spatialDelFoF);
    if strmatch('2013',version('-release'))
        cmap1 = gray;
        cmap2 = jet(size(cmap1,1));
        cmap = [cmap1;cmap2];
        colormap(cmap);
        
%         colormap(ax2,gray)
    else
        colormap(ax2,jet)
    end
    rectangle('Position',rect)
elseif strmatch('.tif',ext)
    ax2 = subplot(2,2,2);
    roiy = roi(1):roi(2);
    roix = roi(3):roi(4);
    imagesc(imgroi(:,:,1));
    colormap(ax2,gray)
    yticks(find(ismember(roiy,roi(1):5:roi(2))))
    yticklabels(string(roiy(ismember(roiy,roi(1):5:roi(2)))))
    xticks(find(ismember(roix,roi(3):5:roi(4))))
    xticklabels(string(roix(ismember(roix,roi(3):5:roi(4)))))
end

subplot(2,2,3)
plot(t,rawintensityROI)
xlabel('Time (s)')
ylabel('Raw Average Pixel Intensity')

subplot(2,2,4)
delF = rawintensityROI./mean(rawintensityROI(t>=p.Results.baseline(1) & t<=p.Results.baseline(2)));
plot(t,delF)
xlabel('Time (s)')
ylabel('\DeltaF')
end
