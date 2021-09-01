function tifFiles = moCorRawF2tifList(tifFiles,moCorrImgNonRigid,moCorROI,rawCatImg)
%gets roiF from corrected and uncorrected tif stacks and adds to tifFiles
%structure
moCorRawFroi = cell(length(tifFiles),1);
if nargin==4
    rawFroi = cell(length(tifFiles),1);
end

endFrame = 0;
for fileNo = 1:length(tifFiles)
    %find start and end frames corresponding to tif file
    [~,imHeader,~,~] = ...
    readSCIMtif([tifFiles(fileNo).folder filesep tifFiles(fileNo).name],'metaOnly');
    tifFiles(fileNo).nFrames = imHeader.nFrames;
    tifFiles(fileNo).frameRate = imHeader.frameRate;
    startFrame = endFrame+1;
    endFrame = startFrame+(imHeader.nFrames-1);
    
    %initialize vars
    moCorRawFroi{fileNo} = zeros(length(moCorROI),imHeader.nFrames);
    if nargin==4
        rawFroi{fileNo} = zeros(length(moCorROI),imHeader.nFrames);
    end
    
    %for each roi get roiF at each of the corresponding frames for given
    %tif
    for nROI = 1:length(moCorROI)
        traceFrame = 0;
        for nFrame = startFrame:endFrame
            traceFrame = traceFrame+1;
            
            fMoCor = moCorrImgNonRigid(:,:,nFrame);
            moCorRawFroi{fileNo}(nROI,traceFrame) = mean(fMoCor(moCorROI(nROI).mask));
            if nargin==4
                fImg = rawCatImg(:,:,nFrame);
                rawFroi{fileNo}(nROI,traceFrame) = mean(fImg(moCorROI(nROI).mask));
                clear fImg
            end
            clear fMoCor
        end
    end
    
    clear startFrame tmpImInfo
end

%fill structure w/ roiF data
if nargin==4
    [tifFiles(:).rawFroi] = rawFroi{:};
end
[tifFiles(:).moCorRawFroi] = moCorRawFroi{:};


