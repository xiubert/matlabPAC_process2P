function writeMoCorTifs(files,moCorImgData,outputPath)
%prevents printing of irrelevant tif header errors
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');

fStart = 1;
for fileNo = 1:length(files)
    tmpImInfo = imfinfo([files(fileNo).folder filesep files(fileNo).name]);
    nFrames = length(tmpImInfo);
    fEnd = fStart+(nFrames-1);
    
    %ensure header is written for scim_openTif compatibility (Scanimage)
    hInputTif = Tiff([files(fileNo).folder filesep files(fileNo).name]);
    tagList = hInputTif.getTagNames;

    tifWrite = Tiff([outputPath filesep strrep(files(fileNo).name,'.tif','_NoRMCorre.tif')],'w');
    
    frameNo = 0;
    for k = fStart:fEnd
        frameNo = frameNo+1;
        for tagID = 1:length(tagList)
            if ~any(contains({'SubFileType','StripOffsets','YCbCrSubSampling',...
                    'NumberOfInks','MinSampleValue','MaxSampleValue'},tagList{tagID}))
                try
                    setTag(tifWrite,tagList{tagID},hInputTif.getTag(tagList{tagID}));
                catch
                end
            end
        end
        setTag(tifWrite,'BitsPerSample',tmpImInfo(1).BitsPerSample);
        setTag(tifWrite,'SampleFormat',Tiff.SampleFormat.Int);
        setTag(tifWrite,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
        
        write(tifWrite,int16(moCorImgData(:,:,k)))
        if frameNo<nFrames
            tifWrite.writeDirectory;
        end        
    end
   
    hInputTif.close();
    tifWrite.close();
    
    fStart = fStart+nFrames;
    
    clear nFrames hTif
end
warning('on','MATLAB:imagesci:tiffmexutils:libtiffWarning')
