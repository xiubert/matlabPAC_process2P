function [rawFroi,nFrames,varargout] = TifROImask2rawFroi(img,roiStruct)

%process tif file if input is not an image matrix
if ischar(img) || isstring(img)
    %always load channel 2 containing green channel
    [img,imHeader] = readSCIMtif(img,2);
    nFrames = imHeader.nFrames;
    varargout{1} = round(imHeader.frameRate);
else
    nFrames = size(img,3);
end

nROI = length(roiStruct);
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