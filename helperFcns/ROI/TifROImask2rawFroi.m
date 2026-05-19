function [rawFroi,nFrames,varargout] = TifROImask2rawFroi(img,roiStruct,channel)

if nargin < 3
    channel = 2;
end

%process tif file if input is not an image matrix
if ischar(img) || isstring(img)
    [img,imHeader] = readSCIMtif(img,channel);
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