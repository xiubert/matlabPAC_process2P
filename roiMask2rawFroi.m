function rawFroi = roiMask2rawFroi(img,roiStruct)
    nROI = length(roiStruct);
    nFrames = size(img,3);
    rawFroi = NaN(nROI,nFrames);

    for roiN = 1:length(roiStruct)
        if roiStruct(roiN).deleted == false
            for nFrame = 1:size(img,3)
                f = img(:,:,nFrame);
                rawFroi(roiN,nFrame) = mean(f(roiStruct(roiN).mask));
                clear f
            end
        end
    end

end