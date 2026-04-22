function RGBout = redGreenMerge(redImg,greenImg)
% REDGREENMERGE  Merge two grayscale images into a red-green RGB composite.
%
%   RGBout = redGreenMerge(redImg, greenImg)
%
%   Inputs:
%     redImg   - 2D grayscale image assigned to the red channel
%     greenImg - 2D grayscale image assigned to the green channel;
%                must be the same size as redImg
%
%   Output:
%     RGBout - MxNx3 RGB image with redImg in channel 1, greenImg in
%              channel 2, and NaN (black) in the blue channel; both
%              input channels are scaled to [0, 1] before merging
    if size(redImg)~=size(greenImg)
        error('Image sizes must be the same')
    end
    RGBout = cat(3,scaleZeroToOne(redImg),scaleZeroToOne(greenImg),nan(size(redImg)));