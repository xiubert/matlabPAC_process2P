function justImgData = justLoadTif(TifFile)
% JUSTLOADTIF  Load all frames from a single-channel .tif, ignoring metadata.
%
%   justImgData = justLoadTif(TifFile)
%
%   Reads every directory (frame) of a .tif using the low-level Tiff API,
%   suppressing libtiff header warnings.  Intended for ScanImage .tif files
%   when metadata is not needed.  Does not support multi-channel .tif files.
%
%   Input:
%     TifFile     - path to a single-channel .tif file
%
%   Output:
%     justImgData - Width x Height x nFrames array of raw pixel values
%
%   See also readSCIMtif, splitTifChans
tifInfo = imfinfo(TifFile);
nFrames = length(tifInfo);
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
warning('off','imageio:tiffmexutils:libtiffWarning');
curtiff = Tiff(TifFile,'r');
justImgData = zeros(tifInfo(1).Width,tifInfo(1).Height,nFrames);

for k = 1:nFrames
    curtiff.setDirectory(k);
    justImgData(:,:,k) = curtiff.read();
end
curtiff.close();

warning('on','MATLAB:imagesci:tiffmexutils:libtiffWarning')
warning('on','imageio:tiffmexutils:libtiffWarning');
