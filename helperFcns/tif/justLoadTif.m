%this function loads a scanimage type tif but ignores metadata
%this doesn't handle 2 channel tifs
function justImgData = justLoadTif(TifFile)
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
