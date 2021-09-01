function fileNames = splitTifChans(mergedChanTif,varargin)
% splitTifChans  Split ScanImage .tif with multiple channels into .tif files with single channels.
%   fileNames = splitTifChans(mergedChanTif,varargin)
%
%   Additional input arguments: 
%       '1' --> saves only tif with specified channel, replace with desired channel number
%             
%   PAC_20200213
%
%   See also readSCIMtif.m

%only for SCIM 5+
%prevents printing of irrelevant tif header errors
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');

imginfo = imfinfo(mergedChanTif);

%get header for writing SCIM compatible outputs
hInputTif = Tiff(mergedChanTif);
tagList = hInputTif.getTagNames;

if ~isempty(varargin) && cellfun(@(c) numel(c)==1 && isnumeric(c),varargin)
    img(1).chanID = varargin{cellfun(@(c) numel(c)==1 && isnumeric(c),varargin)};
    tmp = readSCIMtif(mergedChanTif,img(1).chanID);
    img(1).img = tmp;
    clear tmp
else
    img = readSCIMtif(mergedChanTif);
end

for nChan = 1:length(img)
    fileNames{nChan} = strrep(mergedChanTif,'.tif',['_chan' num2str(img(nChan).chanID) '.tif']);
    tifWrite = Tiff(fileNames{nChan},'w');
    
    for k = 1:size(img(nChan).img,3)
        
        for tagID = 1:length(tagList)
            if ~any(contains({'SubFileType','StripOffsets','YCbCrSubSampling',...
                    'NumberOfInks','MinSampleValue','MaxSampleValue'},tagList{tagID}))
                try
                    setTag(tifWrite,tagList{tagID},hInputTif.getTag(tagList{tagID}));
                catch
                end
            end
        end
        setTag(tifWrite,'BitsPerSample',imginfo(1).BitsPerSample);
        setTag(tifWrite,'SampleFormat',Tiff.SampleFormat.Int);
        setTag(tifWrite,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
        setTag(tifWrite,'Copyright',['splitTif_channel_' num2str(nChan)]);
        
        write(tifWrite,int16(img(nChan).img(:,:,k)))
        if k<size(img(nChan).img,3)
            tifWrite.writeDirectory;
        end
    end
    close(tifWrite)
end
close(hInputTif)
warning('on','MATLAB:imagesci:tiffmexutils:libtiffWarning')

