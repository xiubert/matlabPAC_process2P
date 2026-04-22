function fileNames = splitTifChans(mergedChanTif,varargin)
% SPLITTIFCHANS  Split a multi-channel ScanImage .tif into single-channel .tif files.
%
%   fileNames = splitTifChans(mergedChanTif)
%   fileNames = splitTifChans(mergedChanTif, chanNum)
%
%   Reads a ScanImage 5+ multi-channel .tif and writes one output file per
%   channel, preserving the original ScanImage-compatible TIFF tags and
%   storing frames as int16.  Output filenames are derived from the input
%   by appending '_chan<N>' before the .tif extension.
%
%   Inputs:
%     mergedChanTif - path to a multi-channel ScanImage .tif file
%     chanNum       - (optional) scalar integer; if provided, only the
%                     specified channel is extracted and saved
%
%   Output:
%     fileNames - cell array of output file paths, one per channel written
%
%   See also readSCIMtif

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

