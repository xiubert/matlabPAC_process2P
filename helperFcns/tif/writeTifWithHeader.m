function outFile = writeTifWithHeader(imgData,outFile,srcTif,varargin)
%WRITETIFWITHHEADER  Write a ScanImage-readable .tif, inheriting a source header.
%
%   outFile = writeTifWithHeader(imgData,outFile,srcTif)
%   outFile = writeTifWithHeader(...,'copyright',str,'class','int16')
%
%   Copies every writable tag from SRCTIF (including the ScanImage header
%   carried in the Software / ImageDescription tags) onto a new file, then
%   writes IMGDATA frame by frame.  Same tag-copying trick as
%   splitTifChans and writeMoCorTifs, factored out so that derived files
%   (projections, filtered stacks, registered stacks) stay loadable by
%   readSCIMtif / scim_openTif.
%
%   Inputs
%     imgData - height x width x nFrames numeric array (single channel)
%     outFile - path of the .tif to create
%     srcTif  - path to the tif whose header is inherited, or an already
%               open Tiff object (note: its current directory is moved to 1)
%
%   Name/value pairs
%     'copyright'     - text for the Copyright tag.
%                       default: 'writeTifWithHeader'
%                       readSCIMtif reads a non-empty Copyright as
%                       imHeader.tifFromSplit, which sends the file down the
%                       "one channel, read every directory" path.  That is
%                       what you want for anything derived from a
%                       multi-channel acquisition.
%     'class'         - 'int16' (default, ScanImage-native) or 'single'.
%                       int16 rounds and saturates, which is safe for means
%                       of int16 data; use 'single' if the source is not
%                       signed 16-bit or you need sub-integer precision
%                       (readSCIMtif handles floats, scim_openTif may not).
%     'fixFrameCount' - true (default) rewrites the frame-count fields in
%                       the inherited SI header to size(imgData,3).  Without
%                       it, an inherited logFramesPerFile of 1000 makes
%                       readSCIMtif ask for directory 1000 of a one-frame
%                       projection and error out.
%     'descriptions'  - per-frame ImageDescription text: one char/string for
%                       all frames, or a cell array with one entry per frame.
%                       Feed it ScanImageTiffReader's descriptions() output
%                       when you are rewriting a stack and want the real
%                       per-frame timestamps to survive.
%     'renumberFrames'- true (default).  When 'descriptions' is not given,
%                       directory 1's description is inherited by every frame
%                       and frameNumbers / frameNumberAcquisition are
%                       renumbered 1..nFrames, so frames 2..N stop claiming
%                       to be frame 1 (visible via descriptions()).
%     'bigTiff'       - true opens in 'w8' mode for >4 GB output (default false)
%
%   ScanImageTiffReader is read-only, so writing stays on MATLAB's Tiff
%   class.  The Software tag is preserved verbatim apart from the frame-count
%   patch, so the output still parses through both readSCIMtif and
%   ScanImageTiffReader.metadata().
%
%   Example
%     writeTifWithHeader(meanImg,'mouse01_00001_mean.tif', ...
%         'mouse01_00001.tif','copyright','flattenTif_mean');
%
%   See also readSCIMtif, splitTifChans, flattenTif

ip = inputParser;
ip.FunctionName = mfilename;
ip.addParameter('copyright',mfilename,@(s)ischar(s)||isstring(s));
ip.addParameter('class','int16',@(s)any(strcmpi(s,{'int16','single'})));
ip.addParameter('fixFrameCount',true,@(x)isscalar(x)&&(islogical(x)||isnumeric(x)));
ip.addParameter('renumberFrames',true,@(x)isscalar(x)&&(islogical(x)||isnumeric(x)));
ip.addParameter('descriptions',{},@(c)iscell(c)||ischar(c)||isstring(c));
ip.addParameter('bigTiff',false,@(x)isscalar(x)&&(islogical(x)||isnumeric(x)));
ip.parse(varargin{:});
opt = ip.Results;

if ~isnumeric(imgData) || isempty(imgData)
    error([mfilename ':badData'],'imgData must be a non-empty numeric array')
end
if ndims(imgData) > 3
    error([mfilename ':badData'],'imgData must be 2-D or 3-D (h x w x nFrames)')
end
outFile  = char(outFile);
nFrames  = size(imgData,3);

%prevents printing of irrelevant tif header errors
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
warning('off','imageio:tiffmexutils:libtiffWarning');

%--- source header --------------------------------------------------------
if isa(srcTif,'Tiff')
    hSrc     = srcTif;
    closeSrc = false;
else
    hSrc     = Tiff(char(srcTif),'r');
    closeSrc = true;
end
hSrc.setDirectory(1);

%tags libtiff maintains itself, or that we set explicitly below
skipTags = {'SubFileType','StripOffsets','StripByteCounts','TileOffsets',...
    'TileByteCounts','YCbCrSubSampling','NumberOfInks','MinSampleValue',...
    'MaxSampleValue','ImageLength','ImageWidth','BitsPerSample',...
    'SampleFormat','SamplesPerPixel','PlanarConfiguration','Photometric',...
    'Compression','PageNumber','Copyright','Software','ImageDescription'};

tagList = hSrc.getTagNames;
srcTags = struct();
for tagID = 1:numel(tagList)
    if any(strcmp(tagList{tagID},skipTags)), continue, end
    try
        srcTags.(tagList{tagID}) = hSrc.getTag(tagList{tagID});
    catch
        %tag simply isn't present in this file
    end
end
copyNames = fieldnames(srcTags);

%the SI header itself: SI 5+ puts it in Software, SI 3.x in ImageDescription
hdrSoftware = getTagSafe(hSrc,'Software');
hdrDescr    = getTagSafe(hSrc,'ImageDescription');
if closeSrc, hSrc.close(); end

if opt.fixFrameCount
    if ~isempty(hdrSoftware)                        %ScanImage 5+
        hdrSoftware = setHdrField(hdrSoftware,'SI.hScan2D.logFramesPerFile',nFrames);
        hdrSoftware = setHdrField(hdrSoftware,'SI.hStackManager.framesPerSlice',nFrames);
        hdrSoftware = setHdrField(hdrSoftware,'SI.hFastZ.numFramesPerVolume',nFrames);
    end
    if contains(hdrDescr,'state.acq')               %ScanImage 3.x
        hdrDescr = setHdrField(hdrDescr,'state.acq.numberOfFrames',nFrames);
    end
end

%--- per-frame ImageDescription -------------------------------------------
%SI 5+ puts frameNumbers / timestamps in each directory's ImageDescription.
%Copying directory 1's onto every frame (what splitTifChans and
%writeMoCorTifs do) leaves frames 2..N claiming to be frame 1, which shows
%up as soon as you read them back with ScanImageTiffReader.descriptions().
%Pass 'descriptions' (e.g. straight from descriptions() on the source) to
%write real per-frame text, or let renumberFrames fix up the frame indices.
if ~isempty(opt.descriptions)
    d = opt.descriptions;
    if ischar(d) || isstring(d), d = {char(d)}; end
    if numel(d) == 1
        frameDescr = repmat({char(d{1})},1,nFrames);
    elseif numel(d) == nFrames
        frameDescr = cellfun(@char,reshape(d,1,[]),'UniformOutput',false);
    else
        error([mfilename ':badDescr'],...
            '''descriptions'' needs 1 entry or one per frame (%d), got %d',...
            nFrames,numel(d))
    end
else
    frameDescr = repmat({hdrDescr},1,nFrames);
    if opt.renumberFrames && ~isempty(hdrDescr)
        for k = 1:nFrames
            frameDescr{k} = setHdrField(frameDescr{k},'frameNumbers',k);
            frameDescr{k} = setHdrField(frameDescr{k},'frameNumberAcquisition',k);
        end
    end
end

%--- output class ---------------------------------------------------------
switch lower(opt.class)
    case 'int16'
        imgOut  = int16(imgData);                   %rounds and saturates
        nBits   = 16;
        sampFmt = Tiff.SampleFormat.Int;
    case 'single'
        imgOut  = single(imgData);
        nBits   = 32;
        sampFmt = Tiff.SampleFormat.IEEEFP;
end

%--- write ----------------------------------------------------------------
outDir = fileparts(outFile);
if ~isempty(outDir) && ~isfolder(outDir), mkdir(outDir); end

if opt.bigTiff, mode = 'w8'; else, mode = 'w'; end
hOut = Tiff(outFile,mode);
try
    for k = 1:nFrames
        %tags have to be re-set for every directory.
        %structural tags first and in this order: libtiff refuses
        %SamplesPerPixel until Photometric exists, and PlanarConfiguration
        %until SamplesPerPixel exists.
        hOut.setTag('Photometric',Tiff.Photometric.MinIsBlack);
        hOut.setTag('ImageLength',size(imgOut,1));
        hOut.setTag('ImageWidth',size(imgOut,2));
        hOut.setTag('BitsPerSample',nBits);
        hOut.setTag('SampleFormat',sampFmt);
        hOut.setTag('SamplesPerPixel',1);
        hOut.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
        hOut.setTag('Compression',Tiff.Compression.None);

        %then the inherited metadata
        for tagID = 1:numel(copyNames)
            try
                hOut.setTag(copyNames{tagID},srcTags.(copyNames{tagID}));
            catch
            end
        end
        if ~isempty(hdrSoftware), hOut.setTag('Software',hdrSoftware); end
        if ~isempty(frameDescr{k}), hOut.setTag('ImageDescription',frameDescr{k}); end
        hOut.setTag('Copyright',char(opt.copyright));

        hOut.write(imgOut(:,:,k));
        if k < nFrames
            hOut.writeDirectory;
        end
    end
    hOut.close();
catch ME
    hOut.close();
    warning('on','MATLAB:imagesci:tiffmexutils:libtiffWarning');
    warning('on','imageio:tiffmexutils:libtiffWarning');
    rethrow(ME)
end

warning('on','MATLAB:imagesci:tiffmexutils:libtiffWarning');
warning('on','imageio:tiffmexutils:libtiffWarning');

end
%--------------------------------------------------------------------------
function v = getTagSafe(hTif,tagName)
%returns '' instead of throwing when the tag is absent
try
    v = hTif.getTag(tagName);
catch
    v = '';
end
end
%--------------------------------------------------------------------------
function s = setHdrField(s,name,value)
%replace "name = ..." inside a ScanImage header string, if the field exists.
%the '\s*=' guard stops 'numberOfFrames' from matching 'numberOfFramesFoo'.
pat = [regexptranslate('escape',name) '\s*=[^\r\n]*'];
if ~isempty(regexp(s,pat,'once'))
    s = regexprep(s,pat,[name ' = ' num2str(value)]);
end
end