function [meanImg,outFile] = flattenTif(tifFile,varargin)
%FLATTENTIF  Mean-project a ScanImage tif across frames (3rd dimension).
%
%   [meanImg,outFile] = flattenTif(tifFile)
%   [meanImg,outFile] = flattenTif(tifFile,chanNum)
%   [meanImg,outFile] = flattenTif(...,'reader','auto','outputPath',dir,...
%                                      'suffix','_mean','class','int16')
%
%   Averages every frame of a tif and writes the projection through
%   writeTifWithHeader, so the result still loads with readSCIMtif (and with
%   ScanImageTiffReader) and keeps its acquisition metadata.
%
%   Two read paths:
%     'sitiff' - Vidrio's ScanImageTiffReader mex (via hasSITiffReader).
%                Much faster than walking directories, but pulls the whole
%                stack into memory as int16.
%     'tiff'   - MATLAB's Tiff class, one directory at a time.  Slow on
%                stacks with thousands of directories, but memory stays flat
%                at one double frame.
%   'auto' (default) uses ScanImageTiffReader when it is usable AND the file
%   is at or below 'maxRAM' GB on disk, otherwise streams.  A failure inside
%   ScanImageTiffReader warns and falls back rather than erroring out.
%
%   Inputs
%     tifFile - path to a .tif (ScanImage 3.x, 5.x, or a split
%               single-channel file from splitTifChans)
%     chanNum - optional scalar; index into hChannels.channelSave, same
%               meaning as readSCIMtif's channel argument.  Omit to flatten
%               every saved channel.
%
%   Name/value pairs
%     'reader'     - 'auto' (default) | 'sitiff' | 'tiff'.  'sitiff' errors
%                    if the package is unusable and ignores the RAM guard.
%     'pkgPath'    - folder containing +ScanImageTiffReader
%                    (default /home/pac/Documents/MATLAB)
%     'maxRAM'     - GB of file-on-disk above which 'auto' streams instead of
%                    loading (default 4).  File size is a good proxy for the
%                    int16 stack size and costs nothing to check, unlike
%                    imfinfo on a 10k-directory tif.
%     'permuteXY'  - 'auto' (default) | true | false.  ScanImageTiffReader
%                    returns data in file order, i.e. transposed relative to
%                    Tiff.read/imread.  'auto' settles it by comparing
%                    directory 1 against MATLAB's own read instead of
%                    trusting a hard-coded permute.
%     'outputPath' - folder for the output (default: alongside the input)
%     'suffix'     - appended before .tif (default '_mean')
%     'class'      - 'int16' (default) or 'single', passed to writeTifWithHeader
%     'write'      - false to return the projection without writing it
%
%   Outputs
%     meanImg - height x width double projection, in MATLAB orientation
%               regardless of which reader was used.  When more than one
%               channel is flattened, a struct array with .chanID and .img,
%               matching readSCIMtif / splitTifChans conventions.
%     outFile - path written, or a cell array of paths for multiple
%               channels.  [] when 'write' is false.
%
%   Examples
%     mn    = flattenTif('mouse01_00001.tif',1);              %channel 1 mean
%     [~,f] = flattenTif('mouse01_00001_NoRMCorre.tif');       %mean of a registered stack
%     mn    = flattenTif(bigFile,'reader','tiff');             %force streaming
%
%   See also writeTifWithHeader, hasSITiffReader, readSCIMtif, splitTifChans

ip = inputParser;
ip.FunctionName = mfilename;
ip.addOptional('chanNum',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)));
ip.addParameter('reader','auto',@(s)any(strcmpi(s,{'auto','sitiff','tiff'})));
ip.addParameter('pkgPath','',@(s)ischar(s)||isstring(s));
ip.addParameter('maxRAM',4,@(x)isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('permuteXY','auto',@(x)islogical(x)||isnumeric(x)||...
    (ischar(x)||isstring(x))&&strcmpi(x,'auto'));
ip.addParameter('outputPath','',@(s)ischar(s)||isstring(s));
ip.addParameter('suffix','_mean',@(s)ischar(s)||isstring(s));
ip.addParameter('class','int16',@(s)any(strcmpi(s,{'int16','single'})));
ip.addParameter('write',true,@(x)isscalar(x)&&(islogical(x)||isnumeric(x)));
ip.parse(varargin{:});
opt = ip.Results;

tifFile = char(tifFile);
if ~isfile(tifFile)
    error([mfilename ':noFile'],'Can''t find %s',tifFile)
end

%prevents printing of irrelevant tif header errors
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
warning('off','imageio:tiffmexutils:libtiffWarning');

%--- pick a reader --------------------------------------------------------
switch lower(char(opt.reader))
    case 'tiff'
        useSIR = false;
    case 'sitiff'
        useSIR = hasSITiffReader(opt.pkgPath);
        if ~useSIR
            error([mfilename ':noSIReader'],['ScanImageTiffReader requested but '...
                'not usable (looked for +ScanImageTiffReader under %s)'],...
                char(hasSITiffReaderPath(opt.pkgPath)))
        end
    case 'auto'
        fileGB = dir(tifFile).bytes/2^30;
        useSIR = fileGB <= opt.maxRAM && hasSITiffReader(opt.pkgPath);
end

vol = [];
if useSIR
    try
        [vol,imHeader] = readWithSIR(tifFile);
    catch ME
        warning([mfilename ':sirFailed'],['ScanImageTiffReader failed (%s); '...
            'falling back to MATLAB''s Tiff class'],ME.message)
        useSIR = false;
        vol    = [];
    end
end

if useSIR
    %one cheap directory read gives us both the orientation reference and the
    %Copyright tag, which is how readSCIMtif recognises a split file.  The
    %SI header text alone can't tell us: splitTifChans copies Software
    %verbatim, so a single-channel split file still claims 2 channels there.
    [refFrame,copyTag] = peekFirstDirectory(tifFile);
    if isempty(fieldnames(imHeader))
        [~,imHeader] = readSCIMtif(tifFile,'metaOnly');   %SI 3.x: no metadata block
    end
    if ~isempty(copyTag) && ~isfield(imHeader,'tifFromSplit')
        imHeader.tifFromSplit = copyTag;
    end
    [needsTranspose,matchedRef] = sirNeedsTranspose(vol(:,:,1),refFrame,opt.permuteXY);
    if ~matchedRef
        %frame 1 matches MATLAB's own read in NEITHER orientation, so the mex
        %is not returning the pixels this file actually holds -- seen on tifs
        %written by writeMoCorTifs.  Guessing "it must be transposed" turns
        %that into a silently wrong projection (a mean of int16 data coming
        %back as +-8000), so drop to the reader that agrees with imread.
        warning([mfilename ':sirMismatch'],['ScanImageTiffReader''s frame 1 does '...
            'not match MATLAB''s in either orientation for %s; falling back to '...
            'the Tiff class rather than guessing.'],tifFile)
        useSIR = false;
        vol    = [];
    else
        nDirs = size(vol,3);
        accH  = size(vol,1);
        accW  = size(vol,2);
    end
end

if ~useSIR
    imginfo        = imfinfo(tifFile);
    nDirs          = numel(imginfo);   %actual directories, not what the header claims
    accH           = imginfo(1).Height;
    accW           = imginfo(1).Width;
    [~,imHeader]   = readSCIMtif(tifFile,'metaOnly');
    needsTranspose = false;
end

%--- how are channels laid out? -------------------------------------------
%SI 5+ interleaves channels across directories; a file already split by
%splitTifChans carries a Copyright tag (imHeader.tifFromSplit) and holds one
%channel.  SI 3.x is read sequentially, as readSCIMtif does.
if isfield(imHeader,'hChannels') && ~isfield(imHeader,'tifFromSplit')
    chanSave = imHeader.hChannels.channelSave(:)';
else
    chanSave = 1;
end
nChan = numel(chanSave);

if isfield(imHeader,'acq') && isfield(imHeader.acq,'numberOfChannelsSave') && ...
        imHeader.acq.numberOfChannelsSave > 1
    warning([mfilename ':si3MultiChan'],['SI 3.x file with %d channels saved: '...
        'directories are averaged as-is (no de-interleaving), same as readSCIMtif.'],...
        imHeader.acq.numberOfChannelsSave)
end

if isempty(opt.chanNum)
    chanIdx = 1:nChan;
else
    chanIdx = opt.chanNum;
    if chanIdx < 1 || chanIdx > nChan
        error([mfilename ':badChan'],...
            'Channel %d requested but this file holds %d channel(s)',chanIdx,nChan)
    end
end

%--- accumulate -----------------------------------------------------------
[srcDir,srcName] = fileparts(tifFile);
if isempty(char(opt.outputPath)), outDir = srcDir; else, outDir = char(opt.outputPath); end

if ~useSIR, hSrc = Tiff(tifFile,'r'); end
try
    for iCh = 1:numel(chanIdx)
        dirList = chanIdx(iCh):nChan:nDirs;
        acc = zeros(accH,accW);
        if useSIR
            %frame-at-a-time cast keeps the peak at the int16 stack plus one
            %double frame, instead of double()-ing the whole substack
            for k = dirList
                acc = acc + double(vol(:,:,k));
            end
        else
            for k = dirList
                hSrc.setDirectory(k);
                acc = acc + double(hSrc.read());
            end
        end
        proj = acc./numel(dirList);

        %transposing the projection == transposing every frame, and is free
        if needsTranspose, proj = proj.'; end

        if opt.write
            if nChan > 1
                nameOut = sprintf('%s_chan%d%s.tif',srcName,...
                    chanSave(chanIdx(iCh)),char(opt.suffix));
            else
                nameOut = [srcName char(opt.suffix) '.tif'];
            end
            fOut = fullfile(outDir,nameOut);
            writeTifWithHeader(proj,fOut,tifFile,...
                'class',opt.class,...
                'copyright',sprintf('flattenTif_mean_%dframes_chan%d',...
                                    numel(dirList),chanSave(chanIdx(iCh))));
        else
            fOut = [];
        end

        if numel(chanIdx) > 1
            meanImg(iCh).chanID = chanSave(chanIdx(iCh)); %#ok<AGROW>
            meanImg(iCh).img    = proj;                   %#ok<AGROW>
            outFile{iCh}        = fOut;                   %#ok<AGROW>
        else
            meanImg = proj;
            outFile = fOut;
        end
    end
    if ~useSIR, hSrc.close(); end
catch ME
    if ~useSIR, hSrc.close(); end
    warning('on','MATLAB:imagesci:tiffmexutils:libtiffWarning');
    warning('on','imageio:tiffmexutils:libtiffWarning');
    rethrow(ME)
end

warning('on','MATLAB:imagesci:tiffmexutils:libtiffWarning');
warning('on','imageio:tiffmexutils:libtiffWarning');

end
%--------------------------------------------------------------------------
function [vol,imHeader] = readWithSIR(tifFile)
%whole-stack read through Vidrio's mex, plus the SI header text it carries
hSIR = ScanImageTiffReader.ScanImageTiffReader(tifFile);
try
    vol  = hSIR.data();
    meta = char(hSIR.metadata());
catch ME
    hSIR.close();
    rethrow(ME)
end
hSIR.close();

if isempty(vol)
    error('flattenTif:emptySIR','ScanImageTiffReader returned no data')
end
imHeader = parseSIheaderText(meta);
end
%--------------------------------------------------------------------------
function imHeader = parseSIheaderText(txt)
%evaluate the "SI.xxx = yyy" lines of a ScanImage 5+ header blob, same trick
%readSCIMtif uses on the Software tag.  The trailing RoiGroups JSON has no
%'SI.' prefix, so it drops out on its own.
imHeader = struct();
if isempty(txt), return, end
lines = strsplit(char(txt),'\n');
for k = 1:numel(lines)
    ln = strtrim(lines{k});
    if isempty(ln) || ~strncmp(ln,'SI.',3), continue, end
    try
        evalc(ln);
    catch
    end
end
if exist('SI','var'), imHeader = SI; end
end
%--------------------------------------------------------------------------
function [refFrame,copyTag] = peekFirstDirectory(tifFile)
%read directory 1 with MATLAB's Tiff class: gives an orientation reference
%and the Copyright tag without walking the whole file
hT = Tiff(tifFile,'r');
try
    refFrame = hT.read();
    try
        copyTag = hT.getTag('Copyright');
    catch
        copyTag = '';
    end
catch ME
    hT.close();
    rethrow(ME)
end
hT.close();
end
%--------------------------------------------------------------------------
function [tf,matched] = sirNeedsTranspose(sirFrame,refFrame,mode)
%ScanImageTiffReader hands back data in file (row-major) order, so it is
%normally the transpose of Tiff.read.  Rather than hard-coding
%permute(vol,[2 1 3]), check it against directory 1.
%
%matched says whether the two readers agree at all.  When they do not, the
%caller falls back to the Tiff class: a mismatch means the mex is handing
%back different PIXELS, not merely a different orientation, and no permute
%can rescue that.
if ~(ischar(mode) || isstring(mode))
    tf = logical(mode);
    matched = true;
    return
end
a = double(refFrame);
b = double(sirFrame);
matched = true;
if isequal(size(a),size(b)) && isequal(a,b)
    tf = false;
elseif isequal(size(a),size(b.')) && isequal(a,b.')
    tf = true;
else
    tf = true;
    matched = false;
end
end
%--------------------------------------------------------------------------
function p = hasSITiffReaderPath(pkgParent)
%only for the error message: mirrors hasSITiffReader's default
if isempty(char(pkgParent))
    p = '/home/pac/Documents/MATLAB';
else
    p = char(pkgParent);
end
end