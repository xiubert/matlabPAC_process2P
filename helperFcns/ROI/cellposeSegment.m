function [labelImg,params] = cellposeSegment(img,opts)
% cellposeSegment  Segment one 2-D image with Cellpose, return a label image.
%
%   [labelImg,params] = cellposeSegment(img)
%   [labelImg,params] = cellposeSegment(img,'name','TO0003_all',...)
%
%   Runs the Cellpose CLI on a single frame -- normally the mean projection of
%   a motion-corrected stack -- and hands back the label image it produces
%   (0 = background, one positive integer per cell), ready for
%   labelImg2moCorROI. This is the headless stand-in for drawing ROIs in
%   TIFcatROIgui.
%
%   Inputs
%     img   H x W numeric image, OR a path to an existing .tif to segment.
%
%   Name-value -- where it runs
%     'backend'   'podman-exec' (default) | 'local' | 'none'
%                 podman-exec : exec into an already-running container. The
%                     image is staged into a directory visible from both
%                     sides (see stagingHost/stagingContainer) because the
%                     animal data drive is not bind-mounted into it.
%                 local       : `cellpose` on the host PATH; no staging, the
%                     image is written next to outDir.
%                 none        : do not run anything; read a
%                     <name>_cp_masks.tif that already exists. For re-reading
%                     a previous run, or wiring in another segmenter.
%     'container'     container name for podman-exec. Default 'cellpose'.
%     'containerUser' user to exec as. Default 'cellpose_cc'.
%     'podman'        podman executable. Default 'podman'.
%     'stagingHost'      host-side exchange dir.
%                        Default '/media/DATA/Chris/cellpose2D/mpac_stage'.
%     'stagingContainer' the same dir as seen inside the container.
%                        Default '/data/mpac_stage'.
%     'envPrefix'     command prefix used for every shell-out. Default
%                     'env -u LD_LIBRARY_PATH -u LD_PRELOAD ' -- MATLAB
%                     points LD_LIBRARY_PATH at its own bundled libs, which
%                     breaks system binaries like podman in confusing ways.
%
%   Name-value -- what it runs
%     'model'          --pretrained_model value, a CONTAINER path for
%                      podman-exec. Default
%                      '/data/cellpose_cc/.cellpose/models/cpsam'.
%     'useGpu'         Default true.
%     'normPercentile' [lo hi] for --norm_percentile. Default [1 99].
%                      Pass [] to omit (Cellpose's own default normalisation).
%     'diameter'       --diameter. Default [] = omit (cpsam sizes itself).
%     'flowThreshold'  --flow_threshold. Default [] = omit.
%     'cellprobThreshold' --cellprob_threshold. Default [] = omit.
%     'minSize'        --min_size. Default [] = omit. Prefer minAreaPx in
%                      labelImg2moCorROI so the cut is visible in the ROI QC.
%     'excludeOnEdges' --exclude_on_edges. Default false; prefer
%                      edgeMarginPx in labelImg2moCorROI, which also protects
%                      the neuropil annulus.
%     'extraArgs'      char appended verbatim to the cellpose command.
%
%   Name-value -- files
%     'name'      basename for the staged image. Default: derived from img if
%                 it is a path, else 'cellposeSegment_<timestamp>'.
%     'srcTif'    a tif whose ScanImage header the staged image inherits (via
%                 writeTifWithHeader), so the staged mean stays loadable by
%                 readSCIMtif. Default '' = plain tif.
%     'imgClass'  'int16' (default) or 'single', passed to writeTifWithHeader.
%     'keepFiles' keep the staged image and mask after reading. Default false.
%     'timeout'   seconds to allow the cellpose call. Default 1800.
%
%   Outputs
%     labelImg  H x W double label image.
%     params    provenance, to be saved alongside the ROIs: .backend .command
%               .cellposeVersion .model (+ .modelBytes/.modelMtime when
%               readable) .normPercentile .diameter .dilate-relevant flags
%               .runSeconds .when .imgSize .nLabels
%
%   Cellpose writes <name>_cp_masks.tif as the container user, so on this
%   setup the mask comes back owned by a subordinate uid: readable by the
%   host user but not deletable by it. Cleanup therefore goes back through the
%   container. Nothing downstream depends on the staged files.
%
%   See also labelImg2moCorROI, cellposeROIset, flattenTif, writeTifWithHeader

arguments
    img
    opts.backend           (1,:) char {mustBeMember(opts.backend,{'podman-exec','local','none'})} = 'podman-exec'
    opts.container         (1,:) char = 'cellpose'
    opts.containerUser     (1,:) char = 'cellpose_cc'
    opts.podman            (1,:) char = 'podman'
    opts.stagingHost       (1,:) char = '/media/DATA/Chris/cellpose2D/mpac_stage'
    opts.stagingContainer  (1,:) char = '/data/mpac_stage'
    opts.envPrefix         (1,:) char = 'env -u LD_LIBRARY_PATH -u LD_PRELOAD '
    opts.model             (1,:) char = '/data/cellpose_cc/.cellpose/models/cpsam'
    opts.useGpu            (1,1) logical = true
    opts.normPercentile          double = [1 99]
    opts.diameter                double = []
    opts.flowThreshold           double = []
    opts.cellprobThreshold       double = []
    opts.minSize                 double = []
    opts.excludeOnEdges    (1,1) logical = false
    opts.extraArgs         (1,:) char = ''
    opts.name              (1,:) char = ''
    opts.srcTif            (1,:) char = ''
    opts.imgClass          (1,:) char {mustBeMember(opts.imgClass,{'int16','single'})} = 'int16'
    opts.keepFiles         (1,1) logical = false
    opts.timeout           (1,1) double {mustBePositive} = 1800
end

tRun = tic;

%--- where does the image live, and under what name? ----------------------
imgIsPath = (ischar(img) || isstring(img)) && isfile(char(img));
if imgIsPath
    [~,defName] = fileparts(char(img));
else
    defName = sprintf('cellposeSegment_%s',datestr(now,'yyyymmdd_HHMMSSFFF')); %#ok<TNOW1,DATST>
end
name = opts.name;
if isempty(name), name = defName; end
name = matlab.lang.makeValidName(name);

switch opts.backend
    case 'podman-exec'
        workHost = opts.stagingHost;
        workCont = opts.stagingContainer;
    otherwise
        %local / none: work in place, next to the image when we were given one
        if imgIsPath
            workHost = fileparts(char(img));
            if isempty(workHost), workHost = pwd; end
        else
            workHost = tempdir;
        end
        workCont = workHost;
end
if ~isfolder(workHost)
    error('cellposeSegment:noWorkDir',...
        ['Working directory does not exist: %s\n'...
         'For the podman-exec backend this must be a directory writable by '...
         'BOTH the host user and the container user (see etc/headless_plan.md).'],...
        workHost);
end

imgHost  = fullfile(workHost,[name '.tif']);
imgCont  = [workCont '/' name '.tif'];
maskHost = fullfile(workHost,[name '_cp_masks.tif']);

%--- stage the image ------------------------------------------------------
wroteImg = false;
if strcmp(opts.backend,'none')
    %nothing to stage; the mask is expected to be there already
elseif imgIsPath && strcmp(char(img),imgHost)
    %already in place
else
    if imgIsPath
        srcImg = readImage(char(img));
    else
        srcImg = img;
    end
    if ~ismatrix(srcImg)
        error('cellposeSegment:not2D',...
            'img must be a single 2-D frame (got %s). Project the stack first.',...
            mat2str(size(srcImg)));
    end
    writeStagedImage(double(srcImg),imgHost,opts);
    wroteImg = true;
end

%--- run ------------------------------------------------------------------
cmd = '';
switch opts.backend
    case 'none'
        if ~isfile(maskHost)
            error('cellposeSegment:noMask',...
                'backend ''none'' but %s does not exist.',maskHost);
        end
    otherwise
        cpArgs = buildCellposeArgs(imgCont,opts);
        switch opts.backend
            case 'podman-exec'
                cmd = sprintf('%s%s exec -u %s %s sh -c "%s"',...
                    opts.envPrefix,opts.podman,opts.containerUser,...
                    opts.container,strrep(cpArgs,'"','\"'));
            case 'local'
                cmd = [opts.envPrefix cpArgs];
        end
        cmd = sprintf('timeout %d %s',round(opts.timeout),cmd);

        [st,out] = system(cmd);
        if st ~= 0
            cleanup(imgHost,maskHost,wroteImg,false,opts);
            if st == 124
                error('cellposeSegment:timeout',...
                    'Cellpose timed out after %g s.\ncommand: %s',opts.timeout,cmd);
            end
            error('cellposeSegment:runFailed',...
                'Cellpose failed (exit %d).\ncommand: %s\noutput:\n%s',st,cmd,out);
        end
        if ~isfile(maskHost)
            cleanup(imgHost,maskHost,wroteImg,false,opts);
            error('cellposeSegment:noMaskWritten',...
                ['Cellpose reported success but %s was not written.\n'...
                 'command: %s\noutput:\n%s'],maskHost,cmd,out);
        end
end

%--- read back ------------------------------------------------------------
labelImg = double(imread(maskHost));
if ~ismatrix(labelImg)
    error('cellposeSegment:badMask','%s is not a 2-D label image.',maskHost);
end

params = struct();
params.backend        = opts.backend;
params.command        = cmd;
params.model          = opts.model;
params.normPercentile = opts.normPercentile;
params.diameter       = opts.diameter;
params.useGpu         = opts.useGpu;
params.excludeOnEdges = opts.excludeOnEdges;
params.extraArgs      = opts.extraArgs;
params.imgSize        = size(labelImg);
params.nLabels        = numel(setdiff(unique(labelImg(:)),0));
params.runSeconds     = toc(tRun);
params.when           = datetime('now');
[params.cellposeVersion,params.modelBytes,params.modelMtime] = probeVersion(opts);

cleanup(imgHost,maskHost,wroteImg,true,opts);
end

% ------------------------------------------------------------------------
function s = buildCellposeArgs(imgCont,opts)
s = sprintf('cellpose --image_path %s --pretrained_model %s --save_tif --no_npy --verbose',...
    imgCont,opts.model);
if opts.useGpu,         s = [s ' --use_gpu']; end
if ~isempty(opts.normPercentile)
    if numel(opts.normPercentile) ~= 2
        error('cellposeSegment:badNorm','normPercentile must be [lo hi] or [].');
    end
    s = [s sprintf(' --norm_percentile %g %g',opts.normPercentile(1),opts.normPercentile(2))];
end
if ~isempty(opts.diameter),          s = [s sprintf(' --diameter %g',opts.diameter)]; end
if ~isempty(opts.flowThreshold),     s = [s sprintf(' --flow_threshold %g',opts.flowThreshold)]; end
if ~isempty(opts.cellprobThreshold), s = [s sprintf(' --cellprob_threshold %g',opts.cellprobThreshold)]; end
if ~isempty(opts.minSize),           s = [s sprintf(' --min_size %d',round(opts.minSize))]; end
if opts.excludeOnEdges,              s = [s ' --exclude_on_edges']; end
if ~isempty(opts.extraArgs),         s = [s ' ' opts.extraArgs]; end
end
% ------------------------------------------------------------------------
function writeStagedImage(imgD,imgHost,opts)
%int16 is what the rest of the pipeline writes and what this setup was
%validated against; Cellpose normalises by percentiles anyway, so absolute
%scale does not matter -- but silent saturation would.
if strcmpi(opts.imgClass,'int16')
    lim = double(intmax('int16'));
    if max(imgD(:)) > lim || min(imgD(:)) < double(intmin('int16'))
        warning('cellposeSegment:rescaled',...
            ['Image range [%g %g] exceeds int16; rescaling to 0..%d for '...
             'segmentation. Cellpose normalises internally, so this does not '...
             'change the result, but pass ''imgClass'',''single'' to avoid it.'],...
            min(imgD(:)),max(imgD(:)),lim);
        imgD = (imgD - min(imgD(:)));
        if max(imgD(:)) > 0
            imgD = imgD./max(imgD(:)).*lim;
        end
    end
end

if ~isempty(opts.srcTif) && isfile(opts.srcTif)
    writeTifWithHeader(imgD,imgHost,opts.srcTif,...
        'class',opts.imgClass,'copyright','cellposeSegment_input');
else
    %IMWRITE refuses int16, so both classes go through the Tiff class
    if strcmpi(opts.imgClass,'int16')
        px = int16(imgD); nBits = 16; sampFmt = Tiff.SampleFormat.Int;
    else
        px = single(imgD); nBits = 32; sampFmt = Tiff.SampleFormat.IEEEFP;
    end
    t = Tiff(imgHost,'w');
    try
        t.setTag('Photometric',Tiff.Photometric.MinIsBlack);
        t.setTag('ImageLength',size(px,1));
        t.setTag('ImageWidth',size(px,2));
        t.setTag('BitsPerSample',nBits);
        t.setTag('SampleFormat',sampFmt);
        t.setTag('SamplesPerPixel',1);
        t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
        t.setTag('Compression',Tiff.Compression.None);
        t.write(px);
        t.close();
    catch ME
        t.close();
        rethrow(ME)
    end
end
%the container user has to be able to read it
try, fileattrib(imgHost,'+w +r','a'); catch, end %#ok<CTCH,NOCOM>
end
% ------------------------------------------------------------------------
function img = readImage(p)
try
    [im,~] = readSCIMtif(p);
    if isstruct(im), im = im(1).img; end
    img = double(im(:,:,1));
catch
    img = double(imread(p));
end
end
% ------------------------------------------------------------------------
function [ver,mBytes,mMtime] = probeVersion(opts)
%best-effort provenance: never let a failed probe break a good segmentation
ver = ''; mBytes = NaN; mMtime = '';
try
    switch opts.backend
        case 'podman-exec'
            base = sprintf('%s%s exec -u %s %s sh -c ',...
                opts.envPrefix,opts.podman,opts.containerUser,opts.container);
            [s1,o1] = system([base '"cellpose --version"']);
            if s1==0, ver = strtrim(o1); end
            [s2,o2] = system([base sprintf('"stat -c ''%%s %%Y'' %s"',opts.model)]);
            if s2==0
                v = sscanf(strtrim(o2),'%f %f');
                if numel(v)==2
                    mBytes = v(1);
                    mMtime = char(datetime(v(2),'ConvertFrom','posixtime'));
                end
            end
        case 'local'
            [s1,o1] = system([opts.envPrefix 'cellpose --version']);
            if s1==0, ver = strtrim(o1); end
            d = dir(opts.model);
            if ~isempty(d)
                mBytes = d.bytes;
                mMtime = d.date;
            end
    end
catch
end
end
% ------------------------------------------------------------------------
function cleanup(imgHost,maskHost,wroteImg,removeMask,opts)
if opts.keepFiles, return, end
if wroteImg
    try, delete(imgHost); catch, end %#ok<CTCH,NOCOM>
end
if removeMask && isfile(maskHost) && ~strcmp(opts.backend,'none')
    %written by the container user, so the host user usually cannot unlink it
    ok = false;
    try
        delete(maskHost);
        ok = ~isfile(maskHost);
    catch
    end
    if ~ok && strcmp(opts.backend,'podman-exec')
        try
            system(sprintf('%s%s exec -u %s %s sh -c "rm -f %s"',...
                opts.envPrefix,opts.podman,opts.containerUser,opts.container,...
                [opts.stagingContainer '/' erase(maskHost,[fileparts(maskHost) filesep])]));
        catch
        end
    end
end
end
