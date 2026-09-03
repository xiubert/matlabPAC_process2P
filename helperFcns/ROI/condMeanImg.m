function [meanImg,nFramesTotal] = condMeanImg(dataPath,condTifs,opts)
% condMeanImg  Frame-count-weighted mean of a condition's motion-corrected tifs.
%
%   [meanImg,nFrames] = condMeanImg(dataPath,condTifs)
%   meanImg = condMeanImg(dataPath,condTifs,'moCorrected',false)
%
%   Reproduces exactly what TIFcatROIgui's "Show mean image" toggle displays
%   -- mean(moCorrImgNonRigid,3) over a condition's whole concatenated stack
%   -- WITHOUT holding that stack in memory. Each tif is projected on its own
%   through flattenTif and the projections are combined weighted by frame
%   count, which is identical to projecting the concatenation.
%
%   Use this on the resume path, where the corrected stack would otherwise
%   have to be re-read from disk purely to take its mean. Straight after
%   motion correction, when the stack is already in the workspace,
%   mean(moCorrImgNonRigid.(cond),3) is equivalent and free.
%
%   Inputs
%     dataPath  animal data folder (holds NoRMCorred/).
%     condTifs  struct array of that condition's ORIGINAL tifs (a tifList
%               field); names are mapped to their _NoRMCorre.tif counterparts.
%
%   Name-value
%     'moCorrected'  true (default) reads NoRMCorred/<name>_NoRMCorre.tif;
%                    false reads the raw tifs in place.
%     'verbose'      Default false.
%
%   Outputs
%     meanImg       H x W double projection.
%     nFramesTotal  frames that went into it.
%
%   See also flattenTif, cellposeROIset, loadNoRMCorrNonRigidImgViaTifs

arguments
    dataPath (1,:) char
    condTifs       struct
    opts.moCorrected (1,1) logical = true
    opts.verbose     (1,1) logical = false
end

if isempty(condTifs)
    error('condMeanImg:noTifs','condTifs is empty.');
end

acc = [];
nFramesTotal = 0;
for k = 1:numel(condTifs)
    if opts.moCorrected
        f = fullfile(dataPath,'NoRMCorred',...
            strrep(condTifs(k).name,'.tif','_NoRMCorre.tif'));
    else
        f = fullfile(condTifs(k).folder,condTifs(k).name);
    end
    if ~isfile(f)
        error('condMeanImg:missingTif','Not found: %s',f);
    end

    m = flattenTif(f,'write',false);
    if isstruct(m)
        %by this stage the pipeline has already split channels, so a
        %multi-channel projection means we were pointed at the wrong file --
        %picking one arbitrarily would quietly build the mean from the
        %structural channel
        error('condMeanImg:multiChannel',...
            ['%s still holds %d channels. Motion-corrected tifs are single '...
             'channel; extract the functional channel first (splitTifChans).'],...
            f,numel(m));
    end
    n = numel(imfinfo(f));

    if isempty(acc)
        acc = zeros(size(m));
    elseif ~isequal(size(acc),size(m))
        error('condMeanImg:sizeMismatch',...
            '%s is %s but the condition started at %s.',...
            f,mat2str(size(m)),mat2str(size(acc)));
    end
    acc = acc + double(m).*n;
    nFramesTotal = nFramesTotal + n;

    if opts.verbose
        fprintf('  %s (%d frames)\n',condTifs(k).name,n);
    end
end

meanImg = acc./nFramesTotal;
end
