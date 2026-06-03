function tifFileList = mergeFISSAgroups(tifFileList, fissaDir, fissaScaleFactor)
% mergeFISSAgroups  Merge per-ROI-count FISSA group outputs into tifFileList.
%
%   tifFileList = mergeFISSAgroups(tifFileList, fissaDir, fissaScaleFactor)
%
%   For sessions that mix 256x256 and 256x128 acquisitions, FISSA is run once
%   per ROI-count group (see FISSAviaMatlab_prePostTreatment.py), producing
%   FISSAoutput/groups.json + one g<k>/matlab.mat per group. This attaches the
%   per-tif neuropil-corrected traces to each tifFileList.stim/.map entry,
%   matched BY TIF NAME via the manifest (so ordering never matters), and
%   computes SCALEDfissaFroi. Each tif's trace row count follows its own
%   group's ROI count (e.g. 18 for 256x256 stim tifs, 11 for 256x128 spont),
%   which is exactly what stimParam2ROI's resolveROIset keys on downstream.
%
%   Inputs
%     tifFileList       struct with .stim and/or .map struct arrays; each entry
%                       must carry .name and .moCorRawFroi (nROI x nFrames).
%     fissaDir          path to NoRMCorred/FISSAoutput (contains groups.json).
%     fissaScaleFactor  neuropil scale: SCALED = moCorRawF - s*(moCorRawF - fissaF).
%
%   Used by processAnimal2P section 9 when groups.json exists; the single-group
%   (legacy / no-256x128) case uses FISSAoutput2tifFileList instead.
%
%   See also FISSAoutput2tifFileList, FISSAviaMatlab_prePostTreatment.py

manifest = jsondecode(fileread(fullfile(fissaDir,'groups.json')));
manifest = manifest(:);

% load each group's FISSA result; map moCorr-tif basename -> [groupIdx trialIdx]
grpResult = cell(numel(manifest),1);
key2gt = containers.Map('KeyType','char','ValueType','any');
for gi = 1:numel(manifest)
    Gd = load(fullfile(fissaDir, manifest(gi).matfile),'result');
    grpResult{gi} = Gd.result;
    names = cellstr(manifest(gi).tifNames);
    for ti = 1:numel(names)
        key2gt(names{ti}) = [gi ti];
    end
end

for fld = ["stim","map"]
    f = char(fld);
    if ~isfield(tifFileList,f) || isempty(tifFileList.(f))
        continue
    end
    arr = tifFileList.(f);
    for n = 1:numel(arr)
        key = strrep(arr(n).name,'.tif','_NoRMCorre.tif');   % moCorr basename
        if ~isKey(key2gt,key)
            error('mergeFISSAgroups:noGroup',...
                '%s not found in any FISSA group manifest.',arr(n).name);
        end
        gt  = key2gt(key);
        res = grpResult{gt(1)};
        ti  = gt(2);
        cells = sortByNum(fieldnames(res),'cell');           % cell0..cellN
        rows = cell(numel(cells),1);
        for c = 1:numel(cells)
            trials = res.(cells{c});
            tnames = sortByNum(fieldnames(trials),'trial');  % trial0..trialT
            tr = trials.(tnames{ti});                        % (nNeuropil+1) x nFrames
            rows{c} = tr(1,:);                               % row 1 = ROI signal
        end
        ff = vertcat(rows{:});                               % nROI x nFrames
        arr(n).fissaFroi = ff;
        arr(n).SCALEDfissaFroi = arr(n).moCorRawFroi - ...
            fissaScaleFactor .* (arr(n).moCorRawFroi - ff);
    end
    tifFileList.(f) = arr;
end
end

% ------------------------------------------------------------------------
function s = sortByNum(names, prefix)
% order {'cell0','cell10','cell2',...} numerically by the trailing integer
num = cellfun(@(x) str2double(x(numel(prefix)+1:end)), names);
[~,o] = sort(num);
s = names(o);
end
