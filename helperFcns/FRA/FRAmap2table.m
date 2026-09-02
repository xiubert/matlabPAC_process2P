function T = FRAmap2table(FRA,animal,varargin)
% FRAmap2table  Long-form (ROI x freq x dB) table from an FRAmap struct.
%
%   T = FRAmap2table(FRA,animal)
%   T = FRAmap2table(FRA,animal,'roiID',ids,'treatment',"none FRAmap")
%
%   FRAmap returns a struct whose response data lives in a nDB x nFreq cell
%   array. Every other stim family in this pipeline stores per-cell responses
%   as a long-form table, which is what aggregateStimGroup, validateStimGroup
%   and groupN all consume. This converts one into the other so FRA joins the
%   generic cohort path rather than needing a parallel one.
%
%   One row per (ROI, frequency, level). Peak and significance are per
%   condition, already trial-averaged by FRAmap -- this function derives
%   nothing and re-tests nothing.
%
%   Inputs:
%     FRA    - FRAoutput struct from FRAmap. Must carry .params, i.e. be
%              built by the trial-averaged pipeline; a pre-fix file (no
%              .params) holds per-trial significance and a peak-SQUARED
%              response map and is rejected.
%     animal - animal ID, e.g. "TO0003".
%
%   Name/Value:
%     'roiID'     - ROI labels, one per cell. Default "1".."nCell" as
%                   strings, matching how the BPN/CGC tables label ROIs.
%     'treatment' - treatment label copied onto every row. Default
%                   "none FRAmap".
%     'frameRate' - frame rate stamped on every row. Default inferred from
%                   the aligned time axis.
%
%   Output table columns:
%     animal roiID treatment frameRate   identity
%     freqHz dBSPL                       the stimulus condition
%     pkDFF                              peak of the trial-averaged trace
%     sig                                LOGICAL significance for that peak
%     t_PTrel dFFavg                     onset-aligned time axis and the
%                                        trial-averaged trace (1 x nFrames)
%     nTrials                            presentations contributing
%
%   Note dFFavg is the trial-average FRAmap tested, so a cohort plot that
%   stacks per-cell traces gets exactly the traces significance was computed
%   from.
%
%   See also FRAmap, stimGroupSpec, aggregateStimGroup, tableFRAmetrics.

p = inputParser;
addRequired(p,'FRA',@isstruct)
addRequired(p,'animal',@(x) ischar(x)||isstring(x))
addParameter(p,'roiID',[],@(x) isempty(x)||isstring(x)||iscellstr(x)||isnumeric(x)) %#ok<ISCLSTR>
addParameter(p,'treatment',"none FRAmap",@(x) ischar(x)||isstring(x))
addParameter(p,'frameRate',[],@(x) isempty(x)||isscalar(x))
parse(p,FRA,animal,varargin{:});
animal    = string(p.Results.animal);
treatment = string(p.Results.treatment);
frameRate = p.Results.frameRate;

if ~isfield(FRA,'params')
    error('FRAmap2table:staleInput', ...
        ['FRAmap struct has no .params field, so it predates trial-averaged ' ...
         'significance: its sigPkDFF holds peak VALUES and its response map ' ...
         'is peak SQUARED. Regenerate with processFRA before aggregating.']);
end

nDB   = numel(FRA.dBlist);
nFreq = numel(FRA.freqList);
nCell = numel(FRA.BFuDB);

roiID = p.Results.roiID;
if isempty(roiID)
    roiID = string((1:nCell)');
else
    roiID = string(roiID(:));
    if numel(roiID)~=nCell
        error('FRAmap2table:roiCount','roiID has %d entries, FRA has %d cells', ...
            numel(roiID),nCell);
    end
end

% time axis is shared by every condition; tPTrel is nFrames x nTrials
tRef = [];
for i = 1:numel(FRA.dBFreqMap)
    if ~isempty(FRA.dBFreqMap{i}.tPTrel)
        tRef = FRA.dBFreqMap{i}.tPTrel(:,1);
        break
    end
end
if isempty(tRef)
    error('FRAmap2table:noTrials','no condition in this FRA has any trials');
end
% express relative to tone onset so animals with different absolute pulse
% times share one axis
onsetCol = FRA.params.onsetCol;
tRel = reshape(tRef - tRef(onsetCol),1,[]);
nFrames = numel(tRel);
if isempty(frameRate)
    frameRate = 1/median(diff(tRel));
end

nRow = nCell*nFreq*nDB;
[animalC,roiC,treatC] = deal(strings(nRow,1));
[freqC,dbC,pkC,frC,ntC] = deal(nan(nRow,1));
sigC = false(nRow,1);
[tC,traceC] = deal(cell(nRow,1));

r = 0;
for dbRow = 1:nDB
    for fqCol = 1:nFreq
        c = FRA.dBFreqMap{dbRow,fqCol};
        nTr = size(c.dFFptRel,3);
        if isempty(c.dFFavg)
            avg = nan(nCell,nFrames);          % missing freq/dB pair
        else
            avg = c.dFFavg;
        end
        for k = 1:nCell
            r = r+1;
            animalC(r) = animal;  roiC(r) = roiID(k);  treatC(r) = treatment;
            freqC(r) = FRA.freqList(fqCol);
            dbC(r)   = FRA.dBlist(dbRow);
            pkC(r)   = c.pkDFF(k);
            sigC(r)  = c.sigPkDFF(k);
            frC(r)   = frameRate;
            ntC(r)   = nTr;
            tC{r}    = tRel;
            traceC{r} = reshape(avg(k,:),1,[]);
        end
    end
end

T = table(animalC,roiC,treatC,frC,freqC,dbC,pkC,sigC,tC,traceC,ntC, ...
    'VariableNames',{'animal','roiID','treatment','frameRate', ...
                     'freqHz','dBSPL','pkDFF','sig','t_PTrel','dFFavg','nTrials'});

% sort so rows read cell-major, then level, then frequency
T = sortrows(T,{'animal','roiID','dBSPL','freqHz'});
end %function
