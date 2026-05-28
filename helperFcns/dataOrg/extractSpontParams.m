function stimParams = extractSpontParams(params,pulse)
% extractSpontParams  Parse spontaneous-activity stimulus parameters from
%                     an N-pulse pulseset.
%
%   stimParams = extractSpontParams(params, pulse)
%
%   Parses each pulse(i).pulsename in a Spont pulseset and packs per-pulse
%   parameters into stimParams as N-by-1 cell-array fields. Scalar
%   (set-level) parameters are stored directly on stimParams.
%
%   Inputs:
%     params - struct with fields:
%                stimDelay   - stimulus delay (s)
%                ISI         - inter-stimulus interval (s)
%                totalPulses - number of pulses in the set (N)
%     pulse  - N-element struct array with field:
%                pulsename   - pulse name string encoding params
%                              (e.g. '..._500msTotal_...')
%
%   Output fields:
%     Scalar:
%       trigDelay       - stimulus delay (s)
%       ISI             - inter-stimulus interval (s)
%       totalPulses     - number of pulses (N)
%     Per-pulse (N-by-1 cell):
%       Pulse           - full pulse name string
%       SpontMsStimLen  - total stim duration (ms) from '_###msTotal'
%
%   Naming convention:
%     Per-pulse stim-param fields follow <StimFamily><unit><Quantity> in
%     PascalCase (e.g. SpontMsStimLen). Mirrors the BPN convention in
%     extractBPNStimParams.m and the PT/DRC convention in extractStimParams.m.

stimParams.trigDelay = params.stimDelay;
stimParams.ISI = params.ISI;

stimParams.totalPulses = params.totalPulses;

for i = 1:params.totalPulses
    tmp(i).Pulse = pulse(i).pulsename;
    tmp(i).SpontMsStimLen = str2double(string(regexp(tmp(i).Pulse,'(?<=_)\d+(?=msTotal)','match')));
end

fn = fieldnames(tmp);
vals = cell(size(fn));
for k = 1:numel(fn)
    vals{k} = {tmp.(fn{k})}.';
end
tmp_struct = cell2struct(vals, fn);

for k = 1:numel(fn)
    stimParams.(fn{k}) = tmp_struct.(fn{k});
end

end