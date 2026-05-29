function stimParams = extractBPNStimParams(params,pulse)
% extractBPNStimParams  Parse band-pass-noise (BPN) stimulus parameters
%                       from an N-pulse pulseset.
%
%   stimParams = extractBPNStimParams(params, pulse)
%
%   Parses each pulse(i).pulsename string in a BPN pulseset and packs the
%   per-pulse parameters into stimParams as N-by-1 cell-array fields.
%   Scalar (set-level) parameters are stored directly on stimParams.
%
%   Inputs:
%     params - struct with fields:
%                stimDelay   - stimulus delay (s)
%                ISI         - inter-stimulus interval (s)
%                totalPulses - number of pulses in the set (N)
%     pulse  - N-element struct array with field:
%                pulsename   - pulse name string encoding params
%                              (e.g. 'BPN_5-10kHz_..._50dB_..._100msPulse_..._500msTotal_..._1sOnset')
%
%   Output fields:
%     Scalar:
%       trigDelay        - stimulus delay (s), copied from params.stimDelay
%       ISI              - inter-stimulus interval (s)
%       totalPulses      - number of pulses (N)
%     Per-pulse (N-by-1 cell):
%       Pulse            - full pulse name string
%       BPNsOnset        - tone onset time (s) from '_#sOnset' (supports decimals)
%       BPNdBAmpl        - amplitude (dB SPL) from '_##dB'
%       BPNmsPulseLen    - pulse duration (ms) from '_###msPulse'
%       BPNmsStimLen     - total stim duration (ms) from '_###msTotal'
%       BPNkHzFreqCutoff - frequency cutoffs as 'lo-hi' string (kHz),
%                          parsed from 'BPN_#-##kHz'
%
%   Naming convention:
%     Per-pulse stim-param fields follow <StimFamily><unit><Quantity> in
%     PascalCase (e.g. BPNsOnset, BPNkHzFreqCutoff). The unit is included
%     when the value's unit is not self-evident from the quantity name.
%     This mirrors the PT/DRC convention in extractStimParams.m and the
%     Spont convention in extractSpontParams.m.

stimParams.trigDelay = params.stimDelay;
stimParams.ISI = params.ISI;

nPulse = params.totalPulses;
stimParams.totalPulses = params.totalPulses;

for i = 1:nPulse
    tmp(i).Pulse = pulse(i).pulsename;
    tmp(i).BPNsOnset = str2double(string(regexp(tmp(i).Pulse,'(?<=_)\d+(\.\d+)?(?=sOnset)','match')));

    tmp(i).BPNdBAmpl = str2double(string(regexp(tmp(i).Pulse,'(?<=_)\d+(?=dB)','match')));
    tmp(i).BPNmsPulseLen = str2double(string(regexp(tmp(i).Pulse,'(?<=_)\d+(?=msPulse)','match')));
    tmp(i).BPNmsStimLen = str2double(string(regexp(tmp(i).Pulse,'(?<=_)\d+(?=msTotal)','match')));
    tmp(i).BPNkHzFreqCutoff = regexp(tmp(i).Pulse,'(?<=BPN_)\d+-\d+(?=kHz)','match');
end
% convert 1-by-N structure array into a single struct where each field contains a column cell array
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