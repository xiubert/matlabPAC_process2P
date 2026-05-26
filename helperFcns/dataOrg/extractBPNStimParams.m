function stimParams = extractBPNStimParams(params,pulse)

stimParams.trigDelay = params.stimDelay;
stimParams.ISI = params.ISI;

nPulse = params.totalPulses;
stimParams.totalPulses = params.totalPulses;

for i = 1:nPulse
    tmp(i).Pulse = pulse(i).pulsename;
    tmp(i).sonset = str2double(string(regexp(tmp(i).Pulse,'(?<=_)\d+(\.\d+)?(?=sOnset)','match')));
    
    tmp(i).dBampl = str2double(string(regexp(tmp(i).Pulse,'(?<=_)\d+(?=dB)','match')));
    tmp(i).mspulseLen = str2double(string(regexp(tmp(i).Pulse,'(?<=_)\d+(?=msPulse)','match')));
    tmp(i).msStimLen = str2double(string(regexp(tmp(i).Pulse,'(?<=_)\d+(?=msTotal)','match')));
    tmp(i).BPNfreqcutoff = regexp(tmp(i).Pulse,'(?<=BPN_)\d+-\d+(?=kHz)','match');
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