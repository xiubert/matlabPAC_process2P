function stimParams = extractSpontParams(params,pulse)

stimParams.trigDelay = params.stimDelay;
stimParams.ISI = params.ISI;

stimParams.totalPulses = params.totalPulses;

for i = 1:params.totalPulses
    tmp(i).Pulse = pulse(i).pulsename;
    tmp(i).msStimLen = str2double(string(regexp(tmp(i).Pulse,'(?<=_)\d+(?=msTotal)','match')));
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