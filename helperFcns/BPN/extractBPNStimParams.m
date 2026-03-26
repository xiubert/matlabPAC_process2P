function stimParams = extractBPNStimParams(triggerParams,pulse)
%write in BPN logic

%organizes 'TriggerParams' and 'Stim' into relevant params:
stimParams.trigDelay = triggerParams.stimDelay;
stimParams.Pulse = pulse.pulsename;

stimParams.BPNonsetInPulse = str2double(string(regexp(stimParams.Pulse,'(?<=_)\d+(\.\d+)?(?=sOnset)','match')));
stimParams.BPNfreqcutoff = string(regexp(stimParams.Pulse,'(?<=BPN_)\d+-\d+(?=kHz)','match'));
stimParams.BPNdBampl = str2double(string(regexp(stimParams.Pulse,'(?<=_)\d+(?=dB)','match')));
stimParams.msBPNpulseLen = str2double(string(regexp(stimParams.Pulse,'(?<=_)\d+(?=msPulse)','match')));
stimParams.msStimLen = str2double(string(regexp(stimParams.Pulse,'(?<=_)\d+(?=msTotal)','match')));


end