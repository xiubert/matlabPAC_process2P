function stimParams = extractContrastChangeParams(triggerParams,pulse)

%organizes 'TriggerParams' and 'Stim' into relevant params:
stimParams.trigDelay = triggerParams.stimDelay;
pre = regexp(pulse.pulsename,'_\d$','split');
stimParams.rawPulse = pulse.pulsename; %PulseID
stimParams.Pulse = pre{1};
clear pre

stimParams.dBrange = strrep(regexp(stimParams.Pulse,'\d{2}-\d{2}dB','match'),'dB','');
stimParams.dBdelta = cellfun(@(c) abs(eval(string(c))),stimParams.dBrange,'uni',1);
stimParams.dBrange = strjoin(stimParams.dBrange,',');
stimParams.inc1dec0 = diff(stimParams.dBdelta)>0;
stimParams.kHzBandwidth = strrep(regexp(stimParams.Pulse,'\d-\d{2}kHz','match'),'kHz','');
stimParams.msDRClen = str2double(strrep(regexp(stimParams.Pulse,'\d*msDRC','match'),'msDRC',''));
stimParams.sPulseLenEphus = str2double(strrep(strrep(regexp(...
    stimParams.Pulse,'_\d{1,2}s_','match'),'_',''),'s',''));
if length(stimParams.dBrange)>1 && contains(pulse.pulseset,'sEach')
    stimParams.sEachContrast = stimParams.sPulseLenEphus;
    stimParams.sPulseLenEphus = stimParams.sPulseLenEphus*2;
end
