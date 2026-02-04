function stimParams = extractStimParams(triggerParams,pulse)
%write in DRC logic

%organizes 'TriggerParams' and 'Stim' into relevant params:
stimParams.trigDelay = triggerParams.stimDelay;
pre = regexp(pulse.pulsename,'_\d$','split');
stimParams.rawPulse = pulse.pulsename; %PulseID
stimParams.Pulse = pre{1};
clear pre
stimParams.PTfreq = str2double(strrep(...
    regexp(stimParams.Pulse,'\d{4,5}Hz','match'),'Hz',''));
if isempty(stimParams.PTfreq)
    stimParams.PTfreq = str2double(regexp(regexp(stimParams.Pulse,...
        'stim_\d{1,2}kHz_','match','once'),'\d{1,2}','match','once'))*1000;
end

stimParams.PTampl = str2double(strrep(strrep(...
    regexp(stimParams.Pulse,'_\d{1,2}dB','match'),'dB',''),'_',''));
PTsOnsetPre = strrep(regexp(stimParams.Pulse,'_at_(\ds|\dpt\ds)','match'),'_at_','');
if contains(PTsOnsetPre,'pt')
    stimParams.PTsOnset = str2double(strrep(strrep(PTsOnsetPre,'s',''),'pt','.'));
else
    stimParams.PTsOnset = str2double(strrep(PTsOnsetPre,'s',''));
end
% stimParams.PTmsLen = str2double(strrep(regexp(stimParams.Pulse,'\d{3}ms_at','match'),'ms_at',''));
stimParams.PTmsLen = str2double(strrep(regexp(stimParams.Pulse,'\d{2,3}ms_at','match'),'ms_at',''));

if contains(stimParams.Pulse,'DRC')
    stimParams.dBrange = strrep(regexp(stimParams.Pulse,'\d{2}-\d{2}dB','match'),'dB','');
    if length(stimParams.dBrange)>1
        stimParams.dBdelta = cellfun(@(c) abs(eval(string(c))),stimParams.dBrange,'uni',1);
        stimParams.dBdelta = stimParams.dBdelta(1);
        stimParams.dBrange = strjoin(stimParams.dBrange,',');
    else
        stimParams.dBdelta = abs(eval(string(stimParams.dBrange)));
    end
    
    stimParams.kHzBandwidth = strrep(regexp(stimParams.Pulse,'\d-\d{2}kHz','match'),'kHz','');
    if isempty(stimParams.kHzBandwidth) && (strcmp(stimParams.dBrange,'45-55') || strcmp(stimParams.dBrange,'35-65'))
        stimParams.kHzBandwidth = '5-25';
    elseif isempty(stimParams.kHzBandwidth) && ~(strcmp(stimParams.dBrange,'45-55') || strcmp(stimParams.dBrange,'35-65'))
        error('Could not determine DRC stimulus bandwidth')
    end
    
    stimParams.msDRClen = str2double(strrep(regexp(stimParams.Pulse,'\d*msDRC','match'),'msDRC',''));
    stimParams.sPulseLenEphus = str2double(strrep(strrep(regexp(...
        stimParams.Pulse,'_\d{1,2}s_','match'),'_',''),'s',''));
    if length(stimParams.dBrange)>1 && contains(pulse.pulseset,'sEach')
        stimParams.sPulseLenEphus = stimParams.sPulseLenEphus*2;    
    end
end

end









