%extracts pulse parameters from pulse mat file associated w/ tifs
%containing traces w/ multiple pulses
function [BPNfreqcutoff,BPNdBampl,BPNonsetInPulse,msBPNpulseLen,msStimLen] = extractBPNPulseParams(pulse)

pulseCellStr = {pulse.pulsename}';

BPNonsetInPulse = str2double(string(regexp(pulseCellStr,'(?<=_)\d+(\.\d+)?(?=sOnset)','match')));
BPNfreqcutoff = string(regexp(pulseCellStr,'(?<=BPN_)\d+-\d+(?=kHz)','match'));
BPNdBampl = str2double(string(regexp(pulseCellStr,'(?<=_)\d+(?=dB)','match')));
msBPNpulseLen = str2double(string(regexp(pulseCellStr,'(?<=_)\d+(?=msPulse)','match')));
msStimLen = str2double(string(regexp(pulseCellStr,'(?<=_)\d+(?=msTotal)','match')));

