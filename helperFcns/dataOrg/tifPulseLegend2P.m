function pulseLegend2P = tifPulseLegend2P(tif2Pdir)
if ~isfolder(tif2Pdir)
    warning('directory non-existent')
    tif2Pdir = uigetdir();
end
dList = dir(tif2Pdir);

pulses = dList(contains({dList.name},'_Pulses.mat'));

for f = 1:length(pulses)
    pulses(f).tif = strrep(pulses(f).name,'_Pulses.mat','.tif');
    load(fullfile(pulses(f).folder,pulses(f).name),'pulse')
    if size(pulse,2)>1
        pulses(f).pulseName = pulse(1).pulseset;
        pulses(f).pulseSet = pulse(1).pulseset;
    elseif endsWith(pulse.pulsename,{'_1','_2','_3','_4','_5','_6','_7','_8','_9'})
        pulses(f).pulseName = pulse.pulsename(1:end-2);
        pulses(f).pulseSet = pulse.pulseset;
    else
        pulses(f).pulseName = pulse.pulsename;
        pulses(f).pulseSet = pulse.pulseset;
    end
    
    clear pulse
end

pulseLegend2P = pulses;



