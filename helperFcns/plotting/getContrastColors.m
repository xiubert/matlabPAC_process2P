function colors = getContrastColors()
%COLORS for sound level contrast plots
colors.lohiPre = [0 0.4510 0.7412; ... %light blue: 0    0.4510    0.7412
    0.8510 0.3294 0.1020]; %light red/orange: 0.8510    0.3294    0.1020
colors.lohiPost = [0 0.3020 0.4902;... %dark blue: 0.0039    0.3216    0.5216
    0.5882 0.2314 0.0784]; %dark red/orange: 0.7020    0.2706    0.0824
colors.lohiTracePre = [0.7294 0.8745 1.0000;...
    1.0000 0.6941 0.5412];
colors.lohiTracePost = colors.lohiTracePre-0.2;
colors.ratio = [0.6510 0.6510 0.6510;...
    0.1490 0.1490 0.1490];  
end