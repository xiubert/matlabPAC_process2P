%plotCohortData

%% LOAD DATA
clearvars;close all;clc;

load('D:\Data\CaMKII_dataTable.mat')
load('D:\Data\CaMKII_params.mat')

%% dFF
% %{
clc; clearvars -except Tinput params
params.permIters = 100000;
tBaseDRC = [-1.2 0];
tBasePT{1} = [0.2 0.5];
tBasePT{2} = [1 2];

Tinput.t_F = rowfun(@(F,fr,trigDelay) ...
    {([1:size(F,2)].*(1/fr))-trigDelay},...
    Tinput,'InputVariables',{'F','frameRate','trigDelay'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

%dFF re DRC
Tinput.dFF_DRC = rowfun(@(F,t) ...
    {dFoFcalc(F,[find((t>=tBaseDRC(1) & t<=tBaseDRC(2)),1,'first')...
    find((t>=tBaseDRC(1) & t<=tBaseDRC(2)),1,'last')],1)},...
    Tinput,'InputVariables',{'F','t_F'},...
    'ExtractCellContents',true,'OutputFormat','uniform');
Tinput.t_dFF_DRC = rowfun(@(t) ...
    {t(find((t>=tBaseDRC(1) & t<=tBaseDRC(2)),1,'first'):end)},...
    Tinput,'InputVariables',{'t_F'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

%dFF re PT
Tinput.dFF_PT = rowfun(@(F,t,PTonset) ...
    {dFoFcalc(F,[find((t>=tBasePT{round(PTonset)}(1) & t<=tBasePT{round(PTonset)}(2)),1,'first')...
    find((t>=tBasePT{round(PTonset)}(1) & t<=tBasePT{round(PTonset)}(2)),1,'last')],1)},...
    Tinput,'InputVariables',{'F','t_F','PTsOnset'},...
    'ExtractCellContents',true,'OutputFormat','uniform');
Tinput.t_dFF_PT = rowfun(@(t,PTonset) ...
    {t(find((t>=tBasePT{round(PTonset)}(1) & t<=tBasePT{round(PTonset)}(2)),1,'first'):end)},...
    Tinput,'InputVariables',{'t_F','PTsOnset'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

% %peak dFF @PT (sig)
tmp = rowfun(@(dFF_PT,t_dFF_PT,PTonset,fr) ...
    pkFcalc(dFF_PT,...
    find(t_dFF_PT>=(PTonset+(1/fr)),1,'first'),...
    params.nFramesPostPulse,params.pkPTsigSD),...
    Tinput,'InputVariables',{'dFF_PT','t_dFF_PT','PTsOnset','frameRate'},...
    'ExtractCellContents',true,'OutputFormat','cell','OutputVariableNames',{'pk','sigPk','sig'});
Tinput.pkPT_sig = tmp(:,1);
Tinput.sigPk = tmp(:,2);
clear tmp

%% PN | PRE: PT onset
%{
close all hidden;clc;clearvars -except Tinput params

figSave = false;
treatment = 'pre';
BW = '5-52';

Tplot = Tinput(contains(Tinput.treatment,treatment) & ...
    contains(Tinput.kHzBandwidth,BW),:);

[G,~,~,~,~] = findgroups(Tplot.animal,Tplot.roiID,Tplot.dBdelta,Tplot.PTsOnset);
if ~(max(G)==length(G))
    % splitapply(@(x) {vertcat(x{:})},Tplot.pkPT_sig,G)
    error('compute cell average first')
end

[G,idC,idPT] = findgroups(Tplot.dBdelta,Tplot.PTsOnset);
tmp = splitapply(@(x) {cellfun(@nanmean,x,'uni',1)},Tplot.pkPT_sig,G);

mat = splitapply(@(c,pk) {pk{c==10}./pk{c==30}},idC,tmp,findgroups(idPT));
u = splitapply(@(c,pk) nanmean(pk{c==10}./pk{c==30}),idC,tmp,findgroups(idPT));
sem = splitapply(@(c,pk) SEMcalc(pk{c==10}./pk{c==30}),idC,tmp,findgroups(idPT));
clear tmp

[h,p,normBool] = sigDiffCalc(mat{:});
fprintf('sig diff re onset: h = %d | p = %d \n',h,p);
if logical(h)
    [hOn(1),pOn(1)] = sigDiffCalc(mat{1},1);
    [hOn(2),pOn(2)] = sigDiffCalc(mat{2},1);
end
fprintf('sig diff re 0.5 s: h = %d | p = %d \n',hOn(1),pOn(1));
fprintf('sig diff re 2 s: h = %d | p = %d \n',hOn(2),pOn(2));

%perm test
if normBool==0
    try
        pPerm = permutationTest(mat{:},params.permIters,'exact',1);
        iters = 'exact';
    catch
        pPerm = permutationTest(mat{:},params.permIters);
        iters = num2str(params.permIters);
    end
    fprintf('PermTest:  sig diff re onset: h = %d | p = %d | iters: %s \n',pPerm<0.05,pPerm,iters);
end

figure('Name',[strrep(params.cohort,'combined',strcat(BW,'kHz'))...
    '_PTresp_ratio_re_PTonset_re_' treatment]);
b = bar(u,'FaceColor',params.colors.ratio(1,:),'EdgeColor','none');
hold on
errorbar(u,sem,'k-','LineWidth',1.5,'LineStyle','none')
modPlotForPoster(1);
ylabelLowHighPkRespRatio(params.colors);
xticklabels(string(unique(Tplot.PTsOnset)))
xlabel('DRC duration preceding pure tone (s)')
if figSave
    figSaveAsFigEpsPng(gcf);
end

%}

%% PN | PRE: dFF_DRC
%{
close all hidden;clc;clearvars -except Tinput params

figSave = false;
tSustainedDRC = [6 8];
% tSustainedDRC = [1 2];

PTonset = 2;
treatment = 'pre';
BW = '5-52';
plotBoxplot = true;

if plotBoxplot==1
    Tplot = Tinput(contains(Tinput.treatment,treatment) & ...
        contains(Tinput.kHzBandwidth,BW),:);
else
    Tplot = Tinput(contains(Tinput.treatment,treatment) & ...
        contains(Tinput.kHzBandwidth,BW) & ...
        Tinput.PTsOnset==PTonset,:);
end

[G,~,~,idC] = findgroups(Tplot.animal,Tplot.roiID,Tplot.dBdelta);
if ~(max(G)==length(G))
    tmp = splitapply(@(x) {nanmean(vertcat(x{:}),1)},Tplot.dFF_DRC,G);
    [G,idC] = findgroups(idC);
    uDFF_DRC_cell = splitapply(@(x) {cell2mat(x)},tmp,G);
    warning('compute cell average first')
else
    [G,idC] = findgroups(Tplot.dBdelta);
    uDFF_DRC_cell = splitapply(@(x) {cell2mat(x)},cellfun(@(c) nanmean(c,1),Tplot.dFF_DRC,'uni',0),G);
end

uDFF_DRC_cell_u = cell2mat(cellfun(@(c) nanmean(c,1).*100,uDFF_DRC_cell,'uni',0));
uDFF_DRC_cell_sem = cell2mat(cellfun(@(c) SEMcalc(c,1).*100,uDFF_DRC_cell,'uni',0));

tDRC = Tplot.t_dFF_DRC{1};

if plotBoxplot==0
    figure('Name',[strrep(params.cohort,'combined',strcat(BW,'kHz'))...
        '_sustainedDRC_traces_re_' treatment]);
    subplot(2,1,1) %dFF
    fillSEMplot(tDRC,...
        uDFF_DRC_cell_u(idC==10,:),...
        uDFF_DRC_cell_sem(idC==10,:),...
        params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:),...
        params.colors.(...
        strjoin(['lohiTrace' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:));
    hold on
    fillSEMplot(tDRC,...
        uDFF_DRC_cell_u(idC==30,:),...
        uDFF_DRC_cell_sem(idC==30,:),...
        params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:),...
        params.colors.(...
        strjoin(['lohiTrace' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:));
    ylabel('% \DeltaF/F')
    % xlim([-1.2 max(tDRC)])
    xlim([-1.2 8])
    modPlotForPoster(0)
    
    subplot(2,1,2) %signal
    [sigPTcontrastHi,fs] = inspectSignalObject('plot',false,'signalPath',...
        'C:\Data\Rig Software\250kHzPulses\PC_PTinContrast_5-52kHz_25msDRC_10000Hz70dB_0s_delay\25msDRC_5-52kHz_50-60dB_8s_sin2-cos2_DRCramp_250kHz_stim_10000Hz_70dB_400ms_at_2s_2.signal');
    sigPTcontrastLo = inspectSignalObject('plot',false,'signalPath',...
        'C:\Data\Rig Software\250kHzPulses\PC_PTinContrast_5-52kHz_25msDRC_10000Hz70dB_0s_delay\25msDRC_5-52kHz_40-70dB_8s_sin2-cos2_DRCramp_250kHz_stim_10000Hz_70dB_400ms_at_2s_3.signal');
    tSig = [1:length(sigPTcontrastLo)].*(1/fs);
    plot(tSig,sigPTcontrastLo,'Color',params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:));
    hold on
    plot(tSig,sigPTcontrastHi+15000,'Color',params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:));
    % xlim([-1.2 max(tDRC)])
    xlim([-1.2 8])
    hLine(1) = yline(15000);
    hLine(2) = yline(0);
    hLine(1).Color = params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:);
    hLine(2).Color = params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:);
    set(gca,'yticklabel',[]);
    modPlotForPoster(0);
    xlabel('time (s)');
    ax = gca;
    ax.YAxis.TickLength = [0 0];
    if figSave
        figSaveAsFigEpsPng(gcf);
    end
    
end

if plotBoxplot==1
    %sustained DRC
    uDFF_DRCsustained_cell = cellfun(@(c) nanmean(c(:,tDRC>=tSustainedDRC(1) & ...
        tDRC<=tSustainedDRC(2)),2),uDFF_DRC_cell,'uni',0);
    lo = uDFF_DRCsustained_cell{idC==10};
    hi = uDFF_DRCsustained_cell{idC==30};
    
    [hDRC,pDRC,normBool] = sigDiffCalc(lo,hi);
    fprintf('scatter:  lo>hi %d | hi>lo %d | total %d \n',...
        sum((lo>hi)),sum((hi>lo)),sum(all(~isnan([lo hi]),2)));
    fprintf('sustained DRC diff re contrast:  h = %d | p = %d \n',hDRC,pDRC)
    
    %perm test
    if normBool==0
        try
            pPerm = permutationTest(lo,hi,params.permIters,'exact',1);
            iters = 'exact';
        catch
            pPerm = permutationTest(lo,hi,params.permIters);
            iters = num2str(params.permIters);
        end
        fprintf('PermTest:  sustained DRC diff re contrast: h = %d | p = %d | iters: %s \n',...
            pPerm<0.05,pPerm,iters);
    end
    
    figure('Name',sprintf('%s_sustainedDRC_%d-%ds_boxplot_re_%s',...
        strrep(params.cohort,'combined',strcat(BW,'kHz')),tSustainedDRC,treatment));
    
    tmpG = arrayfun(@(r) ones(length(uDFF_DRCsustained_cell{r}),1).*r,[1 2],'uni',0);
    boxplot(vertcat(uDFF_DRCsustained_cell{:}).*100,vertcat(tmpG{:}),...
        'PlotStyle','compact','Colors',params.colors.ratio(1,:))
    % set(gca,'xtick',[]);
    % h=gca; h.XAxis.TickLength = [0 0];
    lowStr = strjoin(arrayfun(@(x) num2str(x),params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:),'UniformOutput',false),',');
    highStr = strjoin(arrayfun(@(x) num2str(x),params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:),'UniformOutput',false),',');
    set(gca,'xtick',1:2,...
        'xticklabel',{['\color[rgb]{' lowStr '} LOW'],...
        ['\color[rgb]{' highStr '} HIGH']})
    modPlotForPoster(0)
    ylabel('mean % \DeltaF/F')
    
    a = get(get(gca,'children'),'children');   % Get the handles of all the objects
    % t = get(a,'tag');   % List the names of all the objects
    % box1 = a(7);   % The 7th object is the first box
    % set(box1, 'Color', 'g');   % Set the color of the first box to green
    set(a(9),'Color',params.colors.lohiPre(2,:))
    set(a(10),'Color',params.colors.lohiPre(1,:))
    % set(gca,'xtick',1:2,...
    %         'xticklabel',{'',''})
    set(gcf,'Position',[150 244.5000 369 420])
    % if figSave
    %     figSaveAsFigEpsPng(gcf);
    % end
end

%}

%% PN | PRE: dFF_PT | scatterplot and traces re contrast
%{
close all hidden;clc;clearvars -except Tinput params

figSave = false;
treatment = 'pre';
BW = '5-52'; %get rid of this w/ SOM and PV
PTonset = 2;

Tplot = Tinput(contains(Tinput.treatment,treatment) & ...
    contains(Tinput.kHzBandwidth,BW) & ...
    Tinput.PTsOnset==PTonset,:);

[G,~,~,~] = findgroups(Tplot.animal,Tplot.roiID,Tplot.dBdelta);
if ~(max(G)==length(G))
    % splitapply(@(x) {vertcat(x{:})},Tplot.dFF_DRC,G)
    error('compute cell average first')
end

[G,idC] = findgroups(Tplot.dBdelta);

%cell avg pk resp re contrast
uPkCell = splitapply(@(x) {cellfun(@nanmean,x,'uni',1)},Tplot.pkPT_sig,G);
lo = uPkCell{idC==10};
hi = uPkCell{idC==30};

if logical(sigDiffCalc(lo,hi))   
    regPlotColors(1,:) = params.colors.lohiPre((~(nanmean(lo)>...
        nanmean(hi)))+1,:);
    regPlotColors(2,:) = params.colors.lohiTracePre((~(nanmean(lo)>...
        nanmean(hi)))+1,:);
else
    regPlotColors(1,:) = [0,0,0];
    regPlotColors(2,:) = [0.6510    0.6510    0.6510];
end
fprintf('lo>hi %d | hi>lo %d | total %d \n',...
    sum((lo>hi)),sum((hi>lo)),sum(all(~isnan([lo hi]),2)));
[hLOHI,pLOHI,normBool] = sigDiffCalc(lo,hi);
fprintf('lo v. hi sig diff? h = %d | p = %d \n',hLOHI,pLOHI)

%perm test
if normBool~=1
    try
        pPerm = permutationTest(lo,hi,params.permIters,'exact',1);
        iters = 'exact';
    catch
        pPerm = permutationTest(lo,hi,params.permIters);
        iters = num2str(params.permIters);
    end
    fprintf('PermTest:  PTonset diff re contrast: h = %d | p = %d | iters: %s \n',pPerm<0.05,pPerm,iters);
end

%pk resp scatter
[hScatter, hRegLine, hRegLineCI, mdl, ciRegLine, hRef] = ...
    regPlot(lo.*100,...
    hi.*100,...
    'regLine',true,...
    'logScale',false,...
    'intercept',false,...
    'colors',regPlotColors,...
    'fName',[strrep(params.cohort,'combined',strcat(BW,'kHz'))...
    '_PTresp_scatter_re_contrast_re_' treatment]);
xlabel({'peak pure tone % \DeltaF/F','in low contrast'},...
    'Color',params.colors.lohiPre(1,:),'interpreter', 'tex')
ylabel({'peak pure tone % \DeltaF/F','in high contrast'},...
    'Color',params.colors.lohiPre(2,:),'interpreter', 'tex')
legend('off')
modPlotForPoster(0);
if figSave
    figSaveAsFigEpsPng(gcf);
end

mCI = mdl.coefCI;
fprintf('does slope confidence interval (%.2f - %.2f) include 1? %d \n',mdl.coefCI,(1>mCI(1) & 1<mCI(2)))

%pk resp trace
sigTrace = rowfun(@(dFF_PT,sig) ...
    {dFF_PT(sig,:)},Tplot,...
    'InputVariables',{'dFF_PT','sigPk'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

uTrace = cell2mat(splitapply(@(x) ...
    {nanmean(cell2mat(cellfun(@(c) nanmean(c,1),x,'uni',0)),1).*100},sigTrace,G));
semTrace = cell2mat(splitapply(@(x) ...
    {SEMcalc(cell2mat(cellfun(@(c) nanmean(c,1),x,'uni',0)),1).*100},sigTrace,G));
tDFF = Tplot.t_dFF_PT{1};

figure('Name',[strrep(params.cohort,'combined',strcat(BW,'kHz'))...
    '_PTresp_traces_re_contrast_re_' treatment]);
fillSEMplot(tDFF,uTrace(idC==10,:),semTrace(idC==10,:),...
    params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(idC==10,:),...
    params.colors.(strjoin(['lohiTrace' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(idC==10,:));
hold on
fillSEMplot(tDFF,uTrace(idC==30,:),semTrace(idC==30,:),...
    params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(idC==30,:),...
    params.colors.(strjoin(['lohiTrace' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(idC==30,:));
xlim([1.5 3.5])
% xlim([0 2]) for 0.5s onset

modPlotForPoster(0)
xlabel('time (s)')
ylabel('% \DeltaF/F')
if figSave
    figSaveAsFigEpsPng(gcf);
end
%}

%% PN | PRE: dFF_PT | BW, OCT
%{
close all hidden;clc;clearvars -except Tinput params T 

figSave = false;
PTonset = 2;

BW = {'5-25','5-52'};

%get ratio by bandwidth
ratioBW = cell(1,2);
for iBW = 1:2
    Tplot = Tinput(contains(Tinput.kHzBandwidth,BW{iBW}) & ...
        Tinput.PTsOnset==PTonset & ...
        contains(Tinput.treatment,'pre'),:);
    
    [G,iA,iCell,iC] = findgroups(Tplot.animal,Tplot.roiID,Tplot.dBdelta);
    tmp = splitapply(@(x) nanmean(vertcat(x{:})),Tplot.pkPT_sig,G);
    tmp2 = tmp(iC==10)./tmp(iC==30);
        
    ratioBW{iBW} = tmp2;
    fprintf('%s: %d cells, %d mice\n',BW{iBW},length(ratioBW{iBW}),numel(unique(Tplot.animal)))
    clear tmp*
end

%ratio re bandwidth
[hBW,pBW,BWnorm] = sigDiffCalc(ratioBW{1},ratioBW{2});
fprintf('diff re bandwidth: h = %d | p = %d \n',hBW,pBW);
[h,p] = sigDiffCalc(ratioBW{1},1);
fprintf('%s: h = %d | p = %d\n',BW{1},h,p)
[h,p] = sigDiffCalc(ratioBW{2},1);
fprintf('%s: h = %d | p = %d\n',BW{2},h,p)

if BWnorm~=1
    %perm test
    try
        pPerm = permutationTest(ratioBW{1},ratioBW{2},params.permIters,'exact',1);
        iters = 'exact';
    catch
        pPerm = permutationTest(ratioBW{1},ratioBW{2},params.permIters);
        iters = num2str(params.permIters);
    end
    fprintf('PermTest:  diff re bandwidth: h = %d | p = %d | iters: %s \n',pPerm<0.05,pPerm,iters);
end

%box plot
tmp = arrayfun(@(x) ones(length(ratioBW{x}),1).*x,1:2,'uni',0);
figure('Name',[params.cohort '_PTresp_ratio_re_DRCbandwidth']);
boxplot(vertcat(ratioBW{:}),vertcat(tmp{:}),...
    'PlotStyle','compact','Colors',params.colors.ratio(1,:));
hold on
modPlotForPoster(1)
xlabel({'','DRC bandwidth (kHz)'})
set(gca,'xtick',1:length(BW),...
        'xticklabel',BW)
ylabelLowHighPkRespRatio(params.colors);
if figSave
    figSaveAsFigEpsPng(gcf);
end

%ratio re PT/BF  within / outside octave
Tplot = Tinput(contains(Tinput.kHzBandwidth,BW{2}) & ...
        Tinput.PTsOnset==PTonset & ...
        contains(Tinput.treatment,'pre'),:);

withinOct = rowfun(@(fq,BF) abs(log2(fq./BF))>1,...
    Tplot,'ExtractCellContents',1,...
    'OutputFormat','uniform',...
    'InputVariables',{'PTfreq','BFsigUDB'});

[G,idC,idOct] = findgroups(Tplot.dBdelta,withinOct);

tmp = splitapply(@(x) {cellfun(@nanmean,x,'uni',1)},Tplot.pkPT_sig,G);
G2 = findgroups(idOct);
pkReOct = splitapply(@(pk,contrast) {pk{contrast==10}./pk{contrast==30}},tmp,idC,G2);

[hOct,pOct,octNorm] = sigDiffCalc(pkReOct{1},pkReOct{2});
fprintf('diff re oct: h = %d | p = %d \n',hOct,pOct);
[h,p] = sigDiffCalc(pkReOct{1},1);
fprintf('within: h = %d | p = %d\n',h,p)
[h,p] = sigDiffCalc(pkReOct{2},1);
fprintf('outside: h = %d | p = %d\n',h,p)

%perm test
if octNorm~=1
    try
        pPerm = permutationTest(pkReOct{1},pkReOct{2},params.permIters,'exact',1);
        iters = 'exact';
    catch
        pPerm = permutationTest(pkReOct{1},pkReOct{2},params.permIters);
        iters = num2str(params.permIters);
    end
    fprintf('PermTest:  diff re oct: h = %d | p = %d | iters: %s \n',pPerm<0.05,pPerm,iters);
end

tmp = arrayfun(@(x) ones(length(pkReOct{x}),1).*x,1:2,'uni',0);
figure('name',[params.cohort '_PTresp_ratio_re_BF-PToctDif']);
boxplot(vertcat(pkReOct{:}),vertcat(tmp{:}),...
    'PlotStyle','compact','Colors',params.colors.ratio(1,:));
hold on
modPlotForPoster(1)
xlabel({'','octave of cell BF'})
set(gca,'xtick',1:2,...
        'xticklabel',{'within','outside'})
ylabelLowHighPkRespRatio(params.colors);
if figSave
    figSaveAsFigEpsPng(gcf);
end


%}

%% PN | POST: dFF_DRC
%{
close all hidden;clc;clearvars -except Tinput params
figSave = false;
tSustainedDRC = [6 8];
% tSustainedDRC = [1 2];

PTonset = 2;
treatment = 'postZX1';
plotBoxplot = true;

if plotBoxplot==1
    Tplot = Tinput(contains(Tinput.treatment,treatment),:);  %for boxplot
else
    Tplot = Tinput(contains(Tinput.treatment,treatment) & ...
        Tinput.PTsOnset==PTonset,:); %for 2s onset trace plot
end

[G,~,~,idC] = findgroups(Tplot.animal,Tplot.roiID,Tplot.dBdelta);
if ~(max(G)==length(G))
    tmp = splitapply(@(x) {nanmean(vertcat(x{:}),1)},Tplot.dFF_DRC,G);
    [G,idC] = findgroups(idC);
    uDFF_DRC_cell = splitapply(@(x) {cell2mat(x)},tmp,G);
%     error('compute cell average first')
else
    [G,idC] = findgroups(Tplot.dBdelta);
    uDFF_DRC_cell = splitapply(@(x) {cell2mat(x)},cellfun(@(c) nanmean(c,1),Tplot.dFF_DRC,'uni',0),G);
end

uDFF_DRC_cell_u = cell2mat(cellfun(@(c) nanmean(c,1).*100,uDFF_DRC_cell,'uni',0));
uDFF_DRC_cell_sem = cell2mat(cellfun(@(c) SEMcalc(c,1).*100,uDFF_DRC_cell,'uni',0));

tDRC = Tplot.t_dFF_DRC{1};

if plotBoxplot~=1
    figure('Name',[params.cohort '_sustainedDRC_traces_re_' treatment]);
    subplot(2,1,1) %dFF
    fillSEMplot(tDRC,...
        uDFF_DRC_cell_u(idC==10,:),...
        uDFF_DRC_cell_sem(idC==10,:),...
        params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:),...
        params.colors.(...
        strjoin(['lohiTrace' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:));
    hold on
    fillSEMplot(tDRC,...
        uDFF_DRC_cell_u(idC==30,:),...
        uDFF_DRC_cell_sem(idC==30,:),...
        params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:),...
        params.colors.(...
        strjoin(['lohiTrace' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:));
    ylabel('% \DeltaF/F')
    xlim([-1.2 8])
    modPlotForPoster(0)
    
    subplot(2,1,2) %signal
    [sigPTcontrastHi,fs] = inspectSignalObject('plot',false,'signalPath',...
        'C:\Data\Rig Software\250kHzPulses\PC_PTinContrast_5-52kHz_25msDRC_10000Hz70dB_0s_delay\25msDRC_5-52kHz_50-60dB_8s_sin2-cos2_DRCramp_250kHz_stim_10000Hz_70dB_400ms_at_2s_2.signal');
    sigPTcontrastLo = inspectSignalObject('plot',false,'signalPath',...
        'C:\Data\Rig Software\250kHzPulses\PC_PTinContrast_5-52kHz_25msDRC_10000Hz70dB_0s_delay\25msDRC_5-52kHz_40-70dB_8s_sin2-cos2_DRCramp_250kHz_stim_10000Hz_70dB_400ms_at_2s_3.signal');
    tSig = [1:length(sigPTcontrastLo)].*(1/fs);
    plot(tSig,sigPTcontrastLo,'Color',params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:));
    hold on
    plot(tSig,sigPTcontrastHi+15000,'Color',params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:));
    xlim([-1.2 8])
    hLine(1) = yline(15000);
    hLine(2) = yline(0);
    hLine(1).Color = params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:);
    hLine(2).Color = params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:);
    set(gca,'yticklabel',[]);
    modPlotForPoster(0);
    xlabel('time (s)');
    ax = gca;
    ax.YAxis.TickLength = [0 0];
    if figSave
        figSaveAsFigEpsPng(gcf);
    end
end

if plotBoxplot==1
    %sustained DRC
    uDFF_DRCsustained_cell = cellfun(@(c) nanmean(c(:,tDRC>=tSustainedDRC(1) & ...
        tDRC<=tSustainedDRC(2)),2),uDFF_DRC_cell,'uni',0);
    lo = uDFF_DRCsustained_cell{idC==10};
    hi = uDFF_DRCsustained_cell{idC==30};
    
    [hDRC,pDRC,normBool] = sigDiffCalc(lo,hi);
    fprintf('lo>hi %d | hi>lo %d | total %d \n',...
        sum((lo>hi)),sum((hi>lo)),sum(all(~isnan([lo hi]),2)));
    fprintf('sustained DRC diff re contrast:  h = %d | p = %d \n',hDRC,pDRC)
    
    %perm test
    if normBool==0
        try
            pPerm = permutationTest(lo,hi,params.permIters,'exact',1);
            iters = 'exact';
        catch
            pPerm = permutationTest(lo,hi,params.permIters);
            iters = num2str(params.permIters);
        end
        fprintf('PermTest:  sustained DRC diff re contrast: h = %d | p = %d | iters: %s \n',pPerm<0.05,pPerm,iters);
    end
    
    figure('Name',sprintf('%s_sustainedDRC_%d-%ds_boxplot_re_%s',...
        params.cohort,tSustainedDRC,treatment))
    tmpG = arrayfun(@(r) ones(length(uDFF_DRCsustained_cell{r}),1).*r,[1 2],'uni',0);
    boxplot(vertcat(uDFF_DRCsustained_cell{:}).*100,vertcat(tmpG{:}),...
        'PlotStyle','compact','Colors',params.colors.ratio(1,:))
    modPlotForPoster(0)
    ylabel('mean % \DeltaF/F')
    % h=gca; h.XAxis.TickLength = [0 0];
    a = get(get(gca,'children'),'children');   % Get the handles of all the objects
    set(a(9),'Color',params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:))
    set(a(10),'Color',params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:))
    lowStr = strjoin(arrayfun(@(x) num2str(x),params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:),'UniformOutput',false),',');
    highStr = strjoin(arrayfun(@(x) num2str(x),params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:),'UniformOutput',false),',');
    set(gca,'xtick',1:2,...
        'xticklabel',{['\color[rgb]{' lowStr '} LOW'],...
        ['\color[rgb]{' highStr '} HIGH']})
    % set(gca,'xtick',1:2,...
    %         'xticklabel',{'',''})
    set(gcf,'Position',[150 244.5000 369 420])
    % if figSave
    %     figSaveAsFigEpsPng(gcf);
    % end
end

%}

%% PN | PRE --> POST:  dFF_PT
%{
close all hidden;clc;clearvars -except Tinput params

figSave = false;
DRUG = 'ZX1';
% DRUG = 'ACSF';
PTonset = 2;

Tplot = Tinput(Tinput.PTsOnset==PTonset & ...
    contains(Tinput.treatment,DRUG),:);

[G,~,~,idT,idC] = findgroups(Tplot.animal,Tplot.roiID,Tplot.treatment,Tplot.dBdelta);
if ~(max(G)==length(G))
    tmp = splitapply(@(x) {nanmean(vertcat(x{:}),1)},Tplot.pkPT_sig,G); %cell average pkPT
    [G,idC,idT] = findgroups(idC,idT);
    uPkPT = splitapply(@(x) {cell2mat(x)},tmp,G);
    %     error('compute cell average first')
else
    [G,idT,idC] = findgroups(Tplot.treatment,Tplot.dBdelta);
    uPkPT = splitapply(@(x) {cellfun(@nanmean,x)},Tplot.pkPT_sig,G);
end

tmpT = table(idT,idC,uPkPT);
tmpT.n = cellfun(@(c) sum(~isnan(c)),tmpT.uPkPT,'uni',1);

prepost = flipud(unique(tmpT.idT));
pp = strrep(prepost,DRUG,'');
mm = {'min','max'};
lh = {'low','high'};

%peak response stats
dT = array2table(horzcat(tmpT.uPkPT{:}),'VariableNames',...
    strcat(tmpT.idT,'_',lh(findgroups(tmpT.idC))'));
withinStruct = table([1 2 1 2]',[2 2 1 1]','VariableNames',{'Contrast','PrePost'});
withinStruct.Contrast = categorical(withinStruct.Contrast);
withinStruct.PrePost = categorical(withinStruct.PrePost);
%no between subjects, constant term on each w/ a-z~1
rm = fitrm(dT, ['post' DRUG '_low,post' DRUG '_high,pre' DRUG '_low,pre' DRUG '_high~1'],'WithinDesign',withinStruct);
ranovatable = ranova(rm,'WithinModel','Contrast*PrePost');
disp(['ranova main effect p = ' num2str(ranovatable{1,'pValue'})])
disp(['ranova contrast effect p = ' num2str(ranovatable{3,'pValue'})])
disp(['ranova PrePost effect p = ' num2str(ranovatable{5,'pValue'})])
disp(['ranova Interaction p = ' num2str(ranovatable{7,'pValue'})])

statsDrugByContrast = multcompare(rm,'PrePost','By','Contrast','ComparisonType','bonferroni');
statsContrastByDrug = multcompare(rm,'Contrast','By','PrePost','ComparisonType','bonferroni');

%bar plot of peak responses
u = cellfun(@nanmean,tmpT.uPkPT,'uni',1).*100;
sem = cellfun(@SEMcalc,tmpT.uPkPT,'uni',1).*100;
sortIDX = [find(tmpT.idC==min(tmpT.idC) & contains(tmpT.idT,'pre')) ...
    find(tmpT.idC==min(tmpT.idC) & contains(tmpT.idT,'post')); ...
    find(tmpT.idC==max(tmpT.idC) & contains(tmpT.idT,'pre')) ...
    find(tmpT.idC==max(tmpT.idC) & contains(tmpT.idT,'post'))];
u = u(sortIDX);
sem = sem(sortIDX);

figure('Name',[params.cohort '_PTresp_bar_re_contrast_re_' DRUG]);
b = bar(u,'FaceColor','flat','EdgeColor','none');
b(1).CData = params.colors.lohiPre;
b(2).CData = params.colors.lohiPost;
hold on
groupBarPlotErrorBar(u,sem);
lowStr = strjoin(arrayfun(@(x) num2str(x),params.colors.lohiPre(1,:),'UniformOutput',false),',');
highStr = strjoin(arrayfun(@(x) num2str(x),params.colors.lohiPre(2,:),'UniformOutput',false),',');

xticklabels({['\color[rgb]{' lowStr '} LOW'],...
    ['\color[rgb]{' highStr '} HIGH']});
ylabel('% \DeltaF/F')
text(1+b(1).XOffset,2.5,'CTRL',...
    'Rotation',90,'Color','w','FontWeight','bold','FontSize',16)
text(1+(-1.*(b(1).XOffset)),2.5,strrep(prepost{1},'pre',''),...
    'Rotation',90,'Color','w','FontWeight','bold','FontSize',16)
text(2+b(1).XOffset,2.5,'CTRL',...
    'Rotation',90,'Color','w','FontWeight','bold','FontSize',16)
text(2+(-1.*(b(1).XOffset)),2.5,strrep(prepost{1},'pre',''),...
    'Rotation',90,'Color','w','FontWeight','bold','FontSize',16)
modPlotForPoster(0)
set(gca,'Position',[0.2000 0.100 0.6000 0.7150])
if figSave
    figSaveAsFigEpsPng(gcf);
end

%low/high pre and post
for ppID = 1:length(pp)
    clear lo hi
    lo = tmpT{contains(tmpT.idT,pp{ppID}) & ...
        tmpT.idC==10,'uPkPT'}{:};
    hi = tmpT{contains(tmpT.idT,pp{ppID}) & ...
        tmpT.idC==30,'uPkPT'}{:};
    
    if statsContrastByDrug{double(statsContrastByDrug.PrePost)==ppID,{'pValue'}}(1)<0.05
        regPlotColors(1,:) = params.colors.(['lohi' strrep(pp{ppID},'p','P')])((~(nanmean(lo)>nanmean(hi)))+1,:);
        regPlotColors(2,:) = params.colors.(['lohiTrace' strrep(pp{ppID},'p','P')])((~(nanmean(lo)>nanmean(hi)))+1,:);
    else
        regPlotColors(1,:) = [0,0,0];
        regPlotColors(2,:) = [0.6510    0.6510    0.6510];
    end
    
    [~, ~, ~, mdl] = regPlot(lo.*100,...
        hi.*100,...
        'fName',[params.cohort '_PTresp_scatter_' pp{ppID} '_' DRUG '_re_contrast'],...
        'colors',regPlotColors,'logScale',false,'regLine',true);
    xlabel({'peak pure tone % \DeltaF/F','in low contrast'},...
        'Color',params.colors.(['lohi' strrep(pp{ppID},'p','P')])(1,:),'interpreter', 'tex')
    ylabel({'peak pure tone % \DeltaF/F','in high contrast'},...
        'Color',params.colors.(['lohi' strrep(pp{ppID},'p','P')])(2,:),'interpreter', 'tex')
    modPlotForPoster(0);
    fprintf('lo>hi %d | hi>lo %d | total %d \n',...
        sum((lo>hi)),sum((hi>lo)),sum(all(~isnan([lo hi]),2)));
    [hLOHI,pLOHI,normBool] = sigDiffCalc(lo,hi);
    fprintf('lo v. hi sig diff? h = %d | p = %d \n',hLOHI,pLOHI)
    
    %perm test
    if normBool~=1
        try
            pPerm = permutationTest(lo,hi,params.permIters,'exact',1);
            iters = 'exact';
        catch
            pPerm = permutationTest(lo,hi,params.permIters);
            iters = num2str(params.permIters);
        end
        fprintf('%s | PermTest:  PTresp diff re contrast: h = %d | p = %d | iters: %s \n',...
            pp{ppID},pPerm<0.05,pPerm,iters);
    end
    
    if figSave
        figSaveAsFigEpsPng(gcf);
    end
    
    mCI = mdl.coefCI;
    fprintf('%s | does slope confidence interval (%.2f - %.2f) include 1? %d \n',...
        pp{ppID},mdl.coefCI,(1>mCI(1) & 1<mCI(2)))
end


%pre/post low and high
for cID = 1:length(lh)
    clear pre post
    if statsDrugByContrast{double(statsDrugByContrast.Contrast)==cID,'pValue'}(1)<0.05
        regPlotColors(1,:) = params.colors.lohiPre(cID,:);
        regPlotColors(2,:) = params.colors.lohiTracePre(cID,:);
    else
        regPlotColors(1,:) = [0,0,0];
        regPlotColors(2,:) = [0.6510    0.6510    0.6510];
    end
    pre = tmpT{contains(tmpT.idT,'pre') & ...
        eval(['tmpT.idC==' mm{cID} '(tmpT.idC)']),'uPkPT'}{:};
    post = tmpT{contains(tmpT.idT,'post') & ...
        eval(['tmpT.idC==' mm{cID} '(tmpT.idC)']),'uPkPT'}{:};
    
    [~, ~, ~, mdl] = regPlot(pre.*100,...
        post.*100,...
        'fName',[params.cohort '_PTresp_scatter_' upper(lh{cID}) '_contrast_re_' DRUG],...
        'colors',regPlotColors,'logScale',false,'regLine',true);
    xlabel({'peak pure tone % \DeltaF/F',['in ' lh{cID} ' contrast']},...
        'Color',params.colors.lohiPre(cID,:),'interpreter', 'tex');
    ylabel({'peak pure tone % \DeltaF/F',['in ' lh{cID} ' contrast']},...
        'Color',params.colors.lohiPost(cID,:),'interpreter', 'tex');
    modPlotForPoster(0);       
    if figSave
        figSaveAsFigEpsPng(gcf);
    end
    fprintf('%s: pre/post n = %d\n',lh{cID},sum(all(~isnan([pre post]),2)))
    [~,~,normBool] = sigDiffCalc(pre,post);
    
    %perm test
    if normBool~=1
        try
            pPerm = permutationTest(pre,post,params.permIters,'exact',1);
            iters = 'exact';
        catch
            pPerm = permutationTest(pre,post,params.permIters);
            iters = num2str(params.permIters);
        end
        fprintf('%s | PermTest:  PTresp diff re %s: h = %d | p = %d | iters: %s \n',lh{cID},DRUG,pPerm<0.05,pPerm,iters);
    end
    
    mCI = mdl.coefCI;
    fprintf('%s | does slope confidence interval (%.2f - %.2f) include 1? %d \n',lh{cID},mdl.coefCI,(1>mCI(1) & 1<mCI(2)))
end

%RATIO pre & post bar plot
uPkByTreatment = splitapply(@(x) {x},tmpT.uPkPT,findgroups(tmpT.idT));
pkRespT = table(uPkByTreatment,...
    unique(tmpT.idT,'stable'),...
    repmat({unique(tmpT.idC,'stable')},max(findgroups(tmpT.idT)),1),...
    'VariableNames',{'uPkPT','treatment','contrast'});
pkRespT.pkRatio = arrayfun(@(r) r.uPkPT{r.contrast==min(r.contrast)}./...
    r.uPkPT{r.contrast==max(r.contrast)},...
    table2struct(pkRespT),'uni',0);
%sort pre/post (opposite alphabetical)
[~,tmp] = sort(pkRespT.treatment);
pkRespT = pkRespT(flipud(tmp),:);
u = cellfun(@(c) nanmean(c),pkRespT.pkRatio,'uni',1);
sem = cellfun(@(c) SEMcalc(c),pkRespT.pkRatio,'uni',1);

figure('Name',[params.cohort '_PTresp_ratio_re_' DRUG]);
b = bar(u,'FaceColor','Flat','EdgeColor','none');
b.CData = params.colors.ratio;
hold on
errorbar(u,sem,...
    'k-','LineWidth',1.5,'LineStyle','none')
hLabel = ylabelLowHighPkRespRatio(params.colors);
modPlotForPoster(1)
xticklabels('')
if figSave
    figSaveAsFigEpsPng(gcf);
end

fprintf('ratio n = %d \n',sum(all(~isnan(horzcat(pkRespT.pkRatio{:})),2)))
if adtest(vertcat(pkRespT.pkRatio{:}))
    disp(['ratio pre/post p = ' num2str(signrank(pkRespT.pkRatio{:}))])
    disp(['ratio pre p = ' num2str(signrank(pkRespT.pkRatio{1},1))])
    disp(['ratio post p = ' num2str(signrank(pkRespT.pkRatio{2},1))])
end

%perm test
if adtest(vertcat(pkRespT.pkRatio{:}))
    try
        pPerm = permutationTest(pkRespT.pkRatio{:},params.permIters,'exact',1);
        iters = 'exact';
    catch
        pPerm = permutationTest(pkRespT.pkRatio{:},params.permIters);
        iters = num2str(params.permIters);
    end
    fprintf('PermTest:  PTresp ratio re %s: h = %d | p = %d | iters: %s \n',DRUG,pPerm<0.05,pPerm,iters);
end

%TRACES pre/post
[G,gTreatment,gContrast] = findgroups(Tplot.treatment,Tplot.dBdelta);
dFF_PT_sig = rowfun(@(pk,sig) {pk(sig,:)},Tplot,...
    'ExtractCellContents',true,...
    'InputVariables',{'dFF_PT','sigPk'},...
    'OutputFormat','uniform');
tmp = splitapply(@(x) {cell2mat(cellfun(@(c) nanmean(c,1),x,'uni',0))},dFF_PT_sig,G);
n = cellfun(@(c) sum(all(~isnan(c),2)),tmp,'uni',1);
u = cell2mat(cellfun(@(c) nanmean(c,1),tmp,'uni',0)).*100;
sem = cell2mat(cellfun(@(c) SEMcalc(c,1),tmp,'uni',0)).*100;
tPT = Tplot.t_dFF_PT{1};

pp = strrep(prepost,DRUG,'');


figure('Name',[params.cohort '_PTresp_traces_re_' DRUG]);
for ppID = 1:length(pp)
    subplot(1,2,ppID)
    for mmID = 1:length(mm)
        clear idx
        idx = eval(['gContrast==' mm{mmID} '(gContrast)']) & contains(gTreatment,pp{ppID});
        
        fillSEMplot(tPT,u(idx,:),sem(idx,:),...
            params.colors.(['lohi' strrep(pp{ppID},'p','P')])(mmID,:),...
            params.colors.(['lohiTrace' strrep(pp{ppID},'p','P')])(mmID,:));
        hold on
    end
    yl(ppID,:) = ylim;
    xlim([1.5 3.5])
    modPlotForPoster(0)
xlabel('time (s)')
ylabel('% \DeltaF/F')
end
if figSave
    figSaveAsFigEpsPng(gcf);
end




%}

%% PN | PRE: DRC OFFSET (supplemental)
%{
close all hidden;clc;clearvars -except Tinput params

figSave = true;
tBaseDRCoff = [7 8];
drcOFFtBin = [8.2 8.4];
treatment = 'pre';

Tplot = Tinput(contains(Tinput.treatment,treatment) & ...
    Tinput.sPulseLenEphus==8,:);

%dFF re DRC
Tplot.dFF_DRCoff = rowfun(@(F,t) ...
    {dFoFcalc(F,[find((t>=tBaseDRCoff(1) & t<=tBaseDRCoff(2)),1,'first')...
    find((t>=tBaseDRCoff(1) & t<=tBaseDRCoff(2)),1,'last')],1)},...
    Tplot,'InputVariables',{'F','t_F'},...
    'ExtractCellContents',true,'OutputFormat','uniform');
Tplot.t_dFF_DRCoff = rowfun(@(t) ...
    {t(find((t>=tBaseDRCoff(1) & t<=tBaseDRCoff(2)),1,'first'):end)},...
    Tplot,'InputVariables',{'t_F'},...
    'ExtractCellContents',true,'OutputFormat','uniform');
tDRCoff = Tplot.t_dFF_DRCoff{1};

[G,idA,idR,idC] = findgroups(Tplot.animal,Tplot.roiID,Tplot.dBdelta);
tmp = cellfun(@(c) nanmean(c,1),splitapply(@(x) {vertcat(x{:})},Tplot.dFF_DRCoff,G),'uni',0);
[G2,idC2] = findgroups(idC);
cellTrace = splitapply(@(x) {cell2mat(x)},tmp,G2);
cellU = cell2mat(cellfun(@(c) nanmean(c,1),cellTrace,'uni',0)).*100;
cellSEM = cell2mat(cellfun(@(c) SEMcalc(c,1),cellTrace,'uni',0)).*100;

figure('Name',[params.cohort '_DRCoffset_traces_re_contrast_wSignal']);
subplot(2,1,1)
fillSEMplot(tDRCoff,cellU(idC2==10,:),...
    cellSEM(idC2==10,:),params.colors.lohiPre(1,:),...
    params.colors.lohiTracePre(1,:));
hold on
fillSEMplot(tDRCoff,cellU(idC2==30,:),...
    cellSEM(idC2==30,:),params.colors.lohiPre(2,:),...
    params.colors.lohiTracePre(2,:));
modPlotForPoster(0)
xlabel('time (s)')
ylabel('% \DeltaF/F')
xlim([7 10])
xl = xlim;

subplot(2,1,2)
 %signal
[sigPTcontrastHi,fs] = inspectSignalObject('plot',false,'signalPath',...
        'C:\Data\Rig Software\250kHzPulses\PC_PTinContrast_5-52kHz_25msDRC_10000Hz70dB_0s_delay\25msDRC_5-52kHz_50-60dB_8s_sin2-cos2_DRCramp_250kHz_stim_10000Hz_70dB_400ms_at_2s_3.signal');
    sigPTcontrastLo = inspectSignalObject('plot',false,'signalPath',...
        'C:\Data\Rig Software\250kHzPulses\PC_PTinContrast_5-52kHz_25msDRC_10000Hz70dB_0s_delay\25msDRC_5-52kHz_40-70dB_8s_sin2-cos2_DRCramp_250kHz_stim_10000Hz_70dB_400ms_at_2s_1.signal');
tSig = [1:length(sigPTcontrastLo)].*(1/fs);
plot(tSig,sigPTcontrastLo,'Color',params.colors.(...
    strjoin(['lohi' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:));
hold on
plot(tSig,sigPTcontrastHi+15000,'Color',params.colors.(...
    strjoin(['lohi' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:));
xlim(xl)
hLine(1) = yline(15000);
hLine(2) = yline(0);
hLine(1).Color = params.colors.(...
    strjoin(['lohi' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:);
hLine(2).Color = params.colors.(...
    strjoin(['lohi' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:);
set(gca,'yticklabel',[]);
modPlotForPoster(0);
xlabel('time (s)');
if figSave
    figSaveAsFigEpsPng(gcf);
end

%over epoch
cellTraceBin = cellfun(@(c) nanmean(c(:,tDRCoff>=drcOFFtBin(1) & ...
    tDRCoff<=drcOFFtBin(2)),2),cellTrace,'uni',0);
offBinU = cellfun(@nanmean,cellTraceBin,'uni',1).*100;
offBinSEM = cellfun(@SEMcalc,cellTraceBin,'uni',1).*100;

lo = cellTraceBin{idC2==10};
hi = cellTraceBin{idC2==30};
[hEpoch,pEpoch,normBool] = sigDiffCalc(lo,hi);
fprintf('sig offset diff re contrast? h = %d \n',...
    hEpoch);

%perm test
if normBool~=1
    try
        pPerm = permutationTest(lo,hi,params.permIters,'exact',1);
        iters = 'exact';
    catch
        pPerm = permutationTest(lo,hi,params.permIters);
        iters = num2str(params.permIters);
    end
    fprintf('PermTest:  OFFSET diff re contrast: h = %d | p = %d | iters: %s \n',pPerm<0.05,pPerm,iters);
end

%boxplot
figure('Name',[params.cohort '_DRCoffset_boxplot_re_' treatment]);
tmpG = arrayfun(@(r) ones(length(cellTraceBin{r}),1).*r,[1 2],'uni',0);
boxplot(vertcat(cellTraceBin{:}).*100,vertcat(tmpG{:}),...
    'PlotStyle','compact','Colors',params.colors.ratio(1,:))
modPlotForPoster(0)
ylabel('mean % \DeltaF/F')

a = get(get(gca,'children'),'children');   % Get the handles of all the objects
set(a(9),'Color',params.colors.lohiPre(2,:))
set(a(10),'Color',params.colors.lohiPre(1,:))
lowStr = strjoin(arrayfun(@(x) num2str(x),params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:),'UniformOutput',false),',');
highStr = strjoin(arrayfun(@(x) num2str(x),params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:),'UniformOutput',false),',');
set(gca,'xtick',1:2,...
        'xticklabel',{['\color[rgb]{' lowStr '} LOW'],...
    ['\color[rgb]{' highStr '} HIGH']})
set(gcf,'Position',[150 244.5000 369 420])
if figSave
    figSaveAsFigEpsPng(gcf);
end



%}

%% PN | PRE: DRC contrast change dFF (supplemental)
%{
clearvars;close all hidden;clc;
load('D:\Data\CaMKII_dataTable_contrastChange.mat')
load('C:\Data\CaMKII_params.mat')
params.permIters = 100000;
treatment = 'pre';
figSave = false;
tBaseDRC = [7 10];
tDRCchangeEpoch = [11.4 12];
tDRCchangeBar = [10.2 12];

Tinput.t_F = rowfun(@(F,fr,trigDelay) ...
    {([1:size(F,2)].*(1/fr))-trigDelay},...
    Tinput,'InputVariables',{'F','frameRate','trigDelay'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

%dFF re DRC
Tinput.dFF_DRC = rowfun(@(F,t) ...
    {dFoFcalc(F,[find((t>=tBaseDRC(1) & t<=tBaseDRC(2)),1,'first')...
    find((t>=tBaseDRC(1) & t<=tBaseDRC(2)),1,'last')],1)},...
    Tinput,'InputVariables',{'F','t_F'},...
    'ExtractCellContents',true,'OutputFormat','uniform');
Tinput.t_dFF_DRC = rowfun(@(t) ...
    {t(find((t>=tBaseDRC(1) & t<=tBaseDRC(2)),1,'first'):end)},...
    Tinput,'InputVariables',{'t_F'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

tDFFdrc = Tinput.t_dFF_DRC{1};

[G,~,~,~] = findgroups(Tinput.animal,Tinput.roiID,Tinput.stimID);
if ~(max(G)==length(G))
    error('compute cell average first')
else
    [G,idC] = findgroups(Tinput.stimID);
    uCell_trace = splitapply(@(x) {cell2mat(x)},cellfun(@(c) nanmean(c,1),Tinput.dFF_DRC,'uni',0),G);
end

u = cell2mat(cellfun(@(c) nanmean(c,1),uCell_trace,'uni',0)).*100;
SEM = cell2mat(cellfun(@(c) SEMcalc(c,1),uCell_trace,'uni',0)).*100;

tDFFdrc = Tinput.t_dFF_DRC{1};
figure('Name','CaMKII_contrastChange_dFFtrace_re_pre');
subplot(2,1,1)
fillSEMplot(tDFFdrc,u(1,:),SEM(1,:),...
    params.colors.lohiPre(2,:),params.colors.lohiTracePre(2,:));
hold on
fillSEMplot(tDFFdrc,u(2,:),SEM(2,:),...
    params.colors.lohiPre(1,:),params.colors.lohiTracePre(1,:));
modPlotForPoster(0)
xlabel('time (s)')
ylabel('% \DeltaF/F')
xlim([9 14])
xL(1) = xline(10);
xL(1).LineWidth = 2;
hL = legend({'LOW \rightarrow HIGH','HIGH \rightarrow LOW','contrast change'});
hL.Box = 'off';
hL.Position = [0.6202 0.8069 0.2960 0.1692];

subplot(2,1,2) %signal
[sigLoHi,fs] = inspectSignalObject('plot',false,'signalPath',...
    'C:\Data\Rig Software\250kHzPulses\PC_contrastChange_25msDRC_5-52kHz_50-60_40-70dB_10sEach\25msDRC_5-52kHz_50-60dB_to_40-70dB_10s_sin2-cos2_DRCramp_250kHz_2.signal');
sigHiLo = inspectSignalObject('plot',false,'signalPath',...
    'C:\Data\Rig Software\250kHzPulses\PC_contrastChange_25msDRC_5-52kHz_50-60_40-70dB_10sEach\25msDRC_5-52kHz_40-70dB_to_50-60dB_10s_sin2-cos2_DRCramp_250kHz_2.signal');
tSig = [1:length(sigLoHi)].*(1/fs);
subplot(2,1,2)
plot(tSig,sigHiLo,'Color',params.colors.lohiPre(1,:));
hold on
plot(tSig,sigLoHi+20000,'Color',params.colors.lohiPre(2,:));
set(gca,'yticklabel',[])
ax = gca;
ax.YAxis.TickLength = [0 0];
xL(2) = xline(10);
xL(2).LineWidth = 2;
modPlotForPoster(0);
xlim([9 14])
xlabel('time (s)')
set(gcf,'Position',[287 302.5000 674 455])
if figSave
    figSaveAsFigEpsPng(gcf);
end

% 11.6-12 scatter
uDRCchangeEpoch = cellfun(@(c) nanmean(c(:,tDFFdrc>=tDRCchangeEpoch(1) & ...
    tDFFdrc<=tDRCchangeEpoch(2)).*100,2)',uCell_trace,'uni',0);
lo = uDRCchangeEpoch{2}';
hi = uDRCchangeEpoch{1}';

%stats
[hEpoch,pEpoch,normBool] = sigDiffCalc(uDRCchangeEpoch{:});
fprintf('sig diff epoch? h = %d | p = %d \n',hEpoch,pEpoch)
%perm test
if normBool~=1
    try
        pPerm = permutationTest(uDRCchangeEpoch{:},params.permIters,'exact',1);
        iters = 'exact';
    catch
        pPerm = permutationTest(uDRCchangeEpoch{:},params.permIters);
        iters = num2str(params.permIters);
    end
    fprintf('PermTest:  epoch diff re contrast: h = %d | p = %d | iters: %s \n',pPerm<0.05,pPerm,iters);
end

%boxplot
figure('Name',sprintf('CaMKII_contrastChange_boxplot_epoch_%.1f-%.1fs_re_pre',tDRCchangeEpoch));
tmpG = arrayfun(@(r) ones(length(uDRCchangeEpoch{r}),1).*r,[1 2],'uni',0);
boxplot(vertcat(uDRCchangeEpoch{:}),vertcat(tmpG{:}),...
    'PlotStyle','compact','Colors',params.colors.ratio(1,:))
h=gca; h.XAxis.TickLength = [0 0];
modPlotForPoster(0)
ylabel('mean % \DeltaF/F')
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
set(a(9),'Color',params.colors.lohiPre(2,:))
set(a(10),'Color',params.colors.lohiPre(1,:))
lowStr = strjoin(arrayfun(@(x) num2str(x),params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:),'UniformOutput',false),',');
highStr = strjoin(arrayfun(@(x) num2str(x),params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:),'UniformOutput',false),',');
set(gca,'xtick',1:2,...
        'xticklabel',{['\color[rgb]{' lowStr '} LOW'],...
    ['\color[rgb]{' highStr '} HIGH']})
set(gcf,'Position',[150 244.5000 369 420])
if figSave
    figSaveAsFigEpsPng(gcf);
end

%}

%% PN | PRE --> POST: ZX1 re tuning (supplemental)
%{
clearvars;close all hidden;clc;
load('D:\Data\CaMKII_dataTable_prePostMap.mat')
load('D:\Data\CaMKII_params_prePostMap.mat')
params.permIters = 100000;

figSave = false;
DRUG = 'ZX1';
PTonset = 2;
% dFF
tBaseDRC = [-1.2 0];
tBasePT{1} = [0.2 0.5];
tBasePT{2} = [1 2];

Tinput.t_F = rowfun(@(F,fr,trigDelay) ...
    {([1:size(F,2)].*(1/fr))-trigDelay},...
    Tinput,'InputVariables',{'F','frameRate','trigDelay'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

%dFF re DRC
Tinput.dFF_DRC = rowfun(@(F,t) ...
    {dFoFcalc(F,[find((t>=tBaseDRC(1) & t<=tBaseDRC(2)),1,'first')...
    find((t>=tBaseDRC(1) & t<=tBaseDRC(2)),1,'last')],1)},...
    Tinput,'InputVariables',{'F','t_F'},...
    'ExtractCellContents',true,'OutputFormat','uniform');
Tinput.t_dFF_DRC = rowfun(@(t) ...
    {t(find((t>=tBaseDRC(1) & t<=tBaseDRC(2)),1,'first'):end)},...
    Tinput,'InputVariables',{'t_F'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

%dFF re PT
Tinput.dFF_PT = rowfun(@(F,t,PTonset) ...
    {dFoFcalc(F,[find((t>=tBasePT{round(PTonset)}(1) & t<=tBasePT{round(PTonset)}(2)),1,'first')...
    find((t>=tBasePT{round(PTonset)}(1) & t<=tBasePT{round(PTonset)}(2)),1,'last')],1)},...
    Tinput,'InputVariables',{'F','t_F','PTsOnset'},...
    'ExtractCellContents',true,'OutputFormat','uniform');
Tinput.t_dFF_PT = rowfun(@(t,PTonset) ...
    {t(find((t>=tBasePT{round(PTonset)}(1) & t<=tBasePT{round(PTonset)}(2)),1,'first'):end)},...
    Tinput,'InputVariables',{'t_F','PTsOnset'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

% %peak dFF @PT (sig)
%looks 0.2 s after PT onset, consider 3 frames instead
tmp = rowfun(@(dFF_PT,t_dFF_PT,PTonset,fr) ...
    pkFcalc(dFF_PT,...
    find(t_dFF_PT>=(PTonset+(1/fr)),1,'first'),...
    params.pkPTframeBin,params.pkPTsigSD),...
    Tinput,'InputVariables',{'dFF_PT','t_dFF_PT','PTsOnset','frameRate'},...
    'ExtractCellContents',true,'OutputFormat','cell','OutputVariableNames',{'pk','sigPk','sig'});
Tinput.pkPT_sig = tmp(:,1);
Tinput.sigPk = tmp(:,2);
clear tmp

Tplot = Tinput(Tinput.PTsOnset==PTonset & ...
    contains(Tinput.treatment,DRUG) & ...
    Tinput.dBFprepost<=0.25,:); 

[G,~,~,idT,idC] = findgroups(Tplot.animal,Tplot.roiID,Tplot.treatment,Tplot.dBdelta);
if ~(max(G)==length(G))
        error('compute cell average first')
else
    [G,gTreatment,gContrast] = findgroups(Tplot.treatment,Tplot.dBdelta);
    uPkPT = splitapply(@(x) {cellfun(@nanmean,x)},Tplot.pkPT_sig,G);
end

tmpT = table(gTreatment,gContrast,uPkPT);
tmpT.n = cellfun(@(c) sum(~isnan(c)),tmpT.uPkPT,'uni',1);

%RATIO pre & post bar plot
uPkByTreatment = splitapply(@(x) {x},tmpT.uPkPT,findgroups(tmpT.gTreatment));
pkRespT = table(uPkByTreatment,...
    unique(tmpT.gTreatment,'stable'),...
    repmat({unique(tmpT.gContrast,'stable')},max(findgroups(tmpT.gTreatment)),1),...
    'VariableNames',{'uPkPT','treatment','contrast'});
pkRespT.pkRatio = arrayfun(@(r) r.uPkPT{r.contrast==min(r.contrast)}./...
    r.uPkPT{r.contrast==max(r.contrast)},...
    table2struct(pkRespT),'uni',0);
%sort pre/post (opposite alphabetical)
[~,tmp] = sort(pkRespT.treatment);
pkRespT = pkRespT(flipud(tmp),:);
u = cellfun(@(c) nanmean(c),pkRespT.pkRatio,'uni',1);
sem = cellfun(@(c) SEMcalc(c),pkRespT.pkRatio,'uni',1);

figure('Name',['CaMKII_RATIO_re_quarterOctave_re_' DRUG]);
b = bar(u,'FaceColor','Flat','EdgeColor','none');
b.CData = params.colors.ratio;
hold on
errorbar(u,sem,...
    'k-','LineWidth',1.5,'LineStyle','none')
hLabel = ylabelLowHighPkRespRatio(params.colors);
modPlotForPoster(1)
xticklabels('')
if figSave
    figSaveAsFigEpsPng(gcf);
end

fprintf('ratio n = %d \n',sum(all(~isnan(horzcat(pkRespT.pkRatio{:})),2)))
if adtest(vertcat(pkRespT.pkRatio{:}))
    disp(['ratio pre/post p = ' num2str(signrank(pkRespT.pkRatio{:}))])
    disp(['ratio pre p = ' num2str(signrank(pkRespT.pkRatio{1},1))])
    disp(['ratio post p = ' num2str(signrank(pkRespT.pkRatio{2},1))])
end
%perm test
if adtest(vertcat(pkRespT.pkRatio{:}))
    try
        pPerm = permutationTest(pkRespT.pkRatio{:},params.permIters,'exact',1);
        iters = 'exact';
    catch
        pPerm = permutationTest(pkRespT.pkRatio{:},params.permIters);
        iters = num2str(params.permIters);
    end
    fprintf('PermTest:  ratio diff re drug: h = %d | p = %d | iters: %s \n',pPerm<0.05,pPerm,iters);
end

figure('Name',['CaMKII_BFshift_re_' DRUG]);
dBF = Tinput{contains(Tinput.treatment,'pre'),{'dBFprepost'}};
histogram(dBF,'EdgeColor','none','FaceColor',params.colors.ratio(1,:),'HandleVisibility','off')
modPlotForPoster(0)
xlabel({'BF octave shift from control','upon ZX1 injection'})
ylabel('cell count')
hold on
uLine = xline(nanmean(dBF));
uLine.LineWidth = 2;
uLine.Color = 'k';
SEMline(1) = xline(nanmean(dBF) + SEMcalc(dBF));
SEMline(2) = xline(nanmean(dBF) - SEMcalc(dBF));
SEMline(2).HandleVisibility = 'off';
for i = 1:length(SEMline)
    SEMline(i).LineWidth = 1;
    SEMline(i).LineStyle = '--';
    SEMline(i).Color = 'k';
end
hL = legend('mean','SEM');
hL.Box = 'off';
if figSave
    figSaveAsFigEpsPng(gcf);
end

if adtest(dBF)
    pNorm = 0;
    [pDBF,hDBF] = signrank(dBF,0);
else
    pNorm = 1;
    [hDBF,pDBF,~,t_stats] = ttest(dBF,0);
end
fprintf('sig octave shift pre-->post? h = %d | p = %d | norm? %d \n',hDBF,pDBF,pNorm)

%}

%% ZnT3?? | PRE --> POST:  dFF_PT RATIO
%{
close all hidden;clc;clearvars -except Tinput params

figSave = false;
DRUG = 'ZX1';
PTonset = 2;

Tplot = Tinput(Tinput.PTsOnset==PTonset & ...
    contains(Tinput.treatment,DRUG),:);

[G,~,~,idT,idC] = findgroups(Tplot.animal,Tplot.roiID,Tplot.treatment,Tplot.dBdelta);
if ~(max(G)==length(G))
    tmp = splitapply(@(x) {nanmean(vertcat(x{:}),1)},Tplot.pkPT_sig,G); %cell average pkPT
    [G,idC,idT] = findgroups(idC,idT);
    uPkPT = splitapply(@(x) {cell2mat(x)},tmp,G);
    %     error('compute cell average first')
else
    [G,idT,idC] = findgroups(Tplot.treatment,Tplot.dBdelta);
    uPkPT = splitapply(@(x) {cellfun(@nanmean,x)},Tplot.pkPT_sig,G);
end

tmpT = table(idT,idC,uPkPT);
tmpT.n = cellfun(@(c) sum(~isnan(c)),tmpT.uPkPT,'uni',1);

%RATIO pre & post bar plot
uPkByTreatment = splitapply(@(x) {x},tmpT.uPkPT,findgroups(tmpT.idT));
pkRespT = table(uPkByTreatment,...
    unique(tmpT.idT,'stable'),...
    repmat({unique(tmpT.idC,'stable')},max(findgroups(tmpT.idT)),1),...
    'VariableNames',{'uPkPT','treatment','contrast'});
pkRespT.pkRatio = arrayfun(@(r) r.uPkPT{r.contrast==min(r.contrast)}./...
    r.uPkPT{r.contrast==max(r.contrast)},...
    table2struct(pkRespT),'uni',0);
%sort pre/post (opposite alphabetical)
[~,tmp] = sort(pkRespT.treatment);
pkRespT = pkRespT(flipud(tmp),:);
u = cellfun(@(c) nanmean(c),pkRespT.pkRatio,'uni',1);
sem = cellfun(@(c) SEMcalc(c),pkRespT.pkRatio,'uni',1);

figure('Name',[params.cohort '_RATIO_re_' DRUG]);
b = bar(u,'FaceColor','Flat','EdgeColor','none');
b.CData = params.colors.ratio;
hold on
errorbar(u,sem,...
    'k-','LineWidth',1.5,'LineStyle','none')
hLabel = ylabelLowHighPkRespRatio(params.colors);
modPlotForPoster(1)
xticklabels('')
if figSave
    figSaveAsFigEpsPng(gcf);
end

fprintf('ratio n = %d \n',sum(all(~isnan(horzcat(pkRespT.pkRatio{:})),2)))
if adtest(vertcat(pkRespT.pkRatio{:}))
    disp(['ratio pre/post p = ' num2str(signrank(pkRespT.pkRatio{:}))])
    disp(['ratio pre p = ' num2str(signrank(pkRespT.pkRatio{1},1))])
    disp(['ratio post p = ' num2str(signrank(pkRespT.pkRatio{2},1))])
end

%perm test
if adtest(vertcat(pkRespT.pkRatio{:}))
    try
        pPerm = permutationTest(pkRespT.pkRatio{:},params.permIters,'exact',1);
        iters = 'exact';
    catch
        pPerm = permutationTest(pkRespT.pkRatio{:},params.permIters);
        iters = num2str(params.permIters);
    end
    fprintf('PermTest:  ratio diff re drug: h = %d | p = %d | iters: %s \n',pPerm<0.05,pPerm,iters);
end

%}

%% PV & SOM | PRE OR POST: dFF_PT | scatterplot and traces re contrast
%{
close all hidden;clc;clearvars -except Tinput params

figSave = false;
% treatment = 'pre';
treatment = 'postZX1';
PTonset = 2;

Tplot = Tinput(contains(Tinput.treatment,treatment) & ...
    Tinput.PTsOnset==PTonset,:);

[G,idC] = findgroups(Tplot.dBdelta);

%cell avg pk resp re contrast
uPkCell = splitapply(@(x) {cellfun(@nanmean,x,'uni',1)},Tplot.pkPT_sig,G);

if logical(sigDiffCalc(uPkCell{idC==10},uPkCell{idC==30}))   
    regPlotColors(1,:) = params.colors.lohiPre((~(nanmean(uPkCell{idC==10})>...
        nanmean(uPkCell{idC==30})))+1,:);
    regPlotColors(2,:) = params.colors.lohiTracePre((~(nanmean(uPkCell{idC==10})>...
        nanmean(uPkCell{idC==30})))+1,:);
else
    regPlotColors(1,:) = [0,0,0];
    regPlotColors(2,:) = [0.6510    0.6510    0.6510];
end
lo = uPkCell{idC==10};
hi = uPkCell{idC==30};
%pk resp scatter
[hScatter, hRegLine, hRegLineCI, mdl, ciRegLine, hRef] = ...
    regPlot(lo.*100,...
    hi.*100,...
    'regLine',true,...
    'logScale',false,...
    'intercept',false,...
    'colors',regPlotColors,...
    'fName',[params.cohort '_PTresp_scatter_re_contrast_re_' treatment]);
xlabel({'low contrast','peak % \DeltaF/F'},...
    'Color',params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:),'interpreter', 'tex')
ylabel({'high contrast','peak % \DeltaF/F'},...
    'Color',params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:),'interpreter', 'tex')
legend('off')
modPlotForPoster(0);
if figSave
    figSaveAsFigEpsPng(gcf);
end
fprintf('lo>hi %d | hi>lo %d | total %d \n',...
    sum((lo>hi)),sum((hi>lo)),sum(all(~isnan([lo hi]),2)));
[hLOHI,pLOHI,normBool] = sigDiffCalc(lo,hi);
fprintf('lo v. hi sig diff? h = %d | p = %d \n',hLOHI,pLOHI)

%perm test
if normBool~=1
    try
        pPerm = permutationTest(lo,hi,params.permIters,'exact',1);
        iters = 'exact';
    catch
        pPerm = permutationTest(lo,hi,params.permIters);
        iters = num2str(params.permIters);
    end
    fprintf('PermTest:  PTonset diff re contrast: h = %d | p = %d | iters: %s \n',pPerm<0.05,pPerm,iters);
end

mCI = mdl.coefCI;
fprintf('does slope confidence interval (%.2f - %.2f) include 1? %d \n',mdl.coefCI,(1>mCI(1) & 1<mCI(2)))

%pk resp trace
sigTrace = rowfun(@(dFF_PT,sig) ...
    {dFF_PT(sig,:)},Tplot,...
    'InputVariables',{'dFF_PT','sigPk'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

uTrace = cell2mat(splitapply(@(x) ...
    {nanmean(cell2mat(cellfun(@(c) nanmean(c,1),x,'uni',0)),1).*100},sigTrace,G));
semTrace = cell2mat(splitapply(@(x) ...
    {SEMcalc(cell2mat(cellfun(@(c) nanmean(c,1),x,'uni',0)),1).*100},sigTrace,G));
tDFF = Tplot.t_dFF_PT{1};

figure('Name',[params.cohort '_PTresp_traces_re_contrast_re_' treatment]);
fillSEMplot(tDFF,uTrace(idC==10,:),semTrace(idC==10,:),...
    params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(idC==10,:),...
    params.colors.(strjoin(['lohiTrace' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(idC==10,:));
hold on
fillSEMplot(tDFF,uTrace(idC==30,:),semTrace(idC==30,:),...
    params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(idC==30,:),...
    params.colors.(strjoin(['lohiTrace' strrep(regexp(treatment,...
    strjoin({'pre','post'},'|'),'match'),'p','P')],''))(idC==30,:));
xlim([1.5 3.5])
modPlotForPoster(0)
xlabel('time (s)')
ylabel('% \DeltaF/F')
if figSave
    figSaveAsFigEpsPng(gcf);
end

%}

%% PV & SOM | PRE OR POST:  dFF_DRC
%{
close all hidden;clc;clearvars -except Tinput params
figSave = false;
tSustainedDRC = [6 8];
PTonset = 2;
treatment = 'postZX1';
% treatment = 'pre';
plotBoxplot = true;

if plotBoxplot==1
    Tplot = Tinput(contains(Tinput.treatment,treatment),:); %for boxplot
else
    Tplot = Tinput(contains(Tinput.treatment,treatment) & ...
        Tinput.PTsOnset==PTonset,:); %for traces
end

[G,~,~,idC] = findgroups(Tplot.animal,Tplot.roiID,Tplot.dBdelta);
if ~(max(G)==length(G))
    tmp = splitapply(@(x) {nanmean(vertcat(x{:}),1)},Tplot.dFF_DRC,G);
    [G,idC] = findgroups(idC);
    uDFF_DRC_cell = splitapply(@(x) {cell2mat(x)},tmp,G);
%     error('compute cell average first')
else
    [G,idC] = findgroups(Tplot.dBdelta);
    uDFF_DRC_cell = splitapply(@(x) {cell2mat(x)},cellfun(@(c) nanmean(c,1),Tplot.dFF_DRC,'uni',0),G);
end

uDFF_DRC_cell_u = cell2mat(cellfun(@(c) nanmean(c,1).*100,uDFF_DRC_cell,'uni',0));
uDFF_DRC_cell_sem = cell2mat(cellfun(@(c) SEMcalc(c,1).*100,uDFF_DRC_cell,'uni',0));

tDRC = Tplot.t_dFF_DRC{1};

if plotBoxplot ~= 1
    figure('Name',[params.cohort '_sustainedDRC_traces_re_' treatment]);
    subplot(2,1,1) %dFF
    fillSEMplot(tDRC,...
        uDFF_DRC_cell_u(idC==10,:),...
        uDFF_DRC_cell_sem(idC==10,:),...
        params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:),...
        params.colors.(...
        strjoin(['lohiTrace' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:));
    hold on
    fillSEMplot(tDRC,...
        uDFF_DRC_cell_u(idC==30,:),...
        uDFF_DRC_cell_sem(idC==30,:),...
        params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:),...
        params.colors.(...
        strjoin(['lohiTrace' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:));
    ylabel('% \DeltaF/F')
    xlim([-1.2 8])
    modPlotForPoster(0)
    
    subplot(2,1,2) %signal
    [sigPTcontrastHi,fs] = inspectSignalObject('plot',false,'signalPath',...
        'C:\Data\Rig Software\250kHzPulses\PC_PTinContrast_5-52kHz_25msDRC_10000Hz70dB_0s_delay\25msDRC_5-52kHz_50-60dB_8s_sin2-cos2_DRCramp_250kHz_stim_10000Hz_70dB_400ms_at_2s_2.signal');
    sigPTcontrastLo = inspectSignalObject('plot',false,'signalPath',...
        'C:\Data\Rig Software\250kHzPulses\PC_PTinContrast_5-52kHz_25msDRC_10000Hz70dB_0s_delay\25msDRC_5-52kHz_40-70dB_8s_sin2-cos2_DRCramp_250kHz_stim_10000Hz_70dB_400ms_at_2s_3.signal');
    tSig = [1:length(sigPTcontrastLo)].*(1/fs);
    plot(tSig,sigPTcontrastLo,'Color',params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:));
    hold on
    plot(tSig,sigPTcontrastHi+15000,'Color',params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:));
    xlim([-1.2 8])
    hLine(1) = yline(15000);
    hLine(2) = yline(0);
    hLine(1).Color = params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:);
    hLine(2).Color = params.colors.(...
        strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:);
    set(gca,'yticklabel',[]);
    modPlotForPoster(0);
    xlabel('time (s)');
    ax = gca;
    ax.YAxis.TickLength = [0 0];
    if figSave
        figSaveAsFigEpsPng(gcf);
    end
end

if plotBoxplot == 1
    
    %sustained DRC
    uDFF_DRCsustained_cell = cellfun(@(c) nanmean(c(:,tDRC>=tSustainedDRC(1) & ...
        tDRC<=tSustainedDRC(2)),2),uDFF_DRC_cell,'uni',0);
    lo = uDFF_DRCsustained_cell{idC==10};
    hi = uDFF_DRCsustained_cell{idC==30};
    
    [hDRC,pDRC,normBool] = sigDiffCalc(lo,hi);
    fprintf('sustained DRC diff re contrast:  h = %d | p = %d \n',hDRC,pDRC)
    
    %perm test
    if normBool~=1
        try
            pPerm = permutationTest(lo,hi,params.permIters,'exact',1);
            iters = 'exact';
        catch
            pPerm = permutationTest(lo,hi,params.permIters);
            iters = num2str(params.permIters);
        end
        fprintf('PermTest:  sustained DRC diff re contrast: h = %d | p = %d | iters: %s \n',pPerm<0.05,pPerm,iters);
    end
    
    %boxplot
    figure('Name',[params.cohort '_sustainedDRC_boxplot_re_' treatment]);
    tmpG = arrayfun(@(r) ones(length(uDFF_DRCsustained_cell{r}),1).*r,[1 2],'uni',0);
    boxplot(vertcat(uDFF_DRCsustained_cell{:}).*100,vertcat(tmpG{:}),...
        'PlotStyle','compact','Colors',params.colors.ratio(1,:))
    modPlotForPoster(0)
    ylabel('mean % \DeltaF/F')
    h=gca; h.XAxis.TickLength = [0 0];
    a = get(get(gca,'children'),'children');   % Get the handles of all the objects
    set(a(9),'Color',params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:))
    set(a(10),'Color',params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:))
    lowStr = strjoin(arrayfun(@(x) num2str(x),params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(1,:),'UniformOutput',false),',');
    highStr = strjoin(arrayfun(@(x) num2str(x),params.colors.(strjoin(['lohi' strrep(regexp(treatment,...
        strjoin({'pre','post'},'|'),'match'),'p','P')],''))(2,:),'UniformOutput',false),',');
    set(gca,'xtick',1:2,...
        'xticklabel',{['\color[rgb]{' lowStr '} LOW'],...
        ['\color[rgb]{' highStr '} HIGH']})
    set(gcf,'Position',[150 244.5000 369 420])
    if figSave
        figSaveAsFigEpsPng(gcf);
    end
end

%}