function hLabel = ylabelLowHighPkRespRatio(colors,varargin)

ylabel('')

%Positions of elements
posNum = [-0.27 0.5 0.5];
posBar = [-0.235 0.5 0.5];
posDen = [-0.195 0.5 0.5];

fontSize = 16;

if isempty(varargin)
    signal = 'F/F';
    barLine = '$$\rule{6cm}{0.4mm}$$';
else
    signal = varargin{1};
    barLine = '$$\rule{7cm}{0.4mm}$$';
end

num = ['peak low contrast \Delta' signal];
den = ['peak high contrast \Delta' signal];

hLabel.tNum = text(1,1,num,'Units','Normalized','FontSize',fontSize);
hLabel.tBar = text(1,1,barLine,'Interpreter','latex','Units','Normalized');
hLabel.tDen = text(1,1,den,'Units','Normalized','FontSize',fontSize);

set(hLabel.tNum,'Position',posNum,'Rotation',90,'Color',...
    colors.lohiPre(1,:),'HorizontalAlignment','center')
set(hLabel.tBar,'Position',posBar,'Rotation',90,'HorizontalAlignment','center')
set(hLabel.tDen,'Position',posDen,'Rotation',90,'Color',...
    colors.lohiPre(2,:),'HorizontalAlignment','center')
set(gca,'Position',[0.2000 0.18 0.6000 0.7150])

end
