function [imO,imO2d,imHistZ,imHistR,imHistX,imRg,imRg2d,imAR,imAR2d] = AnlyPlots2(Simul,Rtab,Anly)
% For one plolymer and discreate crowders
% DY180614
%%
LineWidth		= 2;
LabelFontSize	= 12;
Interpreter		= 'tex';
%%
L	= 2 * Rtab.L;
D	= 2 * Rtab.R;
ac2	= Rtab.ac2;
Nm  = Rtab.Nm;
Nc2	= Rtab.Nc2;
Lr  = 0.5 * L * Rtab.Lr - 0.5;
% OSh	= Rtab.ac * Rtab.Ro;
SimulNameSSP= strsplit(Simul.Name,'_');
SimulNameCell	= join(SimulNameSSP,'\_');
SimulName	= SimulNameCell{1};
TitleTexts	= {SimulName,['L = ',int2str(L),'    D = ',int2str(D),...
	'    ac = ',num2str(ac2),'    Nc = ',int2str(Nc2),...
	'    OutShell = ',num2str(Rtab.Ro,'%.1f'),' \times ac']};
%%
HistStart	= Anly.HistStart;
HistEnd		= Anly.HistEnd;
HistSet		= HistStart:HistEnd;
HistLength	= HistEnd - HistStart + 1;

RogStart	= Anly.RogStart;
RogEnd		= Anly.RogEnd;
RogSet		= RogStart:RogEnd;
% RogLength	= RogEnd - RogStart + 1;

% cx	= Anly.cx;
% cy	= Anly.cy;
% cz	= Anly.cz;
%% HistZ
ZEs = Simul.HistZEs;
Z	= mean([ZEs(1:end-1);ZEs(2:end)]);
Zm	= sum(Simul.HistZ(HistSet,:,1))./diff(ZEs)/HistLength;
figure(1);
% set(gcf,'Position',[360 228 960 720]);
if Nc2
	rcm = Nm/Nc2;
	Zc = rcm * sum(Simul.HistZ(HistSet,:,3))./diff(ZEs)/HistLength;
	plot(Z,Zm,'r',Z,Zc,'b','LineWidth',LineWidth);
	legend('DNA','Crowders','Location','best');
else
	plot(Z,Zm,'r','LineWidth',LineWidth);
	legend('DNA','Location','best');
end
hold on
title(TitleTexts,'Interpreter',Interpreter);
xlabel('Z','FontSize',LabelFontSize);
ylabel('Count density','FontSize',LabelFontSize);
hold off
Frm=getframe(gcf);
[imHistZ,~] = frame2im(Frm);

%% HistR
REs = Simul.HistREs;
R	= mean([REs(1:end-1);REs(2:end)]);
Rm	= sum(Simul.HistR(HistSet,:,1))./diff(REs)./R/HistLength;
figure(1);
% set(gcf,'Position',[360 228 960 720]);
if Nc2
	rcm = Nm/Nc2;
	Rc = rcm * sum(Simul.HistR(HistSet,:,3))./diff(REs)./R/HistLength;
	plot(R,Rm,'r',R,Rc,'b','LineWidth',LineWidth);
	legend('DNA','Crowders','Location','best');
else
	plot(R,Rm,'r','LineWidth',LineWidth);
	legend('DNA','Location','best');
end
hold on
title(TitleTexts,'Interpreter',Interpreter);
xlabel('R','FontSize',LabelFontSize);
ylabel('Count density','FontSize',LabelFontSize);
hold off
Frm=getframe(gcf);
[imHistR,~] = frame2im(Frm);
%% HistX
XEs = Simul.HistXEs;
X	= mean([XEs(1:end-1);XEs(2:end)]);
Xm	= sum(Simul.HistX(HistSet,:,1))./diff(XEs)/HistLength;
figure(1);
% set(gcf,'Position',[360 228 960 720]);
if Nc2
	rcm = Nm/Nc2;
	Xc = rcm * sum(Simul.HistX(HistSet,:,3))./diff(XEs)/HistLength;
	plot(X,Xm,'r',X,Xc,'b','LineWidth',LineWidth);
	legend('DNA','Crowders','Location','best');
else
	plot(X,Xm,'r','LineWidth',LineWidth);
	legend('DNA','Location','best');
end
hold on
title(TitleTexts,'Interpreter',Interpreter);
xlabel('X','FontSize',LabelFontSize);
ylabel('Count density','FontSize',LabelFontSize);
hold off
Frm=getframe(gcf);
[imHistX,~] = frame2im(Frm);

%% Rg 3D
figure(1);
plot(RogSet,Simul.Rg(RogSet),'r',RogSet,Simul.Re(RogSet),'g',...
	[RogStart,RogEnd],[Lr,Lr],'b');
hold on
title(TitleTexts,'Interpreter',Interpreter);
xlabel('Frame','FontSize',LabelFontSize);
ylabel('Linear size','FontSize',LabelFontSize);
legend('Rg','Re','CL/2 - 0.5','Location','best');
hold off
Frm=getframe(gcf);
[imRg,~] = frame2im(Frm);

%% Rg 2D
figure(1);
plot(RogSet,Simul.Rg2d(RogSet),'r',RogSet,Simul.Re(RogSet),'g',...
	[RogStart,RogEnd],[Lr,Lr],'b');
hold on
title(TitleTexts,'Interpreter',Interpreter);
xlabel('Frame','FontSize',LabelFontSize);
ylabel('Linear size','FontSize',LabelFontSize);
legend('Rg - 2D','Re','CL/2 - 0.5','Location','best');
hold off
Frm=getframe(gcf);
[imRg2d,~] = frame2im(Frm);

%% AR 3D
figure(1);
plot(RogSet,Simul.A(RogSet),'b',RogSet,Simul.Ap(RogSet),'r');
hold on
title(TitleTexts,'Interpreter',Interpreter);
xlabel('Frame','FontSize',LabelFontSize);
ylabel('Aspect ratio','FontSize',LabelFontSize);
legend('L/L','M/M','Location','best');
hold off
Frm=getframe(gcf);
[imAR,~] = frame2im(Frm);

%% AR 2D
figure(1);
plot(RogSet,Simul.A2d(RogSet),'b',RogSet,Simul.Ap(RogSet),'r');
hold on
title(TitleTexts,'Interpreter',Interpreter);
xlabel('Frame','FontSize',LabelFontSize);
ylabel('2D Aspect ratio','FontSize',LabelFontSize);
legend('L/L','M/M','Location','best');
hold off
Frm=getframe(gcf);
[imAR2d,~] = frame2im(Frm);
%%
for c = 3:-1:1
	imO(:,:,c) = [imHistZ(:,:,c),imHistR(:,:,c);imRg(:,:,c),imAR(:,:,c)];
end
for c = 3:-1:1
	imO2d(:,:,c) = [imHistZ(:,:,c),imHistX(:,:,c);imRg2d(:,:,c),imAR2d(:,:,c)];
end