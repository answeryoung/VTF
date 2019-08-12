function AnalyzeVTF2(file)
%Analyze vtf-files
%Read the file line-by-line
% clear all

working_dir=cd;
datadir='D:\Dropbox\MD\Esp\jmannik\';
Np1=150;
Np2=150;
L=100;
Nbin=50;    %Total lenght of histogram is 2Nbin+1; It spans Nbin in both directions from zero

Nexclude=200;
Navg=1000;

cd(datadir);
if exist('file','var')
	file_name = [datadir,file]
else
	[file,path]=uigetfile('*.vtf','Select vtf file to be analyzed');
	file_name=[path,file]
end
fid=fopen(file_name,'r');

istimestep=0;       %Find 1st timestep
LineNumber = 0;
while ~istimestep
    tline = fgetl(fid);
	LineNumber = LineNumber + 1;
    istimestep=strncmp(tline,'timestep',8);
end;

%Remove first Nexclude time steps
for i=1:Nexclude
    istimestep=0;
    while ~istimestep
        tline = fgetl(fid);
		LineNumber = LineNumber + 1;
        istimestep=strncmp(tline,'timestep',8);
    end;
end;

%Initialize histogram array
dL		= double(L)/double(Nbin);
histX	= -L:dL:L;
histY1	= zeros(1,2*Nbin+1);
histY2	= zeros(1,2*Nbin+1);

try
	for i=1:Navg
		for j=1:Np1
			tline = fgetl(fid);
			LineNumber = LineNumber + 1;
			Z=sscanf(tline,'%f%f%f');
			indx=floor((Z(1)-25)/dL)+Nbin+1;
			if indx < 2
				indx = 1;
			elseif indx > 2*Nbin
				indx = 2*Nbin+1;
			end
		   histY1(indx)=histY1(indx)+1;
		end;
		tline = fgetl(fid);     %should read empty line
		LineNumber = LineNumber + 1;
		tline = fgetl(fid);     %should read timestep
		LineNumber = LineNumber + 1;
		for j=1:Np2
			tline = fgetl(fid);
			LineNumber = LineNumber + 1;
			Z=sscanf(tline,'%f%f%f');
			indx=floor((Z(2)-25)/dL)+Nbin+1;
			if indx < 2
				indx = 1;
			elseif indx > 2*Nbin
				indx = 2*Nbin+1;
			end
			histY2(indx)=histY2(indx)+1;
		end;
		tline = fgetl(fid);     %should read empty line
		LineNumber = LineNumber + 1;
		tline = fgetl(fid);     %should read timestep
		LineNumber = LineNumber + 1;
	end;
catch
end
tText = file;
tText(tText=='_') = ' ';
figure(1);
plot(histX,histY1,'r');
hold on
plot(histX,histY2,'b');
plot(histX,histY1+histY2,'b');
title(tText(1:end-4));
hold off

F=getframe(gcf);
[im,~] = frame2im(F);
imwrite(im,[file(1:end-4),'.png'],'png');