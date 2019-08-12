%Analyze vtf-files
%Read the file line-by-line
clear all

working_dir=cd;
datadir='D:\Dropbox\MD\Esp\jmannik';
Np1=150;
Np2=150;
L=25;
Nbin=25;    %Total lenght of histogram is 2Nbin+1; It spans Nbin in both directions from zero

Nexclude=200;
Navg=1000;

cd(datadir);
[file,path]=uigetfile('*.vtf','Select vtf file to be analyzed');
file_name=[path,file];

fid=fopen(file_name,'r');

istimestep=0;       %Find 1st timestep
while ~istimestep
    tline = fgetl(fid);
    istimestep=strncmp(tline,'timestep',8);
end;

%Remove first Nexclude time steps
for i=1:Nexclude
    istimestep=0;
    while ~istimestep
        tline = fgetl(fid);
        istimestep=strncmp(tline,'timestep',8);
    end;
end;

%Initialize histogram array
histY1=double(zeros(2*Nbin+1));
histY2=double(zeros(2*Nbin+1));
histX=double(zeros(2*Nbin+1));
dL=double(L)/double(Nbin);
x=-L;
for i=1:2*Nbin+1
    histX(i)=x;
    x=x+dL;
end;

for i=1:Navg
    for j=1:Np1
       tline = fgetl(fid);
       Z=sscanf(tline,'%f%f%f');
       indx=floor((Z(3)-25)/dL)+Nbin+1;
       histY1(indx)=histY1(indx)+1;
    end;
    for j=1:Np2
       tline = fgetl(fid);
       Z=sscanf(tline,'%f%f%f');
       indx=floor((Z(3)-25)/dL)+Nbin+1;
       histY2(indx)=histY2(indx)+1;
    end;
    tline = fgetl(fid);     %should read empty line
    tline = fgetl(fid);     %should read timestep
end;
    
figure(1);
plot(histX,histY1,'r');
hold on
plot(histX,histY2,'b');
plot(histX,histY1+histY2,'b');
hold off