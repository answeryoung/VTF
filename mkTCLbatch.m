NewIdx = 22;
BatchRange = 385:8064;
ReRand = 1;
template = 0;
Mfile = 'FRCenterDS06.mat';
load(Mfile);
%% read template
header	= [head(Rtab.Polymer{NewIdx}),'_'];
fid		= fopen([header,int2strN(template,4),'.tcl'],'r');
i = 1;
tline	= fgetl(fid);
B{i,1}	= tline;
while ischar(tline)
    i	= i+1;
    tline	= fgetl(fid);
    B{i,1}	= tline;
end
fclose(fid);
nL = i - 1;

for NewIdx = BatchRange
%% assign geometric parameters
B{5}	= ['set name "',header,int2strN(NewIdx,4),'"'];
B{8}	= ['set L	',num2str(Rtab{NewIdx,'L'})];
B{10}	= ['set R	',num2str(Rtab{NewIdx,'R'})];
B{14}	= ['set Nm	',num2str(Rtab{NewIdx,'Nm'})];
B{16}	= ['set a	',num2str(Rtab{NewIdx,'a'})];
B{17}	= ['set ac1	',num2str(Rtab{NewIdx,'ac1'})];
B{18}	= ['set ac2	',num2str(Rtab{NewIdx,'ac2'})];
B{19}	= ['set aW	',num2str(Rtab{NewIdx,'aW'})];
B{21}	= ['set Nc1	',num2str(Rtab{NewIdx,'Nc1'})];
B{22}	= ['set Nc2	',num2str(Rtab{NewIdx,'Nc2'})];

B{30}	= ['set Pr		',num2str(Rtab{NewIdx,'Pr'})];
B{31}	= ['set C_off	',num2str(Rtab{NewIdx,'C_off'})];

B{33}	= ['set Zoff	',num2str(Rtab{NewIdx,'Zoff'})];
B{34}	= ['set FC		',num2str(Rtab{NewIdx,'FC'})];
% B{35}	= ['set ConstrictionStep	',num2str(Rtab{NewIdx,'ConstrictionStep'})];
if ~isnan(Rtab{NewIdx,'PorR0'})
	B{37}	= ['set PorR0	',num2str(Rtab{NewIdx,'PorR0'})];
end
if ~isnan(Rtab{NewIdx,'Cr'})
	B{38}	= ['set Cr		',num2str(Rtab{NewIdx,'Cr'})];
end
if ~isnan(Rtab{NewIdx,'psr'})
	B{39}	= ['set psr		',num2str(Rtab{NewIdx,'psr'})];
end
B{40}	= ['set comp_tau	',num2str(Rtab{NewIdx,'comp_tau'})];
B{41}	= ['set comp_alpha	',num2str(Rtab{NewIdx,'comp_alpha'})];
if ~isnan(Rtab{NewIdx,'Ro'})
	B{46}	= ['set Ro	',num2str(Rtab{NewIdx,'Ro'})];
end
if ~isnan(Rtab{NewIdx,'Outer'})
	B{47}	= ['set Outer	',num2str(Rtab{NewIdx,'Outer'})];
end



%% assign integration parameters
B{57}	= ['set FrameTime	',num2str(Rtab{NewIdx,'FrameTime'})];
B{64}	= ['set LangFrac	',num2str(Rtab{NewIdx,'LangFrac'})];
B{70}	= ['set warm_time	',num2str(Rtab{NewIdx,'warm_time'})];
B{71}	= ['set CompactionRatio	',num2str(Rtab{NewIdx,'CompactionRatio'})];

B{73}	= ['set lag_time_init	',int2str(Rtab{NewIdx,'lag_time_init'})];
% B{74}	= ['set lag_time		',int2str(Rtab{NewIdx,'lag_time'})];

B{76}	= ['set comp_step		',int2str(Rtab{NewIdx,'FrameTime'})];	% comp_step = FrameTime
% B{77}	= ['set comp_time		',int2str(Rtab{NewIdx,'comp_time'})];

B{211}	= ['set lj_eps_mw	',num2str(Rtab{NewIdx,'emw'},'%f')];
B{217}	= ['set lj_eps_c1w	',num2str(Rtab{NewIdx,'ec1w'},'%f')];
B{223}	= ['set lj_eps_c2w	',num2str(Rtab{NewIdx,'ec2w'},'%f')];
B{245}	= ['set fene_k	',num2str(Rtab{NewIdx,'fene_k'})];
B{246}	= ['set fene_r	',num2str(Rtab{NewIdx,'fene_r'})];

%% Random Number Generator Seed shift
if ~Rtab{NewIdx,'rand'} || ReRand
	Rtab{NewIdx,'rand'} = uint16(65535*rand);
end
B{140} = ['set rngoffset	',int2str(Rtab{NewIdx,'rand'})];

%% write new file
fid = fopen([header,int2strN(NewIdx,4),'.tcl'],'w');
for i = 1:nL-1
	fprintf(fid,'%s\n', B{i});
end
fprintf(fid,'%s', B{nL});
fclose(fid);
Rtab.tcl{NewIdx}= [header,int2strN(NewIdx,4),'.tcl'];
end
if exist('SmryTabs','var')
	save(Mfile,'LineNum','Rtab','head','Anly','SmryTabs');
else
	save(Mfile,'LineNum','Rtab','head','Anly');
end