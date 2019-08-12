function SimulOut = AnlyMat(Simul,Anly,L,R)
% Simul.Rg	= Rg;
% Simul.Rh	= Rh;
% Simul.Re	= Re;
% Simul.Sg	= Sg;
% Simul.E	= E;
% Simul.A	= A; = RgZ/sqrt(RgX*RgY);
% Simul.Ap	= Ap;= RgMaj/sqrt(RgMi1*RgMi2);
% Simul.R	= R;
% Simul.CoM	= C;
% 
% Simul.RgLong	= RgLong;
% Simul.RgLatX	= RgLatX;
% Simul.RgLatY	= RgLatY;
% Simul.RgMaj	= RgMaj;
% Simul.RgMi1	= RgMi1;
% Simul.RgMi2	= RgMi2;
% 
% Simul.Rg2d	= Rg2d;
% Simul.Rh2d	= Rh2d;
% Simul.Sg2d	= Sg2d;
% Simul.E2d		= E2d;
% Simul.A2d		= A2d;
% Simul.Ap2d	= Ap2d;
% 
% 2D measures Z-X
% Simul.RgLong2d= RgLong2d;
% Simul.RgLat2d	= RgLat2d;
% Simul.RgMaj2d	= RgMaj2d;
% Simul.RgMin2d	= RgMin2d;
%
% Simul.HistZ	= HistZ;
% Simul.HistR	= HistR;
% Simul.HistX	= HistX;
% Simul.HistZEs	= ZEs;
% Simul.HistREs	= REs;
% Simul.HistXEs	= XEs;
% 
% DY180620
% DY190215
%%
cx	= Anly.cx;
cy	= Anly.cy;
cz	= Anly.cz;
TwoDim	= Anly.TwoDim;

Types	= unique(Simul.Atom.type)'+1;
nType	= max(Types);
Pidx	= cell(nType,1);
for tpe = Types
	Pidx{tpe,1} = find(Simul.Atom.type == tpe-1)';
end
Pmidx = Pidx{1};
if max(Types) > 1
	Pcidx1	= Pidx{2};
else
	Pcidx1	= [];
end
if max(Types) > 2
	Pcidx2	= Pidx{3};
else
	Pcidx2	= [];
end

Simul.V	= 2*pi*L*R^2;
%% 3D measures & Rh2d
Rg	= zeros(Simul.TotalTimeSteps,1);
Rh	= zeros(Simul.TotalTimeSteps,1);
Re	= zeros(Simul.TotalTimeSteps,1);
Rm	= zeros(Simul.TotalTimeSteps,length(Pmidx));
% Rc	= zeros(Simul.TotalTimeSteps,length(Pcidx));
A	= zeros(Simul.TotalTimeSteps,1);
Ap	= zeros(Simul.TotalTimeSteps,1);
Sg  = zeros(3,3,Simul.TotalTimeSteps);
E	= zeros(Simul.TotalTimeSteps,3);
CoM	= zeros(Simul.TotalTimeSteps,3);

RgLong	= zeros(Simul.TotalTimeSteps,1);
RgLatX	= zeros(Simul.TotalTimeSteps,1);
RgLatY	= zeros(Simul.TotalTimeSteps,1);
RgMaj	= zeros(Simul.TotalTimeSteps,1);
RgMi1	= zeros(Simul.TotalTimeSteps,1);
RgMi2	= zeros(Simul.TotalTimeSteps,1);

Rh2d	= zeros(Simul.TotalTimeSteps,1);

M(:,3,:)= Simul.Coords(Pmidx,3,:) - cz;
M(:,2,:)= Simul.Coords(Pmidx,2,:) - cy;
M(:,1,:)= Simul.Coords(Pmidx,1,:) - cx;
Mx		= permute(M(:,1,:),[1,3,2]);
My		= permute(M(:,2,:),[1,3,2]);
Mz		= permute(M(:,3,:),[1,3,2]);
M2		= M(:,TwoDim,:);

if ~isempty(Pcidx1)
	Cz1	= permute(Simul.Coords(Pcidx1,3,:) - cz,[1,3,2]);
	Cy1	= permute(Simul.Coords(Pcidx1,2,:) - cy,[1,3,2]);
	Cx1	= permute(Simul.Coords(Pcidx1,1,:) - cx,[1,3,2]);
end
if ~isempty(Pcidx2)
	Cz2	= permute(Simul.Coords(Pcidx2,3,:) - cz,[1,3,2]);
	Cy2	= permute(Simul.Coords(Pcidx2,2,:) - cy,[1,3,2]);
	Cx2	= permute(Simul.Coords(Pcidx2,1,:) - cx,[1,3,2]);
end

parfor k = 1:Simul.TotalTimeSteps
% for k = 1:Simul.TotalTimeSteps
% 	F = Simul.Coords(Pmidx,:,k);
% 	F(:,1) = F(:,1) - cx;
% 	F(:,2) = F(:,2) - cy;
% 	F(:,3) = F(:,3) - cz;
	[Rg(k,1),Sgt,Et,CoM(k,:)]	= rg(M(:,:,k));
	Rh(k,1)	= rh(M(:,:,k));
	Re(k,1) = range(Mz(:,k))/2;

	RgLong(k,1) = sqrt(Sgt(3,3));
	RgLatX(k,1) = sqrt(Sgt(1,1));
	RgLatY(k,1) = sqrt(Sgt(2,2));	
	RgMaj(k,1)	= sqrt(Et(1,3));
	RgMi1(k,1)	= sqrt(Et(1,1));
	RgMi2(k,1)	= sqrt(Et(1,2));	
% 	A(k,1)	= RgLong(k,1)/RgLatX(k,1);
% 	Ap(k,1)	= RgMaj(k,1)/RgMi2(k,1);
	A(k,1)	= RgLong(k,1)/sqrt(RgLatX(k,1)*RgLatY(k,1));
 	Ap(k,1)	= RgMaj(k,1)/sqrt(RgMi1(k,1)*RgMi2(k,1));
	Rm(k,:) = sqrt(Mx(:,k).^2 + My(:,k).^2)';

	Sg(:,:,k)	= Sgt;
	E(k,:)		= Et;
	
	Rh2d(k,1)	= rh(M2(:,:,k));
end
Simul.Rg	= Rg;
Simul.Rh	= Rh;
Simul.Re	= Re;
Simul.Sg	= Sg;
Simul.E		= E;
Simul.A		= A;
Simul.Ap	= Ap;
Simul.Rm	= Rm;
Simul.CoM	= CoM;

Simul.RgLong= RgLong;
Simul.RgLatX= RgLatX;
Simul.RgLatY= RgLatY;
Simul.RgMaj	= RgMaj;
Simul.RgMi1	= RgMi1;
Simul.RgMi2	= RgMi2;

%% 2D measures Z-X, if TwoDim == [3,1]
[Rg2d,Sg2d,E2d] = SgDimReduc(Sg,TwoDim);

RgLong2d(:,1)	= sqrt(Sg2d(1,1,:));
RgLat2d(:,1)	= sqrt(Sg2d(2,2,:));
RgMaj2d(:,1)	= sqrt(E2d(:,2));
RgMin2d(:,1)	= sqrt(E2d(:,1));
A2d(:,1)		= RgLong2d./RgLat2d;
Ap2d(:,1)		= RgMaj2d./RgMin2d;

Simul.Rg2d	= Rg2d;
Simul.Rh2d	= Rh2d;
Simul.Sg2d	= Sg2d;
Simul.E2d	= E2d;
Simul.A2d	= A2d;
Simul.Ap2d	= Ap2d;

Simul.RgLong2d	= RgLong2d;
Simul.RgLat2d	= RgLat2d;
Simul.RgMaj2d	= RgMaj2d;
Simul.RgMin2d	= RgMin2d;
%% Histograms
Zmax= Anly.Zmax;
dZ	= Anly.ZHistWidth;
ZEs = (-Zmax+dZ/2):dZ:(Zmax-dZ/2);
nBZ	= length(ZEs)-1;

Rmax= Anly.Rmax;
dR	= Anly.RHistWidth;
REs = 0:dR:(Rmax-dR);
nBR	= length(REs)-1;

Xmax= Anly.Xmax;
dX	= Anly.XHistWidth;
XEs = (-Xmax+dX/2):dX:(Xmax-dX/2);
nBX	= length(XEs)-1;

HistZm	= nan(Simul.TotalTimeSteps,nBZ);
HistRm	= nan(Simul.TotalTimeSteps,nBR);
HistXm	= nan(Simul.TotalTimeSteps,nBX);
HistZc1	= nan(Simul.TotalTimeSteps,nBZ);
HistRc1	= nan(Simul.TotalTimeSteps,nBR);
HistXc1	= nan(Simul.TotalTimeSteps,nBX);
HistZc2	= nan(Simul.TotalTimeSteps,nBZ);
HistRc2	= nan(Simul.TotalTimeSteps,nBR);
HistXc2	= nan(Simul.TotalTimeSteps,nBX);
NcIn1	= zeros(Simul.TotalTimeSteps,1);
NcIn2	= zeros(Simul.TotalTimeSteps,1);

parfor k = 1:Simul.TotalTimeSteps
% for k = 1:Simul.TotalTimeSteps
	HistZm(k,:) = histcounts(Mz(:,k),ZEs)';
	HistXm(k,:) = histcounts(Mx(:,k),XEs)';
	HistRm(k,:) = histcounts(Rm(k,:)',REs)';
end

if ~isempty(Pcidx1) && ~isempty(Pcidx2)
	parfor k = 1:Simul.TotalTimeSteps
		HistZc1(k,:) = histcounts(Cz1(:,k),ZEs)';
		HistXc1(k,:) = histcounts(Cx1(:,k),XEs)';
		HistZc2(k,:) = histcounts(Cz2(:,k),ZEs)';
		HistXc2(k,:) = histcounts(Cx2(:,k),XEs)';
		
		Rc1 = sqrt(Cx1(:,k).^2 + Cy1(:,k).^2);
		Rc2 = sqrt(Cx2(:,k).^2 + Cy2(:,k).^2);
		
		HistRc1(k,:) = histcounts(Rc1,REs)';
		HistRc2(k,:) = histcounts(Rc2,REs)';
		
		NcIn1(k,1) = sum(abs(Cz1(:,k))<=L & abs(Rc1)<=R);
		NcIn2(k,1) = sum(abs(Cz2(:,k))<=L & abs(Rc2)<=R);
	end
elseif ~isempty(Pcidx1)
	parfor k = 1:Simul.TotalTimeSteps
		HistZc1(k,:) = histcounts(Cz1(:,k),ZEs)';
		HistXc1(k,:) = histcounts(Cx1(:,k),XEs)';
		
		Rc1 = sqrt(Cx1(:,k).^2 + Cy1(:,k).^2);
		HistRc1(k,:) = histcounts(Rc1,REs)';
		
		NcIn1(k,1) = sum(abs(Cz1(:,k))<=L & abs(Rc1)<=R);
	end
elseif ~isempty(Pcidx2)
	parfor k = 1:Simul.TotalTimeSteps
		HistZc2(k,:) = histcounts(Cz2(:,k),ZEs)';
		HistXc2(k,:) = histcounts(Cx2(:,k),XEs)';
		
		Rc2 = sqrt(Cx2(:,k).^2 + Cy2(:,k).^2);
		HistRc2(k,:) = histcounts(Rc2,REs)';
		
		NcIn2(k,1) = sum(abs(Cz2(:,k))<=L & abs(Rc2)<=R);
	end		
end
if isfield(Simul,'HistZ')
	Simul	= rmfield(Simul,{'HistZ','HistR','HistX'});
end
Simul.HistZ(:,:,3)	= HistZc2;
Simul.HistR(:,:,3)	= HistRc2;
Simul.HistX(:,:,3)	= HistXc2;
Simul.HistZ(:,:,2)	= HistZc1;
Simul.HistR(:,:,2)	= HistRc1;
Simul.HistX(:,:,2)	= HistXc1;
Simul.HistZ(:,:,1)	= HistZm;
Simul.HistR(:,:,1)	= HistRm;
Simul.HistX(:,:,1)	= HistXm;
Simul.HistZEs	= ZEs;
Simul.HistREs	= REs;
Simul.HistXEs	= XEs;

Simul.NcIn1			= NcIn1;
Simul.NcIn2			= NcIn2;
Simul.RhoIn1(:,1)	= NcIn1/Simul.V;
Simul.RhoIn2(:,1)	= NcIn2/Simul.V;

SimulOut		= Simul;