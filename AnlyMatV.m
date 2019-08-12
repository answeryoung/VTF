function SimulOut = AnlyMatV(Simul,Anly)
% Simul.VrgMM
% Simul.VrgLL
% Simul.Vrg
% Simul.Vrh
% 
% Simul.VrgMM2d	
% Simul.VrgLL2d
% Simul.Vrg2d
% Simul.Vrh2d
% 
% Simul.VshpC
% Simul.VshpCM
% Simul.VshpCH
% Simul.AshpC
% Simul.AshpCM
% Simul.AshpCH
% Simul.NshpC
% Simul.NshpCM
% Simul.NshpCM
% Simul.RshpC
% Simul.RshpCM
% Simul.RshpCH
%
% DY180618
% DY190215
%%
% TwoDim	= Anly.TwoDim;
Types	= unique(Simul.Atom.type)'+1;
nType	= max(Types);
Pidx	= cell(nType,1);
for tpe = Types
	Pidx{tpe,1} = find(Simul.Atom.type == tpe-1)';
end
Pmidx	= Pidx{1};
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
Nc1		= length(Pcidx1);
Nc2		= length(Pcidx2);
CAM		= Anly.CriticalAlphaMultiple;
pi34	= 4/3*pi;
HoleThreshold	= Anly.HoleThreshold;
%% 3D Volumes
Simul.VrgMM = pi34*Simul.RgMaj.*Simul.RgMi1.*Simul.RgMi2;
Simul.VrgLL = pi34*Simul.RgLong.*Simul.RgLatX.*Simul.RgLatY;
Simul.Vrg	= pi34*Simul.Rg.^3;
Simul.Vrh	= pi34*Simul.Rh.^3;

%% 2D Volumes
Simul.VrgMM2d	= pi34*Simul.RgMaj2d.*Simul.RgMin2d.*Simul.RgMin2d;
Simul.VrgLL2d	= pi34*Simul.RgLong2d.*Simul.RgLat2d.*Simul.RgLat2d;
Simul.Vrg2d		= pi34*Simul.Rg2d.^3;
Simul.Vrh2d		= pi34*Simul.Rh2d.^3;

%% Alpha Shapes
if Anly.RunAlphaShapes
	VshpC	= zeros(Simul.TotalTimeSteps,1);
	VshpCM	= zeros(Simul.TotalTimeSteps,1);
	VshpCH	= zeros(Simul.TotalTimeSteps,1);
	AshpC	= zeros(Simul.TotalTimeSteps,1);
	AshpCM	= zeros(Simul.TotalTimeSteps,1);
	AshpCH	= zeros(Simul.TotalTimeSteps,1);
	N1shpC	= nan(Simul.TotalTimeSteps,1);
	N1shpCM	= nan(Simul.TotalTimeSteps,1);
	N1shpCH	= nan(Simul.TotalTimeSteps,1);
	N2shpC	= nan(Simul.TotalTimeSteps,1);
	N2shpCM	= nan(Simul.TotalTimeSteps,1);
	N2shpCH	= nan(Simul.TotalTimeSteps,1);

	M		= Simul.Coords(Pmidx,:,:);
	C1		= Simul.Coords(Pcidx1,:,:);
	C2		= Simul.Coords(Pcidx2,:,:);
	parfor k = 1:Simul.TotalTimeSteps
	% for k = 1:Simul.TotalTimeSteps
		shpCH	= alphaShape(M(:,:,k),inf);
		VshpCH(k,1)	= volume(shpCH);
		AshpCH(k,1)	= surfaceArea(shpCH);
		if VshpCH(k,1) > 0
			ca		= criticalAlpha(shpCH,'one-region');
			shpC	= alphaShape(M(:,:,k),ca,'HoleThreshold',HoleThreshold);
			shpCM	= alphaShape(M(:,:,k),CAM*ca,'HoleThreshold',HoleThreshold);
			VshpC(k,1)	= volume(shpC);
			VshpCM(k,1)	= volume(shpCM);
			AshpC(k,1)	= surfaceArea(shpC);
			AshpCM(k,1)	= surfaceArea(shpCM);
			N1shpC(k,1)	= sum(inShape(shpC,C1(:,:,k)));
			N1shpCM(k,1)= sum(inShape(shpCM,C1(:,:,k)));
			N1shpCH(k,1)= sum(inShape(shpCH,C1(:,:,k)));
			N2shpC(k,1)	= sum(inShape(shpC,C2(:,:,k)));
			N2shpCM(k,1)= sum(inShape(shpCM,C2(:,:,k)));
			N2shpCH(k,1)= sum(inShape(shpCH,C2(:,:,k)));
		else		
			VshpC(k,1)	= 0;
			VshpCM(k,1)	= 0;
			AshpC(k,1)	= 0;
			AshpCM(k,1)	= 0;
			N1shpC(k,1)	= 0;
			N1shpCM(k,1)	= 0;
			N1shpCH(k,1)	= 0;
			N2shpC(k,1)	= 0;
			N2shpCM(k,1)	= 0;
			N2shpCH(k,1)	= 0;
		end
	end
	Simul.VshpC	= VshpC;
	Simul.VshpCM= VshpCM;
	Simul.VshpCH= VshpCH;

	Simul.AshpC	= AshpC;
	Simul.AshpCM= AshpCM;
	Simul.AshpCH= AshpCH;
	
	if Nc1
		Simul.N1shpC	= N1shpC;
		Simul.N1shpCM	= N1shpCM;
		Simul.N1shpCH	= N1shpCH;
	else
		Simul.N1shpC	= nan(Simul.TotalTimeSteps,1);
		Simul.N1shpCM	= nan(Simul.TotalTimeSteps,1);
		Simul.N1shpCH	= nan(Simul.TotalTimeSteps,1);
	end
	if Nc2
		Simul.N2shpC	= N2shpC;
		Simul.N2shpCM	= N2shpCM;
		Simul.N2shpCH	= N2shpCH;
	else
		Simul.N2shpC	= nan(Simul.TotalTimeSteps,1);
		Simul.N2shpCM	= nan(Simul.TotalTimeSteps,1);
		Simul.N2shpCH	= nan(Simul.TotalTimeSteps,1);
	end
else	
	N1shpC	= Simul.N1shpC;
	N1shpCM	= Simul.N1shpCM;
	N1shpCH	= Simul.N1shpCH;
	N2shpC	= Simul.N2shpC;
	N2shpCM	= Simul.N2shpCM;
	N2shpCH	= Simul.N2shpCH;
	
	VshpC	= Simul.VshpC;
	VshpCM	= Simul.VshpCM;
	VshpCH	= Simul.VshpCH;
end	
%%
Simul.R1shpC	= N1shpC	./ VshpC;
Simul.R1shpCM	= N1shpCM	./ VshpCM;
Simul.R1shpCH	= N1shpCH	./ VshpCH;
Simul.R2shpC	= N2shpC	./ VshpC;
Simul.R2shpCM	= N2shpCM	./ VshpCM;
Simul.R2shpCH	= N2shpCH	./ VshpCH;

Simul.FN1shpC	= N1shpC	./ Simul.NcIn1;
Simul.FN1shpCM	= N1shpCM	./ Simul.NcIn1;
Simul.FN1shpCH	= N1shpCH	./ Simul.NcIn1;
Simul.FN2shpC	= N2shpC	./ Simul.NcIn2;
Simul.FN2shpCM	= N2shpCM	./ Simul.NcIn2;
Simul.FN2shpCH	= N2shpCH	./ Simul.NcIn2;

Simul.FR1shpC	= Simul.R1shpC	./ Simul.RhoIn1;
Simul.FR1shpCM	= Simul.R1shpCM	./ Simul.RhoIn1;
Simul.FR1shpCH	= Simul.R1shpCH	./ Simul.RhoIn1;
Simul.FR2shpC	= Simul.R2shpC	./ Simul.RhoIn2;
Simul.FR2shpCM	= Simul.R2shpCM	./ Simul.RhoIn2;
Simul.FR2shpCH	= Simul.R2shpCH	./ Simul.RhoIn2;

Simul.RO1shpC	= (Simul.NcIn1 - N1shpC)	./ (Simul.V-VshpC);
Simul.RO1shpCM	= (Simul.NcIn1 - N1shpCM)	./ (Simul.V-VshpCM);
Simul.RO1shpCH	= (Simul.NcIn1 - N1shpCH)	./ (Simul.V-VshpCH);
Simul.RO2shpC	= (Simul.NcIn2 - N2shpC)	./ (Simul.V-VshpC);
Simul.RO2shpCM	= (Simul.NcIn2 - N2shpCM)	./ (Simul.V-VshpCM);
Simul.RO2shpCH	= (Simul.NcIn2 - N2shpCH)	./ (Simul.V-VshpCH);

Simul.Rr1shpC	= Simul.RO1shpC	./ Simul.R1shpC;
Simul.Rr1shpCM	= Simul.RO1shpCM./ Simul.R1shpCM;
Simul.Rr1shpCH	= Simul.RO1shpCH./ Simul.R1shpCH;
Simul.Rr2shpC	= Simul.RO2shpC	./ Simul.R2shpC;
Simul.Rr2shpCM	= Simul.RO2shpCM./ Simul.R2shpCM;
Simul.Rr2shpCH	= Simul.RO2shpCH./ Simul.R2shpCH;

if Nc1
	Simul.Rd1shpC	= Simul.RO1shpC	- Simul.R1shpC;
	Simul.Rd1shpCM	= Simul.RO1shpCM- Simul.R1shpCM;
	Simul.Rd1shpCH	= Simul.RO1shpCH- Simul.R1shpCH;
else
	Simul.Rd1shpC	= zeros(Simul.TotalTimeSteps,1);
	Simul.Rd1shpCM	= zeros(Simul.TotalTimeSteps,1);
	Simul.Rd1shpCH	= zeros(Simul.TotalTimeSteps,1);
end	
if Nc2
	Simul.Rd2shpC	= Simul.RO2shpC	- Simul.R2shpC;
	Simul.Rd2shpCM	= Simul.RO2shpCM- Simul.R2shpCM;
	Simul.Rd2shpCH	= Simul.RO2shpCH- Simul.R2shpCH;
else
	Simul.Rd2shpC	= zeros(Simul.TotalTimeSteps,1);
	Simul.Rd2shpCM	= zeros(Simul.TotalTimeSteps,1);
	Simul.Rd2shpCH	= zeros(Simul.TotalTimeSteps,1);
end	
SimulOut	= Simul;