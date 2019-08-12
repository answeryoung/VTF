function Simul = AnlySimImg(Simul,Rtab,Anly,Dim)
%
% DY180628
%%
if ~isfield(Simul,'Image')
	Simul = SimImg(Simul,Anly,Dim);
end
Nim = size(Simul.Image,4);
nucBgStr	= 'Mean';
% nucBgStr	= 'Fit';
BgMean	= 150;
GenBg	= @(Img) imnoise(Img,'gaussian',1e-5,2e-7);
% GenBg	= @(Img) Img;
global WidthCorr PsfWidth
WidthCorr	= @(sig,psf) 2*sig+0.08;
% WidthCorr	= @(sig,psf) 2*sqrt(sig^2-psf^2);
bgMax = 1.2*BgMean;
%%
global br_im br_long br_lat r_fit hlFitWidth bgbrm
br_im	= 4;
br_long = 2;
br_lat	= 2;
r_fit	= 17;
hlFitWidth	= 17;
bgbrm	= 2;
global calibrf
calibrf			= Anly.calibrf;
[~,CytoPsfWidth]= PSFwidth(0,Anly.RedEm);
[~,NucPsfWidth]	= PSFwidth(0,Anly.GreenEm);

L	= Rtab.L*Anly.Unit/calibrf;
R	= Rtab.R*Anly.Unit/calibrf;
CC	= ([Anly.SimImgY,Anly.SimImgX]+1)/2;
ps	= [CC;CC];
ps(1,2) = ps(1,2)-L;
ps(2,2) = ps(2,2)+L;
[Lpx,k,dV] = PorPS(ps);
dx	= (ceil(Lpx) - Lpx) * 0.5;
Coords	= zeros(ceil(Lpx)+1,2);
for x = (0:ceil(Lpx))-dx
	Coords(x+1+dx,:) = ps(1,:) + x*dV;
end
cytoLength	= Lpx * calibrf;
%%
Rtab = table;
for t = 1:Nim
%%
	imNuc	= uint16(65536*Simul.Image(:,:,1,t));
	imNuc	= imNuc + GenBg(0*imNuc + BgMean);
	if size(Simul.Image,3) == 1
		imCyto1	= 0*imNuc;
		imCyto2	= 0*imNuc;
	elseif size(Simul.Image,3) == 2
		imCyto1	= uint16(65536*Simul.Image(:,:,2,t));
		imCyto2	= 0*imNuc;
	else
		imCyto1	= uint16(65536*Simul.Image(:,:,2,t));
		imCyto2	= uint16(65536*Simul.Image(:,:,3,t));
	end
	imCyto1	= imCyto1 + GenBg(0*imCyto1 + BgMean);
	imCyto2	= imCyto2 + GenBg(0*imCyto2 + BgMean);
	imCyto	= imCyto1 + imCyto2;
	
	[cytoWidthRawSig1,cytoRawBg1]	= FitCytoWidth2(imCyto1,Coords,Lpx,k);
	[cytoWidthRawSig2,cytoRawBg2]	= FitCytoWidth2(imCyto2,Coords,Lpx,k);
	[cytoWidthRawSig,cytoRawBg]		= FitCytoWidth2(imCyto,Coords,Lpx,k);
	
	cytoWidth1	= WidthCorr(mean(cytoWidthRawSig1)*calibrf,CytoPsfWidth);
	cytoVolume1	= pi*cytoWidth1*cytoWidth1*(cytoLength/4-cytoWidth1/12);
	cytoImg1	= imProfileW(imCyto1,Coords,br_im);
	MeanCytoBg1	= mean(cytoRawBg1);
	cytoWidth2	= WidthCorr(mean(cytoWidthRawSig2)*calibrf,CytoPsfWidth);
	cytoVolume2	= pi*cytoWidth2*cytoWidth2*(cytoLength/4-cytoWidth2/12);
	cytoImg2	= imProfileW(imCyto2,Coords,br_im);
	MeanCytoBg2	= mean(cytoRawBg2);
	cytoWidth	= WidthCorr(mean(cytoWidthRawSig)*calibrf,CytoPsfWidth);
	cytoVolume	= pi*cytoWidth*cytoWidth*(cytoLength/4-cytoWidth/12);
	cytoImg		= imProfileW(imCyto,Coords,br_im);
	MeanCytoBg	= mean(cytoRawBg);
	
	PsfWidth	= CytoPsfWidth;
	[cytoVolInt1,cytoR01,cytoI01]	= ...
		VolumeIntensity(cytoImg1,imCyto1,Coords,k,MeanCytoBg1);
	[cytoVolInt2,cytoR02,cytoI02]	= ...
		VolumeIntensity(cytoImg2,imCyto2,Coords,k,MeanCytoBg2);
	[cytoVolInt,cytoR0,cytoI0]	= ...
		VolumeIntensity(cytoImg,imCyto,Coords,k,MeanCytoBg);

		Yim = br_im + 1 + (-br_long:br_long);
		LcytoImg1	= size(cytoImg1,1);
		LcytoImg2	= size(cytoImg2,1);
		LcytoImg	= size(cytoImg,1);
		if LcytoImg1 > 12
			MeanCytoIntensity1	= mean(mean(cytoImg1(6:end-5,Yim)));
		else
			Xim		= floor(LcytoImg1/2):ceil(LcytoImg1/2)+1;
			MeanCytoIntensity1	= mean(mean(cytoImg1(Xim,Yim)));
		end
		if LcytoImg1 > 12
			MeanCytoIntensity2	= mean(mean(cytoImg2(6:end-5,Yim)));
		else
			Xim		= floor(LcytoImg2/2):ceil(LcytoImg2/2)+1;
			MeanCytoIntensity2	= mean(mean(cytoImg2(Xim,Yim)));
		end
		if LcytoImg1 > 12
			MeanCytoIntensity	= mean(mean(cytoImg(6:end-5,Yim)));
		else
			Xim		= floor(LcytoImg/2):ceil(LcytoImg/2)+1;
			MeanCytoIntensity	= mean(mean(cytoImg(Xim,Yim)));
		end
		MI1					= MIplateau(cytoImg1,MeanCytoBg1);
		TotalCytoIntensity1	= sum(cytoImg1(:)) - MeanCytoBg1;
		MeanCytoPlateauIntensity1	= MI1 - MeanCytoBg1;
		MI2					= MIplateau(cytoImg2,MeanCytoBg2);
		TotalCytoIntensity2	= sum(cytoImg2(:)) - MeanCytoBg2;
		MeanCytoPlateauIntensity2	= MI2 - MeanCytoBg2;
		MI					= MIplateau(cytoImg,MeanCytoBg);
		TotalCytoIntensity	= sum(cytoImg(:)) - MeanCytoBg;
		MeanCytoPlateauIntensity= MI - MeanCytoBg;

	cytoAspect1		= cytoLength./cytoWidth1;
	cytoAspect2		= cytoLength./cytoWidth2;
	cytoAspect		= cytoLength./cytoWidth;

% 	cytoLengthInv1	= cytoLength1.^-1;
	cytoVolumeInv1	= cytoVolume1.^-1;
	cytoVolIntInv1	= cytoVolInt1.^-1;
	MeanTotalCytoIntensity1...
		= MeanCytoIntensity1./TotalCytoIntensity1;
% 	cytoLengthInv2	= cytoLength2.^-1;
	cytoVolumeInv2	= cytoVolume2.^-1;
	cytoVolIntInv2	= cytoVolInt2.^-1;
	MeanTotalCytoIntensity2...
		= MeanCytoIntensity2./TotalCytoIntensity2;
	cytoLengthInv	= cytoLength.^-1;
	cytoVolumeInv	= cytoVolume.^-1;
	cytoVolIntInv	= cytoVolInt.^-1;
	MeanTotalCytoIntensity...
		= MeanCytoIntensity./TotalCytoIntensity;

%%
	[nInfOff,~,pb]		= FitNucLength(imNuc,Coords,k);
	if 2*nInfOff*calibrf < 2 * cytoLength
		[~,~,dV] = PorPS(pb);
		CW	= [-1;1]*nInfOff*0.4*dV + [mean(pb);mean(pb)];
% 			putvar(im,CW,k_long,nInfOff);
		[nucWidthRawSig,nucRawBg]	= FitNucWidth2(imNuc,CW,k+pi/2);
		nucLength	= 2*nInfOff*calibrf;
		nucWidth	= WidthCorr(mean(nucWidthRawSig)*calibrf,NucPsfWidth);
		nucVolume	= pi*nucWidth*nucWidth*(nucLength/4-nucWidth/12);
		nucImg		= imProfileW(imNuc,Coords,br_im);
		nucBg.Mean	= mean(nucRawBg);

%%
%			[nucRg,nucRgLong,nucRgLat,nucRgMa,nucRgMi,nucGyrationTensor,MeanNucIntensity]...
%				= rgNuc2(nucImg,MeanNucBg);			
		[nucRg,nucRgLong,nucRgLat,nucRgMa,nucRgMi,nucGyrationTensor,MeanNucIntensity,nucBg.Fit]...
			= rgNuc(nucImg,bgbrm,bgMax);
		[nucVolInt,nucR0,nucI0]	= ...
			VolumeIntensity(nucImg,imNuc,Coords,k,nucBg.(nucBgStr));
		MI	= MIplateau(nucImg,nucBg.(nucBgStr));
		TotalNucIntensity		= sum(nucImg(:)) - nucBg.(nucBgStr);
		MeanNucPlateauIntensity = MI - nucBg.(nucBgStr);

	else
		nucLength	= 0;
		nucWidth	= 0;
		nucVolume	= 0;
		nucWidthRawSig	= 0;			
		nucImg			= imProfileW(imNuc,Coords,br_im);
		nucRawBg		= mean(nucImg(:));
		nucBg.Mean		= mean(nucImg(:));

		nucRg		= nan;
		nucRgLong	= nan;
		nucRgLat	= nan;
		nucRgMa		= nan;
		nucRgMi		= nan;
		nucGyrationTensor	= nan;
		MeanNucIntensity	= 0;
		nucBg.Fit	= mean(nucImg(:));
		nucVolInt	= 0;
		nucR0		= 0;
		nucI0		= 0;
		TotalNucIntensity	= 0;
		MeanNucPlateauIntensity		= 0;
	end
		
	nucAspect		= nucLength./nucWidth;
	RgLAspect		= nucRgLong./nucRgLat;
	RgMAspect		= nucRgMa./nucRgMi;
	
% 	LengthFraction1	= nucLength./cytoLength1;
	WidthFraction1	= nucWidth./cytoWidth1;
	VolumeFraction1	= nucVolume./cytoVolume1;
	VolIntFraction1	= nucVolume./cytoVolInt1;
	I0Ratio1		= nucI0./cytoI01;
	MeanIntensityRatio1	= MeanNucIntensity./MeanCytoIntensity1;
	MeanTotalNucIntensity...
		= MeanNucIntensity./TotalNucIntensity;
	MeanTotalIntensityRatio1...
		= MeanTotalNucIntensity./MeanTotalCytoIntensity1;		

% 	LengthFraction2	= nucLength./cytoLength2;
	WidthFraction2	= nucWidth./cytoWidth2;
	VolumeFraction2	= nucVolume./cytoVolume2;
	VolIntFraction2	= nucVolume./cytoVolInt2;
	I0Ratio2		= nucI0./cytoI02;
	MeanIntensityRatio2	= MeanNucIntensity./MeanCytoIntensity2;
% 	MeanTotalNucIntensity...
% 		= MeanNucIntensity./TotalNucIntensity;
	MeanTotalIntensityRatio2...
		= MeanTotalNucIntensity./MeanTotalCytoIntensity2;
	
	LengthFraction	= nucLength./cytoLength;
	WidthFraction	= nucWidth./cytoWidth;
	VolumeFraction	= nucVolume./cytoVolume;
	VolIntFraction	= nucVolume./cytoVolInt;
	I0Ratio			= nucI0./cytoI0;
	MeanIntensityRatio	= MeanNucIntensity./MeanCytoIntensity;
% 	MeanTotalNucIntensity...
% 		= MeanNucIntensity./TotalNucIntensity;
	MeanTotalIntensityRatio...
		= MeanTotalNucIntensity./MeanTotalCytoIntensity;
%%
	Rtab.cytoLength(t,1)= cytoLength;
	Rtab.cytoWidth(t,1)	= cytoWidth;
	Rtab.cytoR0(t,1)	= cytoR0;
	Rtab.cytoVolume(t,1)= cytoVolume;
	Rtab.cytoVolInt(t,1)= cytoVolInt;
	Rtab.cytoI0(t,1)	= cytoI0;
	Rtab.cytoWidth1(t,1)	= cytoWidth1;
	Rtab.cytoR01(t,1)		= cytoR01;
	Rtab.cytoVolume1(t,1)	= cytoVolume1;
	Rtab.cytoVolInt1(t,1)	= cytoVolInt1;
	Rtab.cytoI01(t,1)		= cytoI01;
	Rtab.cytoWidth2(t,1)	= cytoWidth2;
	Rtab.cytoR02(t,1)		= cytoR02;
	Rtab.cytoVolume2(t,1)	= cytoVolume2;
	Rtab.cytoVolInt2(t,1)	= cytoVolInt2;
	Rtab.cytoI02(t,1)		= cytoI02;		
	
	Rtab.nucLength(t,1)	= nucLength;
	Rtab.nucWidth(t,1)	= nucWidth;
	Rtab.nucR0(t,1)		= nucR0;
	Rtab.nucVolume(t,1) = nucVolume;
	Rtab.nucVolInt(t,1) = nucVolInt;
	Rtab.nucI0(t,1)		= nucI0;
	Rtab.nucRg(t,1)		= nucRg*calibrf;
	Rtab.nucRgLong(t,1)	= nucRgLong*calibrf;
	Rtab.nucRgLat(t,1)	= nucRgLat*calibrf;
	Rtab.nucRgMa(t,1)	= nucRgMa*calibrf;
	Rtab.nucRgMi(t,1)	= nucRgMi*calibrf;
	
	Rtab.MeanCytoIntensity(t,1)			= MeanCytoIntensity - MeanCytoBg;
	Rtab.TotalCytoIntensity(t,1)		= TotalCytoIntensity;
	Rtab.MeanCytoPlateauIntensity(t,1)	= MeanCytoPlateauIntensity;
	Rtab.MeanCytoBg(t,1)				= MeanCytoBg;
	Rtab.MeanCytoIntensity1(t,1)		= MeanCytoIntensity1 - MeanCytoBg1;
	Rtab.TotalCytoIntensity1(t,1)		= TotalCytoIntensity1;
	Rtab.MeanCytoPlateauIntensity1(t,1)	= MeanCytoPlateauIntensity1;
	Rtab.MeanCytoBg1(t,1)				= MeanCytoBg1;
	Rtab.MeanCytoIntensity2(t,1)		= MeanCytoIntensity2 - MeanCytoBg2;
	Rtab.TotalCytoIntensity2(t,1)		= TotalCytoIntensity2;
	Rtab.MeanCytoPlateauIntensity2(t,1)	= MeanCytoPlateauIntensity2;
	Rtab.MeanCytoBg2(t,1)				= MeanCytoBg2;
	
	Rtab.MeanNucIntensity(t,1)			= MeanNucIntensity - nucBg.(nucBgStr);
	Rtab.TotalNucIntensity(t,1)			= TotalNucIntensity;
	Rtab.MeanNucPlateauIntensity(t,1)	= MeanNucPlateauIntensity;
	Rtab.MeanNucBg(t,1)					= nucBg.Mean;
	Rtab.nucBgFit(t,1)					= nucBg.Fit;
	
	Rtab.cytoWidthRawSig{t,1}	= cytoWidthRawSig;
	Rtab.cytoRawBg{t,1}			= cytoRawBg;
	Rtab.cytoImg{t,1}			= cytoImg;
	Rtab.cytoWidthRawSig1{t,1}	= cytoWidthRawSig1;
	Rtab.cytoRawBg1{t,1}		= cytoRawBg1;
	Rtab.cytoImg1{t,1}			= cytoImg1;
	Rtab.cytoWidthRawSig1{t,1}	= cytoWidthRawSig2;
	Rtab.cytoRawBg2{t,1}		= cytoRawBg2;
	Rtab.cytoImg2{t,1}			= cytoImg2;
	Rtab.cytoWidthRawSig{t,1}	= cytoWidthRawSig;
	Rtab.cytoRawBg{t,1}			= cytoRawBg;
	Rtab.cytoImg{t,1}			= cytoImg;
	
	Rtab.nucWidthRawSig{t,1}	= nucWidthRawSig;
	Rtab.nucRawBg{t,1}			= nucRawBg;		
	Rtab.nucImg{t,1}			= nucImg;
	Rtab.nucGyrationTensor{t,1}	= nucGyrationTensor;

	Rtab.cytoAspect(t,1)	= cytoAspect;
	Rtab.cytoAspect1(t,1)	= cytoAspect1;
	Rtab.cytoAspect2(t,1)	= cytoAspect2;
	Rtab.nucAspect(t,1)		= nucAspect;
	Rtab.RgLAspect(t,1)		= RgLAspect;
	Rtab.RgMAspect(t,1)		= RgMAspect;
	
	Rtab.LengthFraction(t,1)	= LengthFraction;
	Rtab.WidthFraction(t,1)		= WidthFraction;
	Rtab.VolumeFraction(t,1)	= VolumeFraction;
	Rtab.VolIntFraction(t,1)	= VolIntFraction;
	Rtab.I0Ratio(t,1)			= I0Ratio;
	Rtab.MeanIntensityRatio(t,1)		= MeanIntensityRatio;
	Rtab.MeanTotalCytoIntensity(t,1)	= MeanTotalCytoIntensity;
	Rtab.WidthFraction1(t,1)	= WidthFraction1;
	Rtab.VolumeFraction1(t,1)	= VolumeFraction1;
	Rtab.VolIntFraction1(t,1)	= VolIntFraction1;
	Rtab.I0Ratio1(t,1)			= I0Ratio1;
	Rtab.MeanIntensityRatio1(t,1)		= MeanIntensityRatio1;
	Rtab.MeanTotalCytoIntensity1(t,1)	= MeanTotalCytoIntensity1;
	Rtab.WidthFraction2(t,1)	= WidthFraction2;
	Rtab.VolumeFraction2(t,1)	= VolumeFraction2;
	Rtab.VolIntFraction2(t,1)	= VolIntFraction2;
	Rtab.I0Ratio2(t,1)			= I0Ratio2;
	Rtab.MeanIntensityRatio2(t,1)		= MeanIntensityRatio2;
	Rtab.MeanTotalCytoIntensity2(t,1)	= MeanTotalCytoIntensity2;
	
	Rtab.cytoLengthInv(t,1) = cytoLengthInv;
	Rtab.cytoVolumeInv(t,1) = cytoVolumeInv;
	Rtab.cytoVolIntInv(t,1) = cytoVolIntInv;
	Rtab.MeanTotalNucIntensity(t,1)		= MeanTotalNucIntensity;
	Rtab.MeanTotalIntensityRatio(t,1)	= MeanTotalIntensityRatio;
% 	Rtab.cytoLengthInv1(t,1) = cytoLengthInv1;
	Rtab.cytoVolumeInv1(t,1) = cytoVolumeInv1;
	Rtab.cytoVolIntInv1(t,1) = cytoVolIntInv1;
% 	Rtab.MeanTotalNucIntensity(t,1)		= MeanTotalNucIntensity;
	Rtab.MeanTotalIntensityRatio1(t,1)	= MeanTotalIntensityRatio1;
% 	Rtab.cytoLengthInv2(t,1) = cytoLengthInv2;
	Rtab.cytoVolumeInv2(t,1) = cytoVolumeInv2;
	Rtab.cytoVolIntInv2(t,1) = cytoVolIntInv2;
% 	Rtab.MeanTotalNucIntensity(t,1)		= MeanTotalNucIntensity;
	Rtab.MeanTotalIntensityRatio2(t,1)	= MeanTotalIntensityRatio2;

end
Simul.ImTab = Rtab;
% Rwrite	= Rtab(tset,{'Time_s','cytoLength','cytoWidth','cytoR0','cytoVolume','cytoVolInt',...
% 	'nucLength','nucWidth','nucR0','nucVolume','nucVolInt',...
% 	'nucRg','nucRgLong','nucRgLat',...
% 	'nucRgMa','nucRgMi',...
% 	'MeanCytoIntensity','TotalCytoIntensity','MeanCytoPlateauIntensity',...
% 	'cytoI0','MeanCytoBg',...
% 	'MeanNucIntensity','TotalNucIntensity','MeanNucPlateauIntensity',...
% 	'nucI0','MeanNucBg','nucBgFit',...
% 	'cytoAspect','nucAspect',...
% 	'RgLAspect','RgMAspect',...
% 	'LengthFraction','WidthFraction',...
% 	'VolumeFraction','VolIntFraction',...
% 	'I0Ratio','MeanIntensityRatio','MeanTotalIntensityRatio',...
% 	'MeanTotalCytoIntensity',...
% 	'cytoLengthInv','cytoVolumeInv','cytoVolIntInv',...
% 	'MeanTotalNucIntensity'});
% Rwrite{:,'Time_s'} = Rwrite{:,'Time_s'}/60;
% Rwrite.Properties.VariableNames{1} = 'Time';
% 
% Frame	= tset';
% Frame	= table(Frame);
% Rwrite	= [Frame,Rwrite];
% LongName= {'Frame','Time','Length','Width','Max Radius','Volume','Volume',...
% 	'Length','Width','Max Radius','Volume','Volume',...
% 	'Radius of gyration', 'Longitudinal radius of gyration','Lateral radius of gyration',...
% 	'Longitudinal radius of gyration', 'Lateral radius of gyration',...
% 	'Mean Cytoplasmic Intensity','Total Cytoplasmic Intensity','Plateau Cytoplasmic Intensity',...
% 	'Max Cytoplasmic Intensity','Mean Cytoplasmic Background',...
% 	'Mean Nuclidic Intensity','Total Nuclidic Intensity','Plateau Nuclidic Intensity',...
% 	'Max Nuclidic Intensity','Mean Nuclidic Background','Fitted Nuclidic Background',...
% 	'Cytoplasmic aspect ratio', 'Nucleoid aspect ratio',...
% 	'Rg aspect ratio L/L','Rg Aspect ratio M/M',...
% 	'Nucleoid length fraction','Nucleoid width fraction',...
% 	'Nucleoid volume fraction','Nucleoid volint fraction',...
% 	'Max intensity ratio','Mean intensity ratio','Mean/Total intensity ratio',...
% 	'Mean/Total intensity',...
% 	'Inverse length','Inverse volume','Inverse volume',...
% 	'Mean/Total intensity'};
% Units	= {'','min','\g(m)m','\g(m)m','\g(m)m','\g(m)m\+(3)','\g(m)m\+(3)',...
% 	'\g(m)m','\g(m)m','\g(m)m','\g(m)m\+(3)','\g(m)m\+(3)',...
% 	'\g(m)m','\g(m)m','\g(m)m',...
% 	'\g(m)m','\g(m)m',...
% 	'arb. unit','arb. unit','arb. unit',...
% 	'arb. unit','arb. unit',...
% 	'arb. unit','arb. unit','arb. unit',...
% 	'arb. unit','arb. unit','arb. unit',...
% 	'','',...
% 	'','',...	
% 	'','',...
% 	'','',...
% 	'','','',...
% 	'',...
% 	'\g(m)m\+(-1)','\g(m)m\+(-3)','\g(m)m\+(-3)',...
% 	''};
% 	cyto= [BactName,' cyto'];
% 	nuc	= [BactName,' nuc'];
% Comments= {BactName,BactName,cyto,cyto,cyto,cyto,[cyto,'Int'],...
% 	nuc,nuc,nuc,nuc,[nuc,'Int'],...
% 	[BactName,' Rg'],[BactName,' Rg_Long'],[BactName,' Rg_Lat'],...
% 	[BactName,' Rg_Major'],[BactName,' Rg_Minor'],...
% 	cyto,cyto,cyto,...
% 	cyto,cyto,...
% 	nuc,nuc,nuc,...
% 	nuc,nuc,nuc,...
% 	cyto,nuc,...
% 	[BactName,' RgL/L'],[BactName,' RgM/M'],...
% 	[BactName,' Length'],[BactName,' Width'],...
% 	[BactName,' Volume'],[BactName,' VolInt'],...
% 	[BactName,' I0'],[BactName,' Mean intensity'],[BactName,' Mean/Total intensity']...
% 	cyto,...
% 	cyto,cyto,[cyto,'Int'],...
% 	nuc};
% 
% Rwrite.Properties.UserData.VarNames		= LongName;
% Rwrite.Properties.VariableUnits			= Units;
% Rwrite.Properties.VariableDescriptions	= Comments;
% Data	= table2array(Rwrite);
%%
% save(RName,'Experiment','Rtab','Rwrite');
% XlswriteTable(Rwrite,xlsName,BactName);
% xlswrite(xlsName,LongName,BactName,'A1');
% xlswrite(xlsName,Units,BactName,'A2');
% xlswrite(xlsName,Comments,BactName,'A3');
% xlswrite(xlsName,Data,BactName,'A4');
