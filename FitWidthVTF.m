function [cytoWidthSig,cytoRawBg,CW] = FitWidthVTF(im,Coords,mLpx,k)
% cytoWidthSig = FitCytoWidth(im,Coords,mLpx,k)
% globals: ScanRange ScanStep r_fit br_lat
% DY180322
%%
global ScanRange ScanStep r_fit br_lat
if isempty(ScanRange)
	ScanRange	= 0.16;
end
if isempty(ScanStep)
	ScanStep	= 0.02;
end
if isempty(r_fit)
	r_fit	= 17;
end
if isempty(br_lat)
	br_lat	= 1;
end
%%
nWFit	= floor((mLpx-10)/4);
LCoords = size(Coords,1);
if nWFit > 0
	xset	= 0.5*(LCoords-mLpx+1) + 5 + (0:nWFit)*((mLpx-10)/nWFit);
	CW		= interp1(Coords,xset);
else
	nWFit	= 0;
	CW		= mean(Coords);
end

cytoWidthSig= zeros(1,nWFit+1);
cytoRawBg	= zeros(1,nWFit+1);
for f = 1:nWFit+1
	xb	= FanScanFit2(@Gauss,im,CW(f,:),r_fit,pi/2 + k,br_lat,ScanRange,ScanStep);
	b	= find(xb(:,3) == min(xb(:,3)));
	cytoWidthSig(1,f)	= xb(b,3);
	cytoRawBg(1,f)		= xb(b,1);
end