function O = CheckMixingShp(CoordsA,CoordsB)
% function O = CheckMixingShp(CoordsA,CoordsB)
% O = true	if mixed
% O = false	if not
% DY190811
%%
nT	= size(CoordsA,3);
O	= true;
for t = nT:-1:1
	shpA = alphaShape(CoordsA(:,:,t));
	shpB = alphaShape(CoordsB(:,:,t));
	if (sum(double(inShape(shpA,CoordsB(:,:,t)))) == 0)...
			&& (sum(double(inShape(shpB,CoordsA(:,:,t)))) == 0)
		O = false;
		break
	end
end
