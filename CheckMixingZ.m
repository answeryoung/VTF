function O = CheckMixingZ(CoordsA,CoordsB)
% function O = CheckMixingZ(CoordsA,CoordsB)
% O = true	if mixed
% O = false	if not
% DY190811
%%
Az	= permute(CoordsA(:,3,:),[1,3,2]);
Bz	= permute(CoordsB(:,3,:),[1,3,2]);

Amin = min(Az);
Bmax = max(Bz);
if ~isempty(find(Amin-Bmax>=0,1))
	O = false;
	return
end

Bmin = min(Bz);
Amax = max(Az);
if ~isempty(find(Bmin-Amax>=0,1))
	O = false;
else
	O = true;
end
