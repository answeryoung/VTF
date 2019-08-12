function Simul = SimImg(Simul,Anly,DimSet)
% DY180703
%%
if ischar(DimSet)
	DimSet = {DimSet};
end
cx	= Anly.cx;
cy	= Anly.cy;
cz	= Anly.cz;
Ix	= Anly.SimImgX;
Iy	= Anly.SimImgY;
calibrf = Anly.calibrf;
Nim = length(Anly.Img);

C(:,1,:)= Simul.Coords(:,1,:) - cx;
C(:,2,:)= Simul.Coords(:,2,:) - cy;
C(:,3,:)= Simul.Coords(:,3,:) - cz;

types = unique(Simul.Atom.type)' + 1;
for tp = types
	pidx{tp,1}	= find(Simul.Atom.type == tp-1);
end

imidx = 0;
Simul.Image = zeros(Ix,Iy,max(types),Nim*length(DimSet));
for Dim = DimSet
	switch Dim{1}
		case 'TwoDim'
			Ct	= C(:,Anly.TwoDim,:)*Anly.Unit;
			Nd	= 2;
		case 'TriDim'
			Ct	= C(:,Anly.TriDim,:)*Anly.Unit;
			Nd	= 3;
	end

	for im = 1:Nim
		imidx = imidx + 1;
		for tp = types
			switch Anly.Img(im).Color(tp)
				case 1 % TagRFP-T
					lem = Anly.RedEm;
				case 2 % mNeonGreen
					lem = Anly.GreenEm;
				case 3 % TagRFP-T
					lem = Anly.RedEm;
			end
			Ctt	= Ct(pidx{tp},:,Anly.Img(im).t0:Anly.Img(im).t1);
			Ctt	= reshape(permute(Ctt,[1,3,2]),[],Nd);
			Simul.Image(:,:,tp,imidx)	= PsfGenImg(Ctt,Ix,Iy,calibrf,lem);
		end
	end
end