function imout = PsfGenImg(C,Ix,Iy,calibrf,lem)
% imout = PsfGenImg(C,Ix,Iy,calibrf,lem)
% Generate image based on a given set of coordinates.
% C is in the unit of um
% Ix and Iy are the dimension of imout
% calibrf is in the unit of um/px
% lem is the emmision wavlengh for the PSF.
% DY180625
%%
if ~exist('lem','var')
	lem = 610;
end

im0		= zeros(Ix,Iy);
lC		= size(C,1);
imout	= im0;
for i = 1:lC
	imout	= imout + PSFgen(im0,calibrf,C(i,:),lem);
end
imout	= imout/sum(sum(imout));