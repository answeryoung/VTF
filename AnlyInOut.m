function Simul	= AnlyInOut(Simul,Anly)
% DY190410
%%
% N	= size(Simul.Atom,1);
% nF	= Simul.TotalTimeSteps;
% Simul.RgMin	= sqrt(Simul.RgMi1.*Simul.RgMi2);
% InOut		= zeros(N,nF);
% 
% Cx	= Simul.Coords(:,1,:) - Anly.cx;
% Cy	= Simul.Coords(:,2,:) - Anly.cy;
% Cz	= Simul.Coords(:,3,:) - Anly.cz;
% 
% x0	= Simul.CoM(:,1);	
% y0	= Simul.CoM(:,2);	
% z0	= Simul.CoM(:,3);	
% 
% Rma = Simul.RgMaj;
% Rmi = Simul.RgMin;
% 
% for f = 1:nF	
% 	InOut(:,f)	= ((Cz(:,:,f) - z0(f,:)).^2)/Rma(f,:)^2 + ...
% 		((Cx(:,:,f) - x0(f,:)).^2 + (Cy(:,:,f) - y0(f,:)).^2)/Rmi(f,:);
% end
% Simul.InOut	= InOut;
%%
N	= size(Simul.Atom,1);
nF	= Simul.TotalTimeSteps;
Simul.RgMin	= sqrt(Simul.RgMi1.*Simul.RgMi2);
Simul.RgLat	= sqrt(Simul.RgLatX.*Simul.RgLatY);
InOut		= zeros(N,nF);

Cx	= Simul.Coords(:,1,:) - Anly.cx;
Cy	= Simul.Coords(:,2,:) - Anly.cy;
Cz	= Simul.Coords(:,3,:) - Anly.cz;

x0	= Simul.CoM(:,1);	
y0	= Simul.CoM(:,2);	
z0	= Simul.CoM(:,3);	

Rma = Simul.RgLong;
Rmi = Simul.RgLat;

for f = 1:nF	
	InOut(:,f)	= ((Cz(:,:,f) - z0(f,:)).^2)/Rma(f,:)^2 + ...
		((Cx(:,:,f) - x0(f,:)).^2 + (Cy(:,:,f) - y0(f,:)).^2)/Rmi(f,:);
end
%%
Simul.VnLL	= 4/3*pi*Rma.*Rmi.^2;
Simul.VnMM	= 4/3*pi*Simul.RgMaj.*Simul.RgMin.^2;
Simul.InOut	= InOut;