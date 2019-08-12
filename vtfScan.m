function Simul = vtfScan(vtf_fn)
% scan vtf files and write mat file
% DY170327
% DY180713
%%
if ~strcmp(vtf_fn(end-3:end),'.vtf')
	SimName = vtf_fn;
	vtf_fn = [vtf_fn,'.vtf'];
else
	SimName = vtf_fn(1:end-4);
end
% MaxNBond = 4;
TimeStepMax = inf;
fid=fopen(vtf_fn,'r');
isTimeStep = 0;       %Find 1st timestep

LineNumber = 0;
head = cell(0,1);
write2head = 1;
while ~isTimeStep
	LineNumber = LineNumber + 1;
	tline = fgetl(fid);
	tSSP = strsplit(tline,{' ',':'});
	switch tSSP{1}
		case 'pbc'
			write2head = 0;
			pbc(1,1) = str2double(tSSP{2});
			pbc(1,2) = str2double(tSSP{3});
			pbc(1,3) = str2double(tSSP{4});
		case 'atom'
			write2head = 0;
% 			ltSSP = length(ltSSP);
			aidx(1) = str2double(tSSP{2}) + 1;
			aidx(2) = str2double(tSSP{3}) + 1;
			aidx = sort(aidx);
			f1 = tSSP{4};
			f2 = tSSP{6};
			f3 = tSSP{8};
			f1v = str2double(tSSP{5});
			f2v = tSSP{7};
			f3v = str2double(tSSP{9});
			try
				f4 = tSSP{10};
				f4v = str2double(tSSP{11});
			catch
				f4	= 'q';
				f4v = 0;
			end	
			for i = aidx(1):aidx(2)
				atom(i).(f1) = f1v;
				atom(i).(f2) = f2v;
				atom(i).(f3) = f3v;
				atom(i).(f4) = f4v;
				if ~isfield(atom(i),'bond')
					atom(i).bond = [];
				end
			end
		case 'bond'
			write2head = 0;
			ai = str2double(tSSP{2}) + 1;
			aj = str2double(tSSP{3}) + 1;
			atom(ai).bond = [atom(ai).bond,aj];
			atom(aj).bond = [atom(aj).bond,ai];
			if ~exist('bond','var')
				bond = [ai,aj];
			else
				bond = [bond;[ai,aj]];
			end
		case 'timestep'
			write2head = 0;
			isTimeStep = 1;
		otherwise
			if write2head
				head{LineNumber,1} = tline;
			end
	end
end

nAtom = size(atom,2);
TimeStep = 0;
coords = zeros(nAtom,3,1);
while TimeStep < TimeStepMax
	TimeStep = TimeStep + 1;
	coords_temp = zeros(nAtom,3);
	try
		for a = 1:nAtom
%			LineNumber = LineNumber + 1;
			tline	= fgetl(fid);
			C	= textscan(tline,'%f',3);
			coords_temp(a,:) = C{1};
		end
		coords(:,:,TimeStep) = coords_temp;
	catch
		TimeStep = TimeStep - 1;
	end

	isTimeStep = 0;
	while ~isTimeStep
% 		LineNumber = LineNumber + 1;
		tline = fgetl(fid);
		if ~isempty(tline)
			if (tline(1) >= 0)
				tSSP = strsplit(tline,{' ',':'});
				if strcmp(tSSP{1},'timestep')
					isTimeStep = 1;
				end
			else
				isTimeStep = 1;
				TimeStepMax = TimeStep;
			end
		end
	end
end
			
if ~isempty(head)
	Simul.head = head;
end
atom = struct2table(atom);
Simul.Name = SimName;
Simul.PBC = pbc;
Simul.Atom = atom;
Simul.Bond = bond;
Simul.Coords = coords;
Simul.TotalTimeSteps = TimeStep;
fclose('all');
save([SimName,'.mat'],'Simul','-v7.3');
