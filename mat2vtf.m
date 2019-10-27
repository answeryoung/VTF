function mat2vtf(mFile,vtf_fn)
% Rewrite .vtf files from Simul in .mat files
% Coordinates are not in original percision
% 
% DY170327
% DY191026
%%
load(mFile,'Simul');
if ~exist('vtf_fn','var') || isempty(vtf_fn)
	if ~strcmp(mFile(end-3:end),'.mat')
		SimulName	= mFile;
	else
		SimulName	= mFile(1:end-4);
	end
	vtf_fn	= [SimulName,'.vtf'];
elseif ~strcmp(vtf_fn(end-3:end),'.vtf')
	SimulName	= vtf_fn;
	vtf_fn		= [vtf_fn,'.vtf'];
else
	SimulName = vtf_fn(1:end-4);
end

%% Checks
if ~strcmp(Simul.Name,SimulName)
	disp('Simulation names do not aggree. Using file name...');
end
if exist(vtf_fn,'file')
	copyfile(vtf_fn,[SimulName,'S.vtf']);
end

%% Write header
vtfID		= fopen(vtf_fn,'w');
LineNumber	= 1;
if isfield(Simul,'head')
	l_head = length(Simul.head);
	for j = 1:lhead
		fprintf(vtfID,'%s\n', Simul.head{j});
	end
	LineNumber = LineNumber + l_head;
end

if isfield(Simul,'PBC')
	fprintf(vtfID,'pbc %f %f %f\n',Simul.PBC);
	LineNumber = LineNumber + 1;
end

if isfield(Simul,'Atom')
	AtomTypes	= unique(Simul.Atom.type);
	nTyp		= length(AtomTypes);
	for typ = 1:nTyp
		idx = find(Simul.Atom.type == AtomTypes(typ));
		idx0= idx(1);
		idx1= idx(end);
		AtomRadius	= Simul.Atom.radius(idx0);
		AtomName	= Simul.Atom.name{idx0};
		AtomQ		= Simul.Atom.q(idx0);
		fprintf(vtfID,'atom %d:%d radius %f name %s type %d q %f\n',...
			idx0,idx1,AtomRadius,AtomName,typ,AtomQ);
	end
	LineNumber	= LineNumber + nTyp;
	
	if ismember('bond',Simul.Atom.Properties.VariableNames)...
			&& isfield(Simul,'Bond')
		nBond	= size(Simul.Bond,1);
		for BondIdx = 1:nBond
			fprintf(vtfID,'bond %d:%d\n',Simul.Bond(BondIdx,:));
		end
		LineNumber	= LineNumber + nBond;
	end
end
fprintf(vtfID,'\n');
LineNumber	= LineNumber + 1;

%% Write coordinates timestep by timestep
for t = 1:Simul.TotalTimeSteps
	fprintf(vtfID,'%s\n','timestep ordered');
	Coords	= Simul.Coords(:,:,t);
	Coords	= Coords(~isnan(Coords(:,1)),:);
	nPart	= size(Coords,1);
	for p = 1:nPart
		fprintf(vtfID,'%f %f %f\n',Coords(p,:));
	end	
	fprintf(vtfID,'\n');
	LineNumber	= LineNumber + nPart + 2;
end

%% Close
fclose(vtfID);
disp(['Written ',int2str(LineNumber),' lines to ',vtf_fn]);
end