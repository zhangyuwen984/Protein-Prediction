function Write_ATK(Ind_No)

% Change: created

global POP_STRUC
global ORG_STRUC

Step    = POP_STRUC.POPULATION(Ind_No).Step;
numIons = POP_STRUC.POPULATION(Ind_No).numIons;
lat     = POP_STRUC.POPULATION(Ind_No).LATTICE;
coord   = POP_STRUC.POPULATION(Ind_No).COORDINATES*POP_STRUC.POPULATION(Ind_No).LATTICE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ATK.in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nothing, nothing] = unix('echo "# Set up lattice" > ATK.in');
[nothing, nothing] = unix(['echo vector_a = [' num2str(lat(1,1)) ', ' num2str(lat(1,2)) ', ' num2str(lat(1,3)) ']*Angstrom >> ATK.in']);
[nothing, nothing] = unix(['echo vector_b = [' num2str(lat(2,1)) ', ' num2str(lat(2,2)) ', ' num2str(lat(2,3)) ']*Angstrom >> ATK.in']);
[nothing, nothing] = unix(['echo vector_c = [' num2str(lat(3,1)) ', ' num2str(lat(3,2)) ', ' num2str(lat(3,3)) ']*Angstrom >> ATK.in']);
[nothing, nothing] = unix(['echo "lattice = UnitCell(vector_a, vector_b, vector_c)" >> ATK.in']);
[nothing, nothing] = unix('echo "# Define elements" >> ATK.in');
s = [];
for i = 1 : length(numIons)
 for j = 1 : numIons(i)
  if isempty(s)
    s = elementFullName(ORG_STRUC.atomType(i));
   else
   s = [s ', ' elementFullName(ORG_STRUC.atomType(i))];
  end
 end
end
[nothing, nothing] = unix(['echo "elements = [ ' s ']" >> ATK.in']);

[nothing, nothing] = unix('echo "# Define coordinates" >> ATK.in');
for i = 1 : sum(numIons)
  prefixS = '';
  suffixS = ',';
  if i == 1
    prefixS = 'cartesian_coordinates = [';
  end
  if i == sum(numIons)
    suffixS = ']*Angstrom';
  end
  [nothing, nothing] = unix(['echo ' prefixS '[' num2str(coord(i,1)) ', ' num2str(coord(i,2)) ', ' num2str(coord(i,3)) ']' suffixS ' >> ATK.in']);
end

[nothing, nothing] = unix(['echo "# Set up configuration" >> ATK.in']);
[nothing, nothing] = unix(['echo "bulk_configuration = BulkConfiguration(" >> ATK.in']);
[nothing, nothing] = unix(['echo "bravais_lattice=lattice," >> ATK.in']);
[nothing, nothing] = unix(['echo "elements=elements," >> ATK.in']);
[nothing, nothing] = unix(['echo "cartesian_coordinates=cartesian_coordinates" >> ATK.in']);
[nothing, nothing] = unix(['echo ")" >> ATK.in']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KPOINTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Kpoints, Error] = Kgrid(lat, ORG_STRUC.Kresol(Step), ORG_STRUC.dimension);
if Error == 1  % This LATTICE is extremely wrong, let's skip it from now
    POP_STRUC.POPULATION(Ind_No).Error = POP_STRUC.POPULATION(Ind_No).Error + 4;
else
    POP_STRUC.POPULATION(Ind_No).K_POINTS(Step,:)=Kpoints;
end

try
 handle = fopen(['ATK_input_' num2str(Step)]);
 handle1 = fopen('ATK.in', 'at'); 
 while 1
   tmp = fgetl(handle); % read the line in the file
   if tmp == -1 % End Of File
     break;
   end
   tmp = strrep(tmp,'%','%%');
   if isempty(findstr(tmp, 'k_point_sampling'))
     fprintf(handle1, [tmp '\n']);
   else
     fprintf(handle1, ['k_point_sampling=(' num2str(Kpoints(1,1)) ', ' num2str(Kpoints(1,2)) ', ' num2str(Kpoints(1,3)) '),\n']);
   end 
 end
 status = fclose(handle);
 status = fclose(handle1);
catch
    cd (ORG_STRUC.homePath)
    [nothing, nothing] = unix(['echo "ATK input is not present for Step ' num2str(Step) '" >> ERROR_ATK']);
    quit
end
