function [symmetrized, newCoords, newNumIons, newLat] = symmetrize(lat, coord, numIons, atomType, tolerance)

% symmetrizes the structure using the Stokes SG determination code

% tolerance = 0.02;


global ORG_STRUC


getPy    =[ORG_STRUC.USPEXPath '/FunctionFolder/getInput.py'];
spgBINDIR=[ORG_STRUC.USPEXPath '/FunctionFolder/spacegroup'];


symmetrized = 0; % if something goes bad - no symmetrization is performed

i = 1;
while i <= length(numIons)
 if numIons(i) == 0
   numIons(i) = [];
   atomType(i) = [];
   i = i - 1;
 end
 i = i + 1;
end

if (size(lat,1) == 3) & (size(lat,2) == 3)
  lat = latConverter(lat);
end

lat(4:6) = lat(4:6)*180/pi; % Stokes code works with degrees

lattice = [num2str(lat(1), 10) ' ' num2str(lat(2), 10) ' ' num2str(lat(3), 10) ' ' num2str(lat(4), 10) ' ' num2str(lat(5), 10) ' ' num2str(lat(6), 10) ' ' ];

[a,b]=unix(['echo symmetrization > sym.in']); % title string
[a,b]=unix(['echo ' num2str(tolerance,10) ' accuracy >> sym.in']);
[a,b]=unix(['echo 2    form of lattice parameters: to be entered as lengths and angles >> sym.in']);
[a,b]=unix(['echo ' lattice '   a,b,c,alpha,beta,gamma >> sym.in']);
[a,b]=unix(['echo 1      form of primitive lattice vectors >> sym.in']);
[a,b]=unix(['echo 1 0 0  primitive lattice vectors >> sym.in']);
[a,b]=unix(['echo 0 1 0 >> sym.in']);
[a,b]=unix(['echo 0 0 1  >> sym.in']);
[a,b]=unix(['echo ' num2str(sum(numIons),10) ' number of atoms in the primitive unit cell >> sym.in']);

s = '';
for ii = 1 : length(numIons)
  for jj = 1 : numIons(ii)
    s = [s ' ' num2str(ii)];
  end
end

[a,b]=unix(['echo ' s ' type of each atom >> sym.in']);

for coordLoop = 1 : size(coord,1)
  [a,b]=unix(['echo ' num2str(coord(coordLoop,:),11) ' >> sym.in']);
end

[a, b] = unix([spgBINDIR '/findsym_new < sym.in > sym.out']);

%[a, sgrp_name] = unix('./getStuff sym.out Space 3');
%[tmp, error_check] = unix('./getStuff sym.out bombed 3');   % The program has bombed! message implies an error
sgrp_name  =python_uspex(getPy, '-f sym.out -b Space -c 3');
error_check=python_uspex(getPy, '-f sym.out -b bombed -c 3');  % The program has bombed! message implies an error

if (strcmp(error_check, 'has')) | isempty(sgrp_name)
  a = 1;
end

if a ~= 0 
   disp(['Error during symmetrization']);
   disp(' ');
else
% create a CIF file of the 'symmetrized' structure
  [a,b]=unix('echo end_of_file >> sym.out'); % to easily identify the end of file
  [a,b]=unix(['cat /dev/null > symmetrized.cif']);
%  unix('echo data_ findsym-output >> symmetrized.cif'); % we will write it in write_OUTPUT to have a correct structure number
  toWrite = 0;
  readNormal = 1;
  currentType = '0';
  typeA = 0;
  handle = fopen('sym.out');
 try
  while 1
    tmp = fgetl(handle);
    if readNormal == 0
      coords = sscanf(tmp,'%*s %s %g %g %g %g');
      if isempty(coords)
       break; % break the while cycle, end of file reached
      end
      if isempty(findstr(coords(1), currentType))
        typeA = typeA + 1;
        currentType = coords(1);
      end
      s = megaDoof(atomType(typeA));
      tmp = [s ' ' s ' ' num2str(coords(2),6) ' ' num2str(coords(3),6) ' ' num2str(coords(4),6) ' ' num2str(coords(5),6)];
    end
    if ~isempty(findstr(tmp, 'end_of_file'))
      break; % end the while cycle
    end
    if toWrite
      tmp = strrep(tmp,'"','\''');
      tmp = strrep(tmp,'(','\(');
      tmp = strrep(tmp,')','\)');
      [a,b]=unix(['echo ' tmp ' >> symmetrized.cif']);
    end
    if ~isempty(findstr(tmp, 'data_')) | ~isempty(findstr(tmp, 'loop_'))
     toWrite = 1 - toWrite;
    end
    if ~isempty(findstr(tmp, '_atom_site_occupancy'))
     readNormal = 0;
    end
  end
 catch
   disp(['Error during symmetrization']);
   disp(' ');
   symmetrized = 0;
 end
  status = fclose(handle);
end

if symmetrized == 0 % some error happened
  newCoords = coord;
  newNumIons = numIons;
end


