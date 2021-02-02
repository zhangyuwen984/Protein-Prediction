function [coord, lat] = Read_ATK_Structure()
handle = fopen('ATK.out');
while 1
  tmp = fgetl(handle);
  if ~ischar(tmp), break, end
  if ~isempty(strfind(tmp, 'Primitive vectors'))
     lat = zeros(3);
     for LL = 1 : 3
       lat_tmp = fgetl(handle); % lattice strings in angstrom units
       lat_1 = sscanf(lat_tmp,'%*s %*s %g %g %g %*s');
       lat(LL,1) = lat_1(1); lat(LL,2) = lat_1(2); lat(LL,3) = lat_1(3);
     end
  elseif ~isempty(strfind(tmp, 'Bulk: Cartesian (Angstrom)'))
      tmp = fgetl(handle); % trash string
      nA1 = fgetl(handle); % number of atoms
      nA  = sum(str2num(nA1));
      tmp = fgetl(handle); % word Bulk
      coords = zeros(nA,3);
      for aa = 1 : nA      % coordinates
          at_tmp = fgetl(handle); % coordinates strings
          at_1 = sscanf(at_tmp,'%*s %g %g %g');
          coords(aa,1) = at_1(1); coords(aa,2) = at_1(2); coords(aa,3) = at_1(3);
      end
  end
end
fclose(handle);
	  
coord=coords/lat;
coord=coord-floor(coord);
