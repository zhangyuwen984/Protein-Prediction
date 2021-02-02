function [coor, lat] = Read_FHIaims_Structure()

command = ['grep "lattice_vector" geometry.in.next_step | cut -c16-62 '];
[nothing, lattice] = unix(command);

lat = [];

if ~isempty(lattice)
   lat = str2num(lattice);
end

command = ['grep "atom" geometry.in.next_step | cut -c6-53 '];
[nothing, coordinates] = unix(command);
coor = str2num(coordinates);

if isempty(lat)
   lat_len = 30;
   lat = [lat_len 0 0; 0 lat_len 0; 0 0 lat_len];
   coor = bsxfun(@minus, coor, mean(coor)); %Vectorized
   coor = coor + lat_len/2;
end

coor = coor/lat;
coor = coor - floor(coor);
