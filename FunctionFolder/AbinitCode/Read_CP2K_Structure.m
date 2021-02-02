function [COORDS, LATT] = Read_CP2K_Structure()



handle = fopen('USPEX-pos-1.xyz');
while 1
    try
        tmp = fgetl(handle); % number of atoms
		NCoords = str2num(tmp);
		coords  = zeros(NCoords, 3);
        tmp = fgetl(handle); % string with energy
        for aa = 1 : NCoords
            at_tmp = fgetl(handle); % coordinates strings
            at_1 = sscanf(at_tmp, '%*s %g %g %g');
            coords(aa,1) = at_1(1); coords(aa,2) = at_1(2); coords(aa,3) = at_1(3);
        end
    catch
        break;
    end
end
status = fclose(handle);

[s1, s2] = unix('tail -n 1 USPEX-1.cell');
lat1 = str2num(s2);
lat = zeros(3,3);
for i = 1 : 3
    for j = 1 : 3
        lat(i,j) = lat1(2 + (i-1)*3 + j);
    end
end

%% end of file reading 


% optimize lattice
[coord1,lat] = optLattice(coords,lat);
coords = coord1/lat;

coords = coords - floor(coords);



LATT   = lat;
COORDS = coords;

