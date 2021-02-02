function [sgrp_num] = findsym_stokes(lat, coord, numIons, atomType, tolerance)

global ORG_STRUC

spgBINDIR=[ORG_STRUC.USPEXPath '/FunctionFolder/spacegroup'];

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
fp = fopen('sym.in','w');
fprintf(fp, 'testing\n');
fprintf(fp, '%6.4f accuracy\n', tolerance);
fprintf(fp, '2  form of lattice parameters: to be entered as lengths and angles\n');
fprintf(fp, '%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f a,b,c,alpha,beta,gamma\n', lat(1:6));
fprintf(fp, '1  form of primitive lattice vectors\n');
fprintf(fp, '1 0 0 primitive lattice vectors\n');
fprintf(fp, '0 1 0\n');
fprintf(fp, '0 0 1\n');
fprintf(fp, '%4d number of atoms in the primitive unit cell\n', sum(numIons));

for ii = 1 : length(numIons)
    for jj = 1 : numIons(ii)
        fprintf(fp, '%2d',ii);
    end
end
fprintf(fp, ' type of each atom\n');

for coordLoop = 1 : size(coord,1)
    fprintf(fp, '%12.6f %12.6f %12.6f\n', coord(coordLoop,:));
end
fclose(fp);

[nothing, nothing] = unix([spgBINDIR '/findsym_new < sym.in > sym.out']);
[nothing, sgrp_num] = unix('grep "_symmetry_Int_Tables_number" sym.out | cut -d" " -f2');
sgrp_num = deblank(sgrp_num);
