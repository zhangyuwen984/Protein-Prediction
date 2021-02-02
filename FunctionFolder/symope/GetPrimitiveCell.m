function [numIons, lat1] = GetPrimitiveCell(sGroup, lat1, numIons, fixLat);

% apparently Stokes code reads conventional lattice and returns non-conventional...
if (sGroup(1) ~= 'P') & (fixLat == 1) % non primitive - have to make real unit cell from primitive
    if sGroup(1) == 'F'
        numIons = round(numIons/4); % primitive cell, built by Stokes code, uses only 1/4 of atoms
        lat1(1:3) = lat1(1:3)/power(4,1/3);
    elseif sGroup(1) == 'R'
        [doedl, IX] = sort(lat1(4:6));
        if (abs(doedl(1)-90) < 1) & (abs(doedl(2)-90) < 1) & (abs(doedl(3)-120) < 1)  % hexagonal unit cell
            hexagonal_Rcell = 1;
            numIons = round(numIons/3); % primitive cell, built by Stokes code, uses only 1/3 of atoms
            lat1(1:3) = lat1(1:3)/power(3,1/3);
        else
            hexagonal_Rcell = 0;   % rhombohedral unit cell, already primitive. But Stokes wants the input in the form of hexagonal cell
            if (size(lat1,1) == 3) & (size(lat1,2) == 3)
                lat2 = lat1;
            else
                lat2 = latConverter(lat1);
            end
            lat2 = [1 -1 0; 0 1 -1; 1 1 1]*lat2; % make a hexagonal unit cell out of primitive rhombohedral
            lat1 = latConverter(lat2);
            lat1(4:6) = lat1(4:6)*180/pi; % fortran works with degrees
        end
    else
        numIons = round(numIons/2); % primitive cell, built by Stokes code, uses only 1/2 of atoms
        lat1(1:3) = lat1(1:3)/power(2,1/3);
    end
end
