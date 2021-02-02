function [lat_conv, coor_conv] = unitCellFromPrimitive(cellType, lat_prim, coor_prim)
% makes real unit cell from primitive
% cellType = A, B, C, F, I, R

if cellType == 'F' % face centered
    U = [-1 +1 +1; +1 -1 +1; +1 +1 -1];
    Duplicate = [0 0 0; 0.5 0 0.5; 0 0.5 0.5; 0.5 0.5 0];
    
elseif cellType == 'I' % body centered
    U = [0 1 1; 1 0 1; 1 1 0];
    Duplicate = [0 0 0; 0.5 0.5 0.5];
    
elseif cellType == 'C'
    U = [1 1 0; 1 -1 0; 0 0 1];
    Duplicate = [0 0 0; 0.5 0.5 0];
    
elseif cellType == 'B' % groups number 38-41
    U = [0 0 1; 1 -1 0; 1 1 0];
    Duplicate = [0 0 0; 0 0.5 0.5];
    
elseif cellType == 'A'
    U = [0 1 -1; 0 1 1; 1 0 0];
    Duplicate = [0 0 0; 0.5 0.5 0];
    
elseif cellType == 'R' % rhombohedral
    U = [1 -1 0; 0 1 -1; 1 1 1];
    Duplicate = [0 0 0; -1/3 1/3 1/3; -2/3 2/3 2/3];
end

N_atom = size(coor_prim,1);
N_Dupl = size(Duplicate, 1);
lat_conv = U*lat_prim;    % lattice in general (non triangular) representation
coor_prim = coor_prim*lat_prim/lat_conv;  %fractional coordinates in lat_conventional
coor_conv = zeros(N_atom*N_Dupl, 3);

%We must add the atoms one by one (e.g., 4Cu+4O -> 8Cu+8O),
for i = 1: N_atom
    for j = 1: N_Dupl
        coor_conv(N_Dupl*(i-1)+j,:) = coor_prim(i,:) + Duplicate(j,:);
    end
end

coor_conv = coor_conv - floor(coor_conv);
lat_conv = latConverter(latConverter(lat_conv)); % back to real lattice in triangular form
