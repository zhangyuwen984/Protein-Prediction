function PROTEINS_STRUC = backbone(xyzfile, backbone_atoms)
% $Rev: 584 $
% $Author: maxim $
% $Date: 2014-08-30 10:27:47 +0400 (Sat, 30 Aug 2014) $

% Get x, y, z coordinates from input .xyz file:
backbone_crd   = backbone_xyz(backbone_atoms, xyzfile);
numIons        = size(backbone_crd,1);

min_value = min(backbone_crd);
max_value = max(backbone_crd);

lat = zeros(3,3);
for i=1:3
    lat(i,i) = max_value(i) - min_value(i);
end

norm_crd = [];
for coordLoop = 1 : numIons
    x = (backbone_crd(coordLoop,1) - min_value(1)) / lat(1,1);
    y = (backbone_crd(coordLoop,2) - min_value(2)) / lat(2,2);
    z = (backbone_crd(coordLoop,3) - min_value(3)) / lat(3,3);
    
    new_crd = [x y z];
    norm_crd = [norm_crd; new_crd];
end

lat0      = lat;
norm_crd0 = norm_crd;

% Add vacuum to all the sides and move molecule to the center of the cell:
vacuumSize = 10;
[lat, norm_crd] = makeCluster(lat0, norm_crd0, vacuumSize);
[norm_crd]      = moveCluster(lat,  norm_crd);

%-------------------------------------------------------------------------------
% Print POSCAR file:
fp = fopen('POSCAR_backbone', 'w');
fprintf(fp, 'EA0000\n');
fprintf(fp, '1.0000\n');

for latticeLoop = 1 : 3
   fprintf(fp, '%12.6f %12.6f %12.6f\n', lat(latticeLoop,:));
end

fprintf(fp, '%4s\n', 'C');
fprintf(fp, '%4d\n', numIons);
fprintf(fp, 'Direct\n');

for coordLoop = 1 : numIons
    x = norm_crd(coordLoop, 1);
    y = norm_crd(coordLoop, 2);
    z = norm_crd(coordLoop, 3);
    
    fprintf(fp, '%12.6f %12.6f %12.6f\n', [x y z]);
end
fclose(fp);
%-------------------------------------------------------------------------------

PROTEINS_STRUC.lattice           = lat;
PROTEINS_STRUC.backbone_crd      = backbone_crd;
PROTEINS_STRUC.backbone_crd_norm = norm_crd;
PROTEINS_STRUC.numIons           = numIons;
PROTEINS_STRUC.atomType          = 6;

end
