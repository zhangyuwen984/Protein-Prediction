function [newLattice, candidate] = reduce_Cluster(lattice,coordinates)
%To reduce from the big cell to the small cell
global ORG_STRUC

for i = 1 : length(ORG_STRUC.atomType)
    Vector(i)= str2num(covalentRadius(ceil(ORG_STRUC.atomType(i))));
end
radii = max(Vector);

AbsoluteCoord = coordinates*lattice;
ma = max(AbsoluteCoord) + radii;
mi = min(AbsoluteCoord) - radii;

newLattice = zeros(3,3);
newLattice(1,1) = ma(1) - mi(1);
newLattice(2,2) = ma(2) - mi(2);
newLattice(3,3) = ma(3) - mi(3);

candidate = bsxfun(@minus, AbsoluteCoord, mi)/newLattice; % Matrix-Vector;
