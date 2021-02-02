function [ isAllConnected, isConnectedToSubstrate ] = surface_connectivity_check( lattice, coord, atom_types, chanAList) 
% isAllConnected: true or false
% isConnectedToSubstrate: an array of flags, true for atoms connected to
% substrate, false for atoms not connected to substrate.

% tolerance is the multiplier to the covalent radii
tolerance = 1.6 ;
natoms = length(atom_types);
isConnectedToSubstrate = zeros(natoms, 1) ;
isNew = ones(natoms, 1) ;

for i=1:natoms
   if chanAList(i)==0
      isConnectedToSubstrate(i) = 1;
   end
end
[ID, rank] = sort(coord(:,3));  %sort from low to high
for i = 1:natoms
    ID1 = rank(i);
    radii1 = str2num(covalentRadius(ceil(atom_types(ID1))));
    if isConnectedToSubstrate(ID1) == 0
        for j = 1:natoms
           ID2 = rank(j);
           radii2 = str2num(covalentRadius(ceil(atom_types(ID2))));
           if distance(lattice, coord(ID1,:), coord(ID2,:)) <= tolerance * (radii1+radii2)
              isConnectedToSubstrate(ID1) = 1 ;
              break;
           end
        end
    end
end

if numel( find( isConnectedToSubstrate == 0 ) ) > 0 
    isAllConnected = false ;
else
    isAllConnected = true ;
end

function dist=distance(lattice, coord1, coord2)

check = coord1-coord2;
for index=1:3
   check(index)=check(index)-round(check(index));
end
dist = norm(check*lattice); 

