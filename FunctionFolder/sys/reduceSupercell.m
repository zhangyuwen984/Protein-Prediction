function [reduced, newCoord, newLat, newNumIons, newSupercell] = reduceSupercell(lat, coord, numIons, supercell, tolerance)

% tries to check whether our cell is a supercell and reduces it, if possible
% only consider supercell along the lattice vectors

maxMultiplicity = 20; % check up to this multiplicity
reduced = 0; % by default we don't reduce
N = sum(numIons);

atomType = zeros(1,N);
counter = 1;
for i = 1 : size(numIons)
  for j = 1 : numIons(i)
    atomType(counter) = i;
    counter = counter + 1;
  end
end

% for every atom we shift according to assumed multiplicity and look
% whether there is an atom within tolerance of that spot
try
 Mult = [1 1 1]; % multiplicity
 for Axis = 1 : 3 % x,y,z
  for m = 2 : maxMultiplicity % check supercells from 2 till 20 along x, y, z
   % only consider multiplicities compatible with the supercell
   if mod(supercell(Axis), m) > 0
     continue
   end
   pairFound = zeros(1,N);
   for a1 = 1 : N
     vec_a1 = coord(a1,:);
% shift atom a1 along x-axis according to multiplicity
     vec_a1(Axis) = vec_a1(Axis) + (1/m);
     vec_a1 = vec_a1 - floor(vec_a1); % in case we move beyond cell boundaries
     for a2 = 1 : N
% we transform vector into absolute units
       diff_vec = (coord(a2,:) - vec_a1)*lat;       
       if (norm(diff_vec) <= tolerance) & (atomType(a1) == atomType(a2)) 
           pairFound(a1) = 1;
           break;
       end
     end
   end
   if (sum(pairFound) == N)
     Mult(Axis) = m;  % we only care about biggest multiplicity!
   end
  end
 end

 
 if Mult(1)*Mult(2)*Mult(3) > 1
   reduced = 1;
   newLat(1,:) = lat(1,:)/Mult(1);
   newLat(2,:) = lat(2,:)/Mult(2);
   newLat(3,:) = lat(3,:)/Mult(3);
   newSupercell(1) = supercell(1)/Mult(1);
   newSupercell(2) = supercell(2)/Mult(2);
   newSupercell(3) = supercell(3)/Mult(3);
   newNumIons = numIons/(Mult(1)*Mult(2)*Mult(3));
   newCoord = zeros(sum(newNumIons),3); % for octave compatibility
   counter = 1;
   c = coord - floor(coord); % so that 1.0 become 0.0
   for a = 1 : N
     if (c(a,1) < 1/Mult(1)) & (c(a,2) < 1/Mult(2)) & (c(a,3) < 1/Mult(3))
       newCoord(counter,:) = c(a,:);
       newCoord(counter,1) = newCoord(counter,1) * Mult(1);
       newCoord(counter,2) = newCoord(counter,2) * Mult(2);
       newCoord(counter,3) = newCoord(counter,3) * Mult(3);
       counter = counter + 1;
     end
   end
% final check
   if sum(newNumIons) ~= size(newCoord,1)
       reduced = 0; % couldn't extract the atoms of the smaller cell properly
   end
 else
   newCoord = coord;
   newLat = lat;
   newNumIons = numIons;
   newSupercell = supercell;
 end
 
% don't do anything if something bad happens
catch
  reduced = 0;
  newCoord = coord;
  newLat = lat;
  newNumIons = numIons;
  newSupercell = supercell;
end