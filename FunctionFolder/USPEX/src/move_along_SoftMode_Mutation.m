function [lat, new_Coord, deviation]= move_along_SoftMode_Mutation(coord, numIons, lattice, eigvector, mut_degree)

% USPEX Version 8.0.0
% atom mutation based on soft modes: all atoms moved along eigenvector corresponding to the softest mode

global ORG_STRUC
global POP_STRUC

N = sum(numIons);
vec = zeros(1,3);
new_Coord = coord;
lat = lattice;

close_enough = -0.37*log(ORG_STRUC.goodBonds);        % this part needed for clusters
vect = zeros(1,3);
R_val = zeros(1,length(ORG_STRUC.atomType));
for i = 1 : length(ORG_STRUC.atomType)
 s = covalentRadius(ceil(ORG_STRUC.atomType(i)));
 R_val(i) = str2num(s);
end
if ORG_STRUC.dimension==0 %cluster == 1
  coord = coord + 0.0001;
  [lat, candidate] = makeCluster(lattice, coord, ORG_STRUC.vacuumSize(1));
  [candidates] = moveCluster(lat, candidate);
  lattice = lat;
  coord = candidates;
end
at_types = zeros(1,N);
for k = 1 : N
   tmp = k;
   while tmp > 0
     at_types(k) = at_types(k) + 1;
     tmp = tmp - numIons(at_types(k));
   end
end                                                  % END this part needed for clusters

if length(lattice) == 6
 lattice = latConverter(lattice);
end

coef = 0;
for i = 1 : N
  vec(1) = eigvector((i-1)*3+1);
  vec(2) = eigvector((i-1)*3+2);
  vec(3) = eigvector((i-1)*3+3);
  if norm(vec) > coef
   coef = norm(vec);
  end
end

normfac = mut_degree*ORG_STRUC.howManyMut/coef;   % this way the maximal displacement for an atom is equal to howManyMut*mut_degree

%normfac = mut_degree/coef; 
%normfac = mut_degree*sqrt(N); % this way mut_degree = sqrt([L1^2 + L2^2 + ... +Ln^2]/n) where Li - displacement for atom i


notDone = 1;
step = 0;

while notDone

 for i = 1 : N                      % shift the atoms
   vec(1) = eigvector((i-1)*3+1);
   vec(2) = eigvector((i-1)*3+2);
   vec(3) = eigvector((i-1)*3+3);

%  new_Coord(i,:) = coord(i,:) + (norm(vec)*normfac)*vec*inv(lattice);
   new_Coord(i,:) = coord(i,:) + normfac*vec*inv(lattice);
   %new_Coord(i,1) = new_Coord(i,1) - floor(new_Coord(i,1));
   %new_Coord(i,2) = new_Coord(i,2) - floor(new_Coord(i,2));
   %new_Coord(i,3) = new_Coord(i,3) - floor(new_Coord(i,3));
   deviation(i) = norm(normfac*vec);
 end

 if ORG_STRUC.dimension==0 %cluster == 1    % if atoms move away form cluster - mutate them less
  badAtoms = ones(1,N);
  dist = 10*ones(N,N);
  min_dist = zeros(1,N);
  for i = 1 : N
   for j = i+1 : N 
     vect(1) = new_Coord(i,1) - new_Coord(j,1);
     vect(2) = new_Coord(i,2) - new_Coord(j,2);
     vect(3) = new_Coord(i,3) - new_Coord(j,3);
     delta = sqrt(sum((vect*lattice).^2)) - R_val(at_types(i)) - R_val(at_types(j));
     dist(i,j) = delta;
     dist(j,i) = delta;
     if delta < close_enough(at_types(i), at_types(j))
       badAtoms(i) = 0;
       badAtoms(j) = 0;
     end 
   end
   min_dist(i) = min(dist(i,:));
  end
  if sum(badAtoms) == 0
   notDone = 0;
  else
    step = step + 1;
    for i = 1 : N
     if badAtoms(i)
       eigvector((i-1)*3+1) = eigvector((i-1)*3+1)*0.9;
       eigvector((i-1)*3+2) = eigvector((i-1)*3+2)*0.9;
       eigvector((i-1)*3+3) = eigvector((i-1)*3+3)*0.9;
     end
    end
  end
 else
  notDone = 0;
 end

 if step > 10
%  dist = 10*ones(N,N);
%  min_dist = zeros(1,N);
%  for i = 1 : N
%   for j = i+1 : N 
%     vect(1) = coord(i,1) - coord(j,1);
%     vect(2) = coord(i,2) - coord(j,2);
%     vect(3) = coord(i,3) - coord(j,3);
%     delta = sqrt(sum((vect*lattice).^2)) - R_val(at_types(i)) - R_val(at_types(j));
%     dist(i,j) = delta;
%     dist(j,i) = delta;
%   end
%   min_dist(i) = min(dist(i,:));
%  end
%  dist
%  min_dist
%  coord
%  lattice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%unix(['echo main_structure >> POSCAR_cluster']);
%unix(['echo ' '1.0' ' >> POSCAR_cluster']);
%for latticeLoop = 1:3
%  unix(['echo ' num2str(lattice(latticeLoop,:)) ' >> POSCAR_cluster']);
%end
%unix(['echo ' num2str(numIons) ' >> POSCAR_cluster']);
%unix(['echo ' 'Direct' ' >> POSCAR_cluster']);
%for coordLoop = 1:N
%  unix(['echo ' num2str(coord(coordLoop,:)) ' >> POSCAR_cluster']);
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quit

  notDone = 0;
 end

end
