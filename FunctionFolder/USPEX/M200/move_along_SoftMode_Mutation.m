function [lat, new_Coord, deviation]= move_along_SoftMode_Mutation(coord, numIons, lattice, eigvector, mut_degree)
%all atoms moved along eigenvector corresponding to the softest mode

global ORG_STRUC

N = sum(numIons);
vec = zeros(1,3);
new_Coord = coord;

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

for i = 1 : N                      % shift the atoms
  vec(1) = eigvector((i-1)*3+1);
  vec(2) = eigvector((i-1)*3+2);
  vec(3) = eigvector((i-1)*3+3);
  new_Coord(i,:) = coord(i,:) + (norm(vec)*normfac)*vec*inv(lattice);
  new_Coord(i,:) = coord(i,:) + normfac*vec*inv(lattice);
  new_Coord(i,1) = new_Coord(i,1) - floor(new_Coord(i,1));
  new_Coord(i,2) = new_Coord(i,2) - floor(new_Coord(i,2));
  new_Coord(i,3) = new_Coord(i,3) - floor(new_Coord(i,3));
  deviation(i) = norm(normfac*vec);
end
