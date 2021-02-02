function [eigVector, coords, lat, numIons] = calcEigenvectorK(eigVectorK, supercell, lat0, coords0, numIons0)

% creates the proper eigenvector out of 'k' one for varcomp

N = size(coords0, 1);
k1 = supercell(1);
k2 = supercell(2);
k3 = supercell(3);
kVector = zeros(1,3);
if k1 > 1
   kVector(1) = 1/k1;
end
if k2 > 1
   kVector(2) = 1/k2;
end
if k3 > 1
   kVector(3) = 1/k3;
end

numIons = numIons0*k1*k2*k3;
lat(1,:) = lat0(1,:)*k1;
lat(2,:) = lat0(2,:)*k2;
lat(3,:) = lat0(3,:)*k3;
coords1 = coords0; % build first mini-cell
for a = 1 : N
  coords1(a,1) = coords0(a,1)/k1;
  coords1(a,2) = coords0(a,2)/k2;
  coords1(a,3) = coords0(a,3)/k3;
end
coords = zeros(k1*k2*k3*N,3);
for i = 1 : k1
 for j = 1 : k2
  for k = 1 : k3
   for a = 1 : N
     ind = (i-1)*k2*k3*N + (j-1)*k3*N + (k-1)*N + a;
     coords(ind, 1) = coords1(a,1) + (i-1)/k1;
     coords(ind, 2) = coords1(a,2) + (j-1)/k2;
     coords(ind, 3) = coords1(a,3) + (k-1)/k3;
   end
  end
 end
end

% make real eigenvectors (or rather displacements) out of 'k' ones
% it looks like we can simply take the sign of the real part
rec_lat = zeros(3,3);
rec_lat(1,:) = 2*pi*cross(lat0(2,:), lat0(3,:))/det(lat0);
rec_lat(2,:) = 2*pi*cross(lat0(3,:), lat0(1,:))/det(lat0);
rec_lat(3,:) = 2*pi*cross(lat0(1,:), lat0(2,:))/det(lat0);
kVectorA = kVector*rec_lat;
coordsA = coords*lat;
eigVector = zeros(3*N*k1*k2*k3, 3*N);
col = 1;
for eVec = 1 : size(eigVectorK,2)
 eigenV = eigVectorK(:,eVec);
 coord_ind = 1;
 ind = 1;
 for i = 1 : 3*N*k1*k2*k3
   eigVector(i, col) = real(eigenV(ind)*(cos(dot(kVectorA, coordsA(coord_ind,:)))+sqrt(-1)*sin(dot(kVectorA, coordsA(coord_ind,:)))));
   ind = ind + 1;
   if i/3 == round(i/3)
     coord_ind = coord_ind + 1;
   end
   if ind > 3*N
     ind = 1;
   end
 end
 norm_factor = norm(eigVector(:,col));
 for i = 1 : 3*N*k1*k2*k3
   eigVector(i,col) = eigVector(i,col)/norm_factor;
 end
 col = col + 1;
end
