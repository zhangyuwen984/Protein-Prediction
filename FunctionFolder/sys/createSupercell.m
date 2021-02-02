function [lat, coords, numIons] = createSupercell(lat0, coords0, numIons0, i,j,k)

% creates ixjxk supercell

numIons = numIons0*i*j*k;
lat = lat0;
lat(1,:) = lat0(1,:)*i;
lat(2,:) = lat0(2,:)*j;
lat(3,:) = lat0(3,:)*k;
coords1 = coords0; % build first mini-cell
for a = 1 : N
  coords1(a,1) = coords0(a,1)/i;
  coords1(a,2) = coords0(a,2)/j;
  coords1(a,3) = coords0(a,3)/k;
end
coords = zeros(i*j*k*N,3);
for i1 = 1 : i
 for j1 = 1 : j
  for k1 = 1 : k
   for a = 1 : N
     ind = (i1-1)*j*k*N + (j1-1)*k*N + (k1-1)*N + a;
     coords(ind, 1) = coords1(a,1) + (i1-1)/i;
     coords(ind, 2) = coords1(a,2) + (j1-1)/j;
     coords(ind, 3) = coords1(a,3) + (k1-1)/k;
   end
  end
 end
end

