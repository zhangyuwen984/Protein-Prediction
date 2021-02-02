function connected = checkConnectivity(coords, lat, composition)

% USPEX 8.2.7
% checks if the cluster is connected (threshold - good_bonds)

global ORG_STRUC

N = size(coords,1);
vect = zeros(1,3);

R_val = zeros(1,length(ORG_STRUC.atomType));
for i = 1 : length(ORG_STRUC.atomType)
 s = covalentRadius(ORG_STRUC.atomType(i));
 R_val(i) = str2num(s);
end

atomTypes = zeros(1,N);
color = zeros(1,N);
for i = 1 : N
  color(i) = i;
  tmp = i;
  while tmp > 0 
   atomTypes(i) = atomTypes(i) + 1;
   tmp = tmp - composition(atomTypes(i));
  end
end 

for i = 1 : N-1
 for j = i+1 : N
  vect(1) = coords(i,1) - coords(j,1); 
  vect(2) = coords(i,2) - coords(j,2);
  vect(3) = coords(i,3) - coords(j,3);
  delta = sqrt(sum((vect*lat).^2)) - R_val(atomTypes(i)) - R_val(atomTypes(j)); 
  if delta <= -0.37*log(ORG_STRUC.goodBonds(atomTypes(i), atomTypes(j)))    % atoms i and j - connected
   c = color(i);
   for k = 1 : N
    if color(k) == c
     color(k) = color(j);
    end
   end
  end
 end
end

connected = 1;
for i = 1 : N
 if color(i) ~= color(1)
  connected = 0;
 end
end