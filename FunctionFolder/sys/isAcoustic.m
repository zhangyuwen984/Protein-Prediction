function ac = isAcoustic(eigVector, coord, superCell)

% we assume that mode is accoustic, if all atoms within the basic cell move in the same direction

coord = coord - floor(coord); % transform 1.0 into 0.0
N = size(coord, 1);
eigVector = round(eigVector*10^12)/10^12;   % to make machine zeros like 10^-15 real zeros


% determine the base cells for atoms
cells = zeros(N,3);
for a = 1 : N
 for i = 1 : superCell(1)
  if sum(cells(a,:)) > 0 
    break
  end
  for j = 1 : superCell(2)
   if sum(cells(a,:)) > 0 
     break
   end
   for k = 1 : superCell(3)
     goodX = (coord(a,1) < i/superCell(1)) & (coord(a,1) >= (i-1)/superCell(1));
     goodY = (coord(a,2) < j/superCell(2)) & (coord(a,2) >= (j-1)/superCell(2));
     goodZ = (coord(a,3) < k/superCell(3)) & (coord(a,3) >= (k-1)/superCell(3));
     if goodX & goodY & goodZ 
       cells(a,1) = i;
       cells(a,2) = j;
       cells(a,3) = k;
       break;
     end
   end
  end
 end
end

visited = zeros(superCell(1), superCell(2), superCell(3)); % what basic cells were visited
inPhase = 1;
for i = 1 : N
 a1 = eigVector((i-1)*3+1 :(i-1)*3+3); % displacement for atom i
 if visited(cells(i,1), cells(i,2), cells(i,3))   % basic cell already visited
    continue
 end
 for j = 1 : N
  if norm(cells(i,:)-cells(j,:)) > 0  % different basic cells
    continue
  end
  a2 = eigVector((j-1)*3+1 :(j-1)*3+3); % displacement for atom j
  isZero = ((norm(a1) == 0) | (norm(a2)==0));
  s = dot(a1, a2);
  if (sign(s) < 0) | ((sign(s) == 0) & (isZero == 0))
    inPhase = 0;
    break
  end
 end
 visited(cells(i,1), cells(i,2), cells(i,3)) = 1; % basic cell already visited
end
if inPhase
%  disp('Seems like accoustical');
  ac = 1;
else
  ac = 0;
end

