function [candidate, newLattice, errorS] = symope_001(nsym, numIonsFull, lattice, minDistMatrix)

% creates the cluster satisfying symmetry nsym

% IMPORTANT: we test the generation of clusters in the ellipse inside the lattice. Making it more 'bulk'

ellipse_mode = 1;

%%%%%%%%%%%%%%%%%%%%%% Symmetry operations %%%%%%%%%%%%%%%%%%%
 E = [1.0 0 0; 0 1 0; 0 0 1];                                    % equivalence
 I = [-1.0, 0.0, 0.0; 0.0, -1.0, 0.0; 0.0, 0.0, -1.0];       % inversion   
 C2x = [-1 0 0; 0 1 0; 0 0 1]; % two fold axes along YOZ
 C2y = [1 0 0; 0 -1 0; 0 0 1]; % two fold axes along XOZ
 C2z = [1 0 0; 0 1 0; 0 0 -1]; % two fold axes along XOY
 Hz = [1 0 0; 0 1 0; 0 0 -1];  % reflection in the mirror plane XOY

% Cnz = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1] % rotate by Z axis 
% Cnx = [1 0 0 ; 0 cos(beta) -sin(beta); 0 sin(beta) cos(beta)]    % rotate by X axis
% Cny = [cos(gama) 0 sin(gama); 0 1 0; -sin(gama) 0 cos(gama)]     % rotate by Y axis
% For Cnv symmetry, rotations (angle = 2pi*k/n) and reflections (across a plane that makes an angle pi*k/n with the axis) have the following form:
% Rk = [cos(2pi*k/n) -sin(2pi*k/n) 0; sin(2pi*k/n) cos(2pi*k/n) 0; 0 0 1]
% Sk = [cos(2pi*k/n) sin(2pi*k/n) 0; sin(2pi*k/n) -cos(2pi*k/n) 0; 0 0 1]
% For Dn symmetry:
% Rk = [cos(2pi*k/n) -sin(2pi*k/n) 0; sin(2pi*k/n) cos(2pi*k/n) 0; 0 0 1]
% R2k = [cos(2pi*k/n) sin(2pi*k/n) 0; sin(2pi*k/n) -cos(2pi*k/n) 0; 0 0 -1]

%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% important! We assign Z as a symmetry axe and then set principal axes as coordinate axes, highest moment of inertia corresponds to Z, lowest - to X

errorS = 0;
newLattice = lattice;
volLat = det(lattice);
numIons = sum(numIonsFull);
nI = numIons;

minDistance = minDistMatrix(1,1); % CHANGE TO WHOLE ION DISTANCES MATRIX!!!!!!!!!!!!!!!!!!!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(nsym, 'E')
  candidate = rand(sum(numIons),3) - 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(nsym, 'Ci') | strcmp(nsym, 'S2')  
  candidate = [];
  while 1
    if numIons == 1;   % odd number of atoms, thus one should be in the dead center
      tmp = [0 0 0];          
      candidate = cat(1,candidate,tmp); 
      break;
    end
    if numIons == 0 
      break
    end
    tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)

    while ellipse_mode & (norm(tmp) > 0.5)
      tmp = rand(1,3) - 0.5; 
    end

    if sum(abs(tmp)) < 0.01
      continue
    end;
    nextAtom = tmp*I;
    cand = cat(1,tmp,nextAtom);

    if isempty(candidate)
      candidate = cand; 
    else
      candidate = cat(1,candidate,cand); 
    end
    numIons = numIons - 2;
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif findstr(nsym, 'I') | findstr(nsym, 'i')   % icosahedral symmetries
% easiest way to implement I - icosahedron centered at origin (2-fold axes - X, Z; 3-fold axis - (1,1,1))
% http://en.wikipedia.org/wiki/File:Icosahedron-golden-rectangles.svg
  newLattice = [volLat^(1/3) 0 0; 0 volLat^(1/3) 0; 0 0 volLat^(1/3)];
  candidate = [];
  if mod(numIons,2) == 1 % put in the center
      candidate = [0 0 0];
      numIons = numIons - 1;
  end
  if (numIons < 60) & (mod(numIons,12) ~= 0) & (numIons ~= 20) & (numIons ~= 30) & (numIons ~= 32)
   if (numIons < 40) | (numIons == 58) | (numIons == 46)
     status = ['Impossible to build the cluster with ' num2str(numIons) ' atoms that has symmetry group ' nsym];
     errorS = 1;
%    unix(['echo ' status ' > error_cluster_symmetry']);
%    quit;
   end
  end
  n60 = 0;
  numIons1 = numIons;
  while 1
    numIons1 = numIons1 - 60;
    if (numIons1 ~= 0) & (numIons1 < 60)
     if (numIons1 ~= 12) & (numIons1 ~= 20) & (numIons1 ~= 24) & (numIons1 ~= 30) & (numIons1 ~= 32) & (numIons1 ~= 36) 
      if (numIons1 < 40) | (numIons1 == 58) | (numIons1 == 46)
       break;
      end
     end
    else
      n60 = n60 + 1;
    end
  end
  GR = (sqrt(5)+1)/2;  % Golden Ratio
  alpha = atan(1/GR); % an angle to turn the axis X around Z to become 5-fold axis!
  R1 = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
  R5 = [1 0 0 ; 0 cos(2*pi/5) -sin(2*pi/5); 0 sin(2*pi/5) cos(2*pi/5)]; % turning around 5-fold X axis
  n = 1;
  while n <= n60  % do a 'full' or 'half-full' symmetry point
     tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)

     while ellipse_mode & (norm(tmp) > 0.5)
       tmp = rand(1,3) - 0.5; 
     end

%    tmp = [GR -0.33 0];   for ideal buckyball
     if (strcmp(nsym, 'Ih') | strcmp(nsym, 'ih')) % Ih (full icosahedral) symmetry
       if (mod(n60,2) == 1) & (n == n60)  
         % do only 60 poins :) (by generating a point on mirror plane) 
         tmp(3) = 0;   
       else 
         % do 120 points
         tmp = cat(1,tmp,tmp*[1 0 0; 0 1 0; 0 0 -1]);
       end
     end       
     cand = tmp;
     % rotations around the 5-fold axis (GR, -1, 0) 
     for i = 1 : 4
       nextAtom = tmp*inv(R1)*(R5^i)*R1; 
       cand = cat(1,cand,nextAtom);
     end
     % rotations around the main cube diagonal (1,1,1) (3-fold axis) 
     nextAtoms1 = cand*[0 1 0; 0 0 1; 1 0 0]; % +++ diagonal, -120 degree
     nextAtoms2 = cand*[0 0 1; 1 0 0; 0 1 0]; % -240 degree
     cand = cat(1,cand,nextAtoms1);
     cand = cat(1,cand,nextAtoms2);
     % rotations around the main axis Z (2-fold axis) 
     nextAtoms = cand*[-1 0 0; 0 -1 0; 0 0 1]; % +++ diagonal, -120 degree
     cand = cat(1,cand,nextAtoms);
     % rotations around the main axis X (2-fold axis) 
     nextAtoms = cand*[1 0 0; 0 -1 0; 0 0 -1]; % +++ diagonal, -120 degree
     cand = cat(1,cand,nextAtoms);
     if (strcmp(nsym, 'Ih') | strcmp(nsym, 'ih')) & ~((mod(n60,2) == 1) & (n == n60))  
       numIons = numIons - 120;
       n = n + 2;
     else
       numIons = numIons - 60;
       n = n + 1;
     end        
     if isempty(candidate)
      candidate = cand; 
     else
      candidate = cat(1,candidate,cand); 
     end
     if n > n60 
         break
     end
  end
  if numIons > 0  % numIons = 12x+20y+30z
    for x = 0 : 10
     for y = 0 : 6
      for z = 0 : 4
        if 12*x + 20*y + 30*z == numIons
          break
        end
      end
        if 12*x + 20*y + 30*z == numIons
          break
        end
     end
        if 12*x + 20*y + 30*z == numIons
          break
        end
    end
    for i = 1 : x % 12-atoms, thus put them on 5-fold axis
     tmp = [1 -1/GR 0]*rand*0.5;
     while ellipse_mode & (norm(tmp) > 0.5)
       tmp = [1 -1/GR 0]*rand*0.5;
     end

     cand = tmp;
     % rotations around the main cube diagonal (1,1,1) (3-fold axis) 
     cand = cat(1,cand,tmp*[0 1 0; 0 0 1; 1 0 0]); % +++ diagonal, -120 degree
     cand = cat(1,cand,tmp*[0 0 1; 1 0 0; 0 1 0]); % -240 degree
     % rotations around the main axis Z (2-fold axis) 
     nextAtoms = cand*[-1 0 0; 0 -1 0; 0 0 1]; % +++ diagonal, -120 degree
     cand = cat(1,cand,nextAtoms);
     % rotations around the main axis X (2-fold axis) 
     nextAtoms = cand*[1 0 0; 0 -1 0; 0 0 -1]; % +++ diagonal, -120 degree
     cand = cat(1,cand,nextAtoms);        
     if isempty(candidate)
      candidate = cand; 
     else
      candidate = cat(1,candidate,cand); 
     end     
    end
    for i = 1 : y % 20-atoms, thus put them on 3-fold axis
     tmp1 = rand*0.5;
     while abs(tmp1) < 0.05
      tmp1 = rand - 0.5;
     end
     tmp = [tmp1 -tmp1 tmp1];

     while ellipse_mode & (norm(tmp) > 0.5)
      tmp1 = rand*0.5;
      while abs(tmp1) < 0.05
       tmp1 = rand - 0.5;
      end
      tmp = [tmp1 -tmp1 tmp1];
     end


     cand = tmp;
     % rotations around the main axis Z (2-fold axis) 
     nextAtoms = cand*[-1 0 0; 0 -1 0; 0 0 1]; % +++ diagonal, -120 degree
     cand = cat(1,cand,nextAtoms);
     % rotations around the main axis Y (2-fold axis) 
     nextAtoms = cand*[-1 0 0; 0 1 0; 0 0 -1]; % +++ diagonal, -120 degree
     cand = cat(1,cand,nextAtoms);
     % rotations around the 5-fold axis (GR, -1, 0) 
     candTMP = cand;
     for i = 1 : 4
       nextAtom = candTMP*inv(R1)*(R5^i)*R1; 
       cand = cat(1,cand,nextAtom);
     end
     if isempty(candidate)
      candidate = cand; 
     else
      candidate = cat(1,candidate,cand); 
     end        
    end
    for i = 1 : z % 30-atoms, thus put them on 2-fold axis
     tmp1 = rand - 0.5;
     while abs(tmp1) < 0.05
      tmp1 = rand - 0.5;
     end
     tmp = [0 0 tmp1]; % Z axis, thus rotate around 2-fold X
     cand = tmp;
     % rotations around the main axis X (2-fold axis) 
     nextAtoms = cand*[1 0 0; 0 -1 0; 0 0 -1]; % +++ diagonal, -120 degree
     cand = cat(1,cand,nextAtoms); 
     % rotations around the main cube diagonal (1,1,1) (3-fold axis) 
     nextAtoms1 = cand*[0 1 0; 0 0 1; 1 0 0]; % +++ diagonal, -120 degree
     nextAtoms2 = cand*[0 0 1; 1 0 0; 0 1 0]; % -240 degree
     cand = cat(1,cand,nextAtoms1);
     cand = cat(1,cand,nextAtoms2);
     % rotations around the 5-fold axis (GR, -1, 0) 
     candTMP = cand;
     for i = 1 : 4
       nextAtom = candTMP*inv(R1)*(R5^i)*R1; 
       cand = cat(1,cand,nextAtom);
     end     
     if isempty(candidate)
      candidate = cand; 
     else
      candidate = cat(1,candidate,cand); 
     end 
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isempty(findstr(nsym, 'S')) | ~isempty(findstr(nsym, 's')) % S2n group - 2n rotoreflection (combination rotation+reflection)
  newLattice = [sqrt(lattice(1,1)*lattice(2,2)) 0 0; 0 sqrt(lattice(1,1)*lattice(2,2)) 0; 0 0 lattice(3,3)];
  candidate = [];
  rotRank = str2num(nsym(2:end));   % should be EVEN (=2n), odd is equivalent to Cnh 
  while 1
    if numIons == 1;   
      tmp = [0 0 0];          
      candidate = cat(1,candidate,tmp); 
      break;
    end
    if numIons < rotRank        % put atoms on rotation axis
      cand = [0 0 rand-0.5];
      for i = 2 : numIons
       tmp = [0 0 rand-0.5]; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
       cand = cat(1,cand,tmp); 
      end 
      numIons = 0;
    else
      tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
      while ellipse_mode & (norm(tmp) > 0.5)
        tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
      end

      cand = tmp;
      tooClose = 0;      
      for i = 1 : rotRank-1
       angle = 2*i*pi/rotRank;
       opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  % rotate around Z axis
       nextAtom = tmp*opemat*(Hz^i);
       if norm((nextAtom-tmp)*newLattice) < minDistance
           tooClose = 1;
       end
       cand = cat(1,cand,nextAtom);
      end

      if tooClose
          numIons = numIons - 1;
          cand = [0 0 tmp(3)];
      else
        numIons = numIons - rotRank;
      end
    end

    if isempty(candidate)
      candidate = cand; 
    else
      candidate = cat(1,candidate,cand); 
    end

    if numIons == 0 
      break
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(nsym, 'Th') | strcmp(nsym, 'th')   % Pyritohedral (Th) symmetry
% easiest way to implement T - include tetrahonal into the cube with O in the center and XYZ parallel to the edges
  newLattice = [volLat^(1/3) 0 0; 0 volLat^(1/3) 0; 0 0 volLat^(1/3)];
  candidate = [];
  if (numIons ~= 1) & (numIons ~= 6) & (numIons ~= 7) & (numIons ~= 8) & (numIons ~= 9) & (numIons < 12)
     status = ['Impossible to build the cluster with ' num2str(numIons) ' atoms that has symmetry group ' nsym];
     [nothing, nothing] = unix(['echo ' status ' > error_cluster_symmetry']);
     errorS = 1;
  end
  if mod(numIons,2) == 1 % put in the center
      candidate = [0 0 0];
      numIons = numIons - 1;
  end
  n12 = 0;
  numIons1 = numIons;
  while 1
    numIons1 = numIons1 - 12;
    if (numIons1 ~= 6) & (numIons1 ~= 8) & (numIons1 ~= 0) & (numIons1 < 12)
      break;
    else
      n12 = n12 + 1;
    end
  end
  n = 1;
  while n <= n12  % do a 'full' or 'half-full' symmetry point
     tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
     while ellipse_mode & (norm(tmp) > 0.5)
        tmp = rand(1,3) - 0.5;
     end

     if (mod(n12,2) == 1) & (n == n12) 
       % do only 12 poins :) (by generating a point on mirror plane) 
         tmp(ceil(rand*3)) = 0;  
     end   
     cand = tmp;
   % 2-fold rotations around X, Y, Z (O is in the cube center, Z goes up)
     nextAtom = tmp*[-1 0 0; 0 -1 0; 0 0 1]; % rotate around X axis
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[1 0 0; 0 -1 0; 0 0 -1];  % rotate around X axis
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[-1 0 0; 0 1 0; 0 0 -1];  % rotate around Y axis
     cand = cat(1,cand,nextAtom);
   % rotations around the main cube diagonals (3-fold axes) 
     nextAtom = tmp*[0 1 0; 0 0 1; 1 0 0]; % +++ diagonal, -120 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 0 1; 1 0 0; 0 1 0]; % -240 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 1 0; 0 0 -1; -1 0 0]; % ++- diagonal, 120 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 0 -1; 1 0 0; 0 -1 0]; % 240 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 -1 0; 0 0 1; -1 0 0]; % +-- diagonal, 120 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 0 -1; -1 0 0; 0 1 0]; % 240 degree     
     cand = cat(1,cand,nextAtom);     
     nextAtom = tmp*[0 -1 0; 0 0 -1; 1 0 0]; % +-+ diagonal, 120 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 0 1; -1 0 0; 0 -1 0]; % 240 degree     
     cand = cat(1,cand,nextAtom); 
     if ~((mod(n12,2) == 1) & (n == n12))
       numIons = numIons - 24;
       n = n + 2;
       cand1 = -1*cand; % add an inversion
       cand = cat(1,cand,cand1);
     else
       numIons = numIons - 12;
       n = n + 1;
     end     
     if isempty(candidate)
      candidate = cand; 
     else
      candidate = cat(1,candidate,cand); 
     end
     if n > n12 
         break
     end
  end
  if numIons > 0  % numIons can be described as 6x+8y
    n1 = ceil(numIons/8);
    n2 = floor(numIons/6);
    n = n1 + round(rand*(n2-n1)); % numIons = (4n-numIons/2)*6 + (numIons/2-3n)*8
    x = 4*n - round(numIons/2);
    y = round(numIons/2) - 3*n;
    for i = 1 : x   % 6 atoms per go - atoms on the XYZ axes
      tmp = rand - 0.5;
     while abs(tmp) < 0.05
      tmp = rand - 0.5;
     end
      candidate = cat(1, candidate, [0 0 tmp]);
      candidate = cat(1, candidate, [0 0 -1*tmp]);
      candidate = cat(1, candidate, [0 tmp 0]);
      candidate = cat(1, candidate, [0 -1*tmp 0]);
      candidate = cat(1, candidate, [tmp 0 0]);
      candidate = cat(1, candidate, [-1*tmp 0 0]);
    end
    for i = 1 : y  % 8 atoms per go - atoms on the main cube diagonals
      tmp = rand - 0.5;
      while abs(tmp) < 0.05 | (ellipse_mode & (norm([tmp tmp tmp]) > 0.5))
       tmp = rand - 0.5;
      end
      candidate = cat(1, candidate, [tmp tmp tmp]);
      candidate = cat(1, candidate, [tmp tmp -1*tmp]);
      candidate = cat(1, candidate, [tmp -1*tmp tmp]);
      candidate = cat(1, candidate, [tmp -1*tmp -1*tmp]);
      candidate = cat(1, candidate, [-1*tmp tmp tmp]);
      candidate = cat(1, candidate, [-1*tmp tmp -1*tmp]);
      candidate = cat(1, candidate, [-1*tmp -1*tmp tmp]);
      candidate = cat(1, candidate, [-1*tmp -1*tmp -1*tmp]);
    end
  end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isempty(findstr(nsym, 'T')) | ~isempty(findstr(nsym, 't'))   % Tetrahedral (T, Td) symmetries
% easiest way to implement T - include tetrahonal into the cube with O in the center and XYZ parallel to the edges
  newLattice = [volLat^(1/3) 0 0; 0 volLat^(1/3) 0; 0 0 volLat^(1/3)];
  candidate = [];
  if (numIons ~= 1) & (numIons < 4)
     status = ['Impossible to build the cluster with ' num2str(numIons) ' atoms that has symmetry group ' nsym];
     [nothing, nothing] = unix(['echo ' status ' > error_cluster_symmetry']);
     errorS = 1;
  end
  if mod(numIons,2) == 1 % put in the center
      candidate = [0 0 0];
      numIons = numIons - 1;
  end
  n12 = 0;
  numIons1 = numIons;
  while 1
    numIons1 = numIons1 - 12;
    if (numIons1 ~= 0) & (numIons1 < 4)
      break;
    else
      n12 = n12 + 1;
    end
  end
  n = 1;
  while n <= n12  % do a 'full' symmetry point
     tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
     while ellipse_mode & (norm(tmp) > 0.5)
        tmp = rand(1,3) - 0.5;
     end

     if (strcmp(nsym, 'Td') | strcmp(nsym, 'td')) & (mod(n12,2) == 1) & (n == n12)  % Td (full tetrahedral) symmetry
       % do only 12 poins :) (by generating a point on mirror plane) 
       tmp(1) = tmp(2);   
     end   
     cand = tmp;
   % 2-fold rotations around X, Y, Z (O is in the cube center, Z goes up)
     nextAtom = tmp*[-1 0 0; 0 -1 0; 0 0 1]; % rotate around X axis
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[1 0 0; 0 -1 0; 0 0 -1];  % rotate around X axis
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[-1 0 0; 0 1 0; 0 0 -1];  % rotate around Y axis
     cand = cat(1,cand,nextAtom);
   % rotations around the main cube diagonals (3-fold axes) 
     nextAtom = tmp*[0 1 0; 0 0 1; 1 0 0]; % +++ diagonal, -120 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 0 1; 1 0 0; 0 1 0]; % -240 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 1 0; 0 0 -1; -1 0 0]; % ++- diagonal, 120 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 0 -1; 1 0 0; 0 -1 0]; % 240 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 -1 0; 0 0 1; -1 0 0]; % +-- diagonal, 120 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 0 -1; -1 0 0; 0 1 0]; % 240 degree     
     cand = cat(1,cand,nextAtom);     
     nextAtom = tmp*[0 -1 0; 0 0 -1; 1 0 0]; % +-+ diagonal, 120 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 0 1; -1 0 0; 0 -1 0]; % 240 degree     
     cand = cat(1,cand,nextAtom); 
     if (strcmp(nsym, 'Td') | strcmp(nsym, 'td')) & ~((mod(n12,2) == 1) & (n == n12))  % Td (full tetrahedral) symmetry
       numIons = numIons - 24;
       n = n + 2;
       cand1 = cand; % we add a vertical reflection plane that goes through Z axis and one of the tetrahedra edges
       cand1(:,1) = -1*cand(:,2);
       cand1(:,2) = -1*cand(:,1);
       cand = cat(1,cand,cand1);
     else
       numIons = numIons - 12;
       n = n + 1;
     end     
     if isempty(candidate)
      candidate = cand; 
     else
      candidate = cat(1,candidate,cand); 
     end
     if n > n12 
         break
     end
  end
  if numIons > 0  % numIons can be described as 4x+6y
    n1 = ceil(numIons/6);
    n2 = floor(numIons/4);
    n = n1 + round(rand*(n2-n1)); % numIons = (3n-numIons/2)*4 + (numIons/2-2n)*6
    x = 3*n - round(numIons/2);
    y = round(numIons/2) - 2*n;
    for i = 1 : x   % 4 atoms per go - atoms on the main cube diagonals
      tmp = rand - 0.5;
     while abs(tmp) < 0.05  | (ellipse_mode & (norm([tmp tmp tmp]) > 0.5))
      tmp = rand - 0.5;
     end
      candidate = cat(1, candidate, [tmp tmp tmp]);
      candidate = cat(1, candidate, [tmp -1*tmp -1*tmp]);
      candidate = cat(1, candidate, [-1*tmp -1*tmp tmp]);
      candidate = cat(1, candidate, [-1*tmp tmp -1*tmp]);
    end
    for i = 1 : y  % 6 atoms per go - atoms on the XYZ axes
      tmp = rand - 0.5;
     while abs(tmp) < 0.05
      tmp = rand - 0.5;
     end
      candidate = cat(1, candidate, [0 0 tmp]);
      candidate = cat(1, candidate, [0 0 -1*tmp]);
      candidate = cat(1, candidate, [0 tmp 0]);
      candidate = cat(1, candidate, [0 -1*tmp 0]);
      candidate = cat(1, candidate, [tmp 0 0]);
      candidate = cat(1, candidate, [-1*tmp 0 0]);
    end
  end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isempty(findstr(nsym, 'Cv')) | ~isempty(findstr(nsym, 'cv')) % Cnv group - Cn axis and n vertical mirror planes
  rotRank = str2num(nsym(3:end));   % ceil(rand(1)*nsym);
  if rotRank > 2
   newLattice = [sqrt(lattice(1,1)*lattice(2,2)) 0 0; 0 sqrt(lattice(1,1)*lattice(2,2)) 0; 0 0 lattice(3,3)];
  end
  candidate = [];
  while 1
    if numIons == 1;   
      tmp = [0 0 0];          
      candidate = cat(1,candidate,tmp); 
      break;
    end
    if numIons < 2*rotRank 
     if numIons < rotRank      % put atoms on axis
      cand = [0 0 rand-0.5];
      for i = 2 : numIons
       tmp = [0 0 rand-0.5]; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
       cand = cat(1,cand,tmp); 
      end
      numIons = 0;
     else                   % put atoms in mirror planes: generate atom randomly then 'rotate' it till it meets the plane
      tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
      while ellipse_mode & (norm(tmp) > 0.5)
        tmp = rand(1,3) - 0.5;
      end

      r = norm(tmp);
      tmp = [sqrt(r^2-tmp(3)^2) 0 tmp(3)];  % mirror plan ZOX
      cand = tmp;
      for i = 1 : rotRank-1
       angle = 2*i*pi/rotRank;
       opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  % rotate around Z axis
       nextAtom = tmp*opemat;
       cand = cat(1,cand,nextAtom);
      end
      numIons = numIons - rotRank;
     end
    else
      tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
      while ellipse_mode & (norm(tmp) > 0.5)
        tmp = rand(1,3) - 0.5;
      end

      cand = tmp;
      if (tmp(1)*newLattice(1,1))^2 + (tmp(2)*newLattice(2,2))^2 < minDistance^2  % too close to the main axis
        numIons = numIons - 1;
        cand = [0 0 tmp(3)];          
      else
       tooClose = 0;
       for i = 1 : rotRank        % check if atom is too close to some mirror plane (in this case - put atom on a plane and only rotate)
        angle = 2*i*pi/rotRank;
        opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 1];  % reflection across a plane that goes through axis Z and makes an angle pi*i/n with the axis X
        nextAtom = tmp*opemat;
        if norm((nextAtom-tmp)*newLattice) < minDistance 
          tooClose = 1;
          tmp = (tmp + nextAtom)/2;
          cand = tmp;
          break;
        end
       end
     
       for i = 1 : rotRank-1
        angle = 2*i*pi/rotRank;
        opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  % rotate around Z axis
        nextAtom = tmp*opemat;
        cand = cat(1,cand,nextAtom);
       end
       for i = 1 : rotRank
        if tooClose
           break;
        end
        angle = 2*i*pi/rotRank;
        opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 1];  % reflection across a plane that goes through axis Z and makes an angle pi*i/n with the axis X
        nextAtom = tmp*opemat;
        cand = cat(1,cand,nextAtom);
       end 

       if tooClose
         numIons = numIons - rotRank;
       else
         numIons = numIons - 2*rotRank;
       end
      end
    end

    if isempty(candidate)
      candidate = cand; 
    else
      candidate = cat(1,candidate,cand); 
    end

    if numIons == 0 
      break
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isempty(findstr(nsym, 'Ch')) | ~isempty(findstr(nsym, 'ch')) % Cnh group - Cn axis and horisontal mirror plane
  rotRank = str2num(nsym(3:end));   % ceil(rand(1)*nsym);
  if rotRank > 2
   newLattice = [sqrt(lattice(1,1)*lattice(2,2)) 0 0; 0 sqrt(lattice(1,1)*lattice(2,2)) 0; 0 0 lattice(3,3)];
  end
  candidate = [];
  while 1
    if numIons == 1;   
      tmp = [0 0 0];          
      candidate = cat(1,candidate,tmp); 
      break;
    end
    if numIons < 2*rotRank 
     if numIons < rotRank      % put atoms on axis
      cand = [0 0 rand-0.5];
      while abs(cand(3)) < 0.01
        cand = [0 0 rand-0.5]; 
      end 
      cand = cat(1,cand,cand*Hz);    
      numIons = numIons - 2;
     else                   % put atoms in mirror plane XOY
      tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
      while ellipse_mode & (norm([tmp(1) tmp(2)]) > 0.5)
        tmp = rand(1,3) - 0.5;
      end
      tmp(3) = 0;
      cand = tmp;
      for i = 1 : rotRank-1
       angle = 2*i*pi/rotRank;
       opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  % rotate around Z axis
       nextAtom = tmp*opemat;
       cand = cat(1,cand,nextAtom);
      end
      numIons = numIons - rotRank;
     end
    else
     tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
     while ellipse_mode & (norm(tmp) > 0.5)
       tmp = rand(1,3) - 0.5;
     end
     if ((tmp(1)*newLattice(1,1))^2 + (tmp(2)*newLattice(2,2))^2 < minDistance^2) & (rotRank > 1)  % too close to the main axis
       numIons = numIons - 2;
       cand = [0 0 tmp(3); 0 0 -1*tmp(3)];             
     else      
      if abs(tmp(3)*newLattice(3,3)) < minDistance/2
        tooClose = 1;
        tmp(3) = 0;
      else
        tooClose = 0;
      end 
      cand = tmp;
      for i = 1 : rotRank-1
       angle = 2*i*pi/rotRank;
       opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  % rotate around Z axis
       nextAtom = tmp*opemat;
       cand = cat(1,cand,nextAtom);
      end
      if tooClose == 0
       cand1 = cand*Hz;
       cand = cat(1,cand,cand1);  % reflected atoms in the mirror plane
      end

      if tooClose      
       numIons = numIons - rotRank;
      else
       numIons = numIons - 2*rotRank;
      end
     end
    end

    if isempty(candidate)
      candidate = cand; 
    else
      candidate = cat(1,candidate,cand); 
    end

    if numIons == 0 
      break
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif  ~isempty(findstr(nsym, 'Dh')) | ~isempty(findstr(nsym, 'dh')) % dihedral symmetry with horisontal mirror plane XOY
  rotRank = str2num(nsym(3:end));   % ceil(rand(1)*nsym);
  if rotRank > 2
   newLattice = [sqrt(lattice(1,1)*lattice(2,2)) 0 0; 0 sqrt(lattice(1,1)*lattice(2,2)) 0; 0 0 lattice(3,3)];
  end
  candidate = [];
  while 1
    if numIons == 1;   
      tmp = [0 0 0];          
      candidate = cat(1,candidate,tmp); 
      break;
    end
    if numIons < 4*rotRank 
     if numIons > 2*rotRank      % put atoms in the mirror plane, another possibility is to put them above/below the rotation axis (so that rotation = mirroring)
                                 % (toDo: ALTERNATE OR CHOOSE METHOD RANDOMLY!)
      tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
      while ellipse_mode & (norm([tmp(1) tmp(2)]) > 0.5)
        tmp = rand(1,3) - 0.5;
      end
      tmp(3) = 0;
      cand = tmp;
      for i = 1 : rotRank-1
       angle = 2*i*pi/rotRank;
       opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  % rotate around Z axis
       nextAtom = tmp*opemat;
       cand = cat(1,cand,nextAtom);
      end
      for i = 1 : rotRank
       angle = 2*i*pi/rotRank;
       opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 -1];  % rotation by Pi around axis that lies in XOY and makes an angle pi*i/n with the axis X
       nextAtom = tmp*opemat;
       cand = cat(1,cand,nextAtom);
      end
      numIons = numIons - 2*rotRank;
     elseif numIons < rotRank      % put atoms on axis
      cand = [0 0 rand-0.5];
      cand = cat(1,cand,cand*Hz); % because of 2fold rotational axes and mirror plane
      numIons = numIons - 2;
     else                   % put atoms on 2nd order rotation axes
      tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
      r = norm(tmp);
      tmp = [r 0 0];  % axe X
      cand = tmp;
      for i = 1 : rotRank-1
       angle = 2*i*pi/rotRank;
       opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  % rotate around Z axis
       nextAtom = tmp*opemat;
       cand = cat(1,cand,nextAtom);
      end
      numIons = numIons - rotRank;
     end
    else
     tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
     while ellipse_mode & (norm(tmp) > 0.5)
       tmp = rand(1,3) - 0.5;
     end
     if ((tmp(1)*newLattice(1,1))^2 + (tmp(2)*newLattice(2,2))^2 < minDistance^2) & (rotRank > 1)  % too close to the main axis
       numIons = numIons - 2;
       cand = [0 0 tmp(3); 0 0 -1*tmp(3)];              
     else
      if abs(tmp(3)*newLattice(3,3)) < minDistance/2
        tooClose1 = 1;
        tmp(3) = 0;
      else
        tooClose1 = 0;
      end 
      tooClose2 = 0;
      for i = 1 : rotRank    % check if atom is too close to some rotation axis (in this case - put atom on this axis)
       angle = 2*i*pi/rotRank;
       opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 -1];  % rotation by Pi around axis that lies in XOY and makes an angle pi*i/n with the axis X
       nextAtom = tmp*opemat;
       if norm((nextAtom-tmp)*newLattice) < minDistance 
         tooClose2 = 1;
         tooClose1 = 1;      % since average atom is in the mirror plane (mirror plane contains rotation axes axes)
         tmp = (tmp + nextAtom)/2;
         cand = tmp;
         break;
       end
      end
     
      cand = tmp;
      for i = 1 : rotRank-1
       angle = 2*i*pi/rotRank;
       opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  % rotate around Z axis
       nextAtom = tmp*opemat;
       cand = cat(1,cand,nextAtom);
      end
      if tooClose2 == 0
       for i = 1 : rotRank
        angle = 2*i*pi/rotRank;
        opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 -1];  % rotation by Pi around axis that lies in XOY and makes an angle pi*i/n with the axis X
        nextAtom = tmp*opemat;
        cand = cat(1,cand,nextAtom);
       end
      end
      if tooClose1 == 0
       cand1 = cand*Hz;
       cand = cat(1,cand,cand1);  % reflected atoms in the mirror plane
      end

      if tooClose2      
       numIons = numIons - rotRank;
      elseif tooClose1
       numIons = numIons - 2*rotRank;
      else
       numIons = numIons - 4*rotRank;
      end
     end
    end

    if isempty(candidate)
      candidate = cand; 
    else
      candidate = cat(1,candidate,cand); 
    end

    if numIons == 0 
      break
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif  ~isempty(findstr(nsym, 'Dd')) | ~isempty(findstr(nsym, 'dd')) | ~isempty(findstr(nsym, 'Dv')) | ~isempty(findstr(nsym, 'dv')) % dihedral symmetry with vertical mirror planes (Dnd, Dnv)
  rotRank = str2num(nsym(3:end));   % ceil(rand(1)*nsym);
  if rotRank > 1
   newLattice = [sqrt(lattice(1,1)*lattice(2,2)) 0 0; 0 sqrt(lattice(1,1)*lattice(2,2)) 0; 0 0 lattice(3,3)];
  end
  candidate = [];
  while 1
    if numIons == 1;   
      tmp = [0 0 0];          
      candidate = cat(1,candidate,tmp); 
      break;
    end
    if numIons < 4*rotRank 
     if numIons < 2*rotRank      % put atoms on axis
      cand = [0 0 rand-0.5];
      cand = cat(1,cand,cand*Hz); % because of 2fold rotational axes
      numIons = numIons - 2;
     else                   % put atoms on 2nd order rotation axes, can also put on mirror planes (toDo: ALTERNATE OR CHOOSE METHOD RANDOMLY!)
      tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
      r = norm(tmp);
      tmp = [r 0 0];  % axe X
      cand = tmp;
      angleM = pi/rotRank;
      opematM = [cos(angleM) sin(angleM) 0; sin(angleM) -cos(angleM) 0; 0 0 1];  % mirror plane
      mirrorAtom = tmp*opematM;
      cand = cat(1,cand,mirrorAtom);
      for i = 1 : rotRank-1
       angle = 2*i*pi/rotRank;
       opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  % rotate around Z axis
       nextAtom = tmp*opemat;
       cand = cat(1,cand,nextAtom);
       nextMirrorAtom = mirrorAtom*opemat;
       cand = cat(1,cand,nextMirrorAtom);
      end
      numIons = numIons - 2*rotRank;
     end
    else
     tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
     while ellipse_mode & (norm(tmp) > 0.5)
       tmp = rand(1,3) - 0.5;
     end
     if (tmp(1)*newLattice(1,1))^2 + (tmp(2)*newLattice(2,2))^2 < minDistance^2  % too close to the main axis
       numIons = numIons - 1;
       cand = [0 0 tmp(3)];          
     else
      cand = tmp;
      tooClose1 = 0;
      for i = 1 : rotRank        % check if atom is too close to some mirror plane (in this case - put atom on that plane)
       angle = (2*i*pi/rotRank) + pi/rotRank;    
       opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 1];  % reflection across a plane that goes through axis Z and makes an angle pi*i/n with the axis X
       nextAtom = tmp*opemat;
       if norm((nextAtom-tmp)*newLattice) < minDistance 
         tooClose1 = 1;
         tmp = (tmp + nextAtom)/2;
         cand = tmp;
         break;
       end
      end
      tooClose2 = 0;
      for i = 1 : rotRank    % check if atom is too close to some rotation axis (in this case - put atom on this axis)
       angle = 2*i*pi/rotRank;
       opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 -1];  % rotation by Pi around axis that lies in XOY and makes an angle pi*i/n with the axis X
       nextAtom = tmp*opemat;
       if norm((nextAtom-tmp)*newLattice) < minDistance 
         tooClose2 = 1;
         tmp = (tmp + nextAtom)/2;
         cand = tmp;
         break;
       end
      end

      for i = 1 : rotRank-1
       angle = 2*i*pi/rotRank;
       opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  % rotate around Z axis
       nextAtom = tmp*opemat;
       cand = cat(1,cand,nextAtom);
      end
      if tooClose2 == 0
       for i = 1 : rotRank
        angle = 2*i*pi/rotRank;
        opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 -1];  % rotation by Pi around axis that lies in XOY and makes an angle pi*i/n with the axis X
        nextAtom = tmp*opemat;
        cand = cat(1,cand,nextAtom);
       end
      end

      if (tooClose2 == 1) | ((tooClose2 == 0) & (tooClose1 == 0))
       angleM = pi/rotRank;
       opematM = [cos(angleM) sin(angleM) 0; sin(angleM) -cos(angleM) 0; 0 0 1];  % mirror plane between two rotation axes (axis X being one of them)
       mirrorAtom = tmp*opematM;
       cand = cat(1,cand,mirrorAtom);
       for i = 1 : rotRank-1
        angle = 2*i*pi/rotRank;
        opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  % rotate around Z axis
        nextMirrorAtom = mirrorAtom*opemat;
        cand = cat(1,cand,nextMirrorAtom);
       end
       if tooClose2 == 0
        for i = 1 : rotRank
         angle = 2*i*pi/rotRank;
         opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 -1];  % rotation by Pi around axis that lies in XOY and makes an angle pi*i/n with the axis X
         nextMirrorAtom = mirrorAtom*opemat;
         cand = cat(1,cand,nextMirrorAtom);
        end
       end
      end

      if tooClose2 == 1       % atoms on axes
       numIons = numIons - 2*rotRank;
      elseif tooClose1 == 1   % atoms on mirror planes
       numIons = numIons - 2*rotRank;
      else
       numIons = numIons - 4*rotRank;
      end
     end
    end

    if isempty(candidate)
      candidate = cand; 
    else
      candidate = cat(1,candidate,cand); 
    end

    if numIons == 0 
      break
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif  ~isempty(findstr(nsym, 'D')) | ~isempty(findstr(nsym, 'd')) % dihedral symmetry
  rotRank = str2num(nsym(2:end));   % ceil(rand(1)*nsym);
  if rotRank > 2
   newLattice = [sqrt(lattice(1,1)*lattice(2,2)) 0 0; 0 sqrt(lattice(1,1)*lattice(2,2)) 0; 0 0 lattice(3,3)];
  end
  candidate = [];
  while 1
    if numIons == 1;   
      tmp = [0 0 0];          
      candidate = cat(1,candidate,tmp); 
      break;
    end
    if numIons < 2*rotRank 
     if numIons < rotRank      % put atoms on axis
      cand = [0 0 rand-0.5];
      cand = cat(1,cand,cand*Hz); % because of 2fold rotational axes 
      numIons = numIons - 2;
     else                   % put atoms on 2nd order rotation axes
      tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
      r = norm(tmp);
      tmp = [r 0 0];  % axe X
      cand = tmp;
      for i = 1 : rotRank-1
       angle = 2*i*pi/rotRank;
       opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  % rotate around Z axis
       nextAtom = tmp*opemat;
       cand = cat(1,cand,nextAtom);
      end
      numIons = numIons - rotRank;
     end
    else
     tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
     while ellipse_mode & (norm(tmp) > 0.5)
       tmp = rand(1,3) - 0.5;
     end
     if (tmp(1)*newLattice(1,1))^2 + (tmp(2)*newLattice(2,2))^2 < minDistance^2  % too close to the main axis
       numIons = numIons - 1;
       cand = [0 0 tmp(3)];          
     else
      cand = tmp;
      tooClose = 0;
      for i = 1 : rotRank    % check if atom is too close to some rotation axis (in this case - put atom on this axis)
       angle = 2*i*pi/rotRank;
       opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 -1];  % rotation by Pi around axis that lies in XOY and makes an angle pi*i/n with the axis X
       nextAtom = tmp*opemat;
       if norm((nextAtom-tmp)*newLattice) < minDistance 
         tooClose = 1;
         tmp = (tmp + nextAtom)/2;
         cand = tmp;
         break;
       end
      end

      for i = 1 : rotRank-1
       angle = 2*i*pi/rotRank;
       opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  % rotate around Z axis
       nextAtom = tmp*opemat;
       cand = cat(1,cand,nextAtom);
      end
      if tooClose == 0
       for i = 1 : rotRank
        angle = 2*i*pi/rotRank;
        opemat = [cos(angle) sin(angle) 0; sin(angle) -cos(angle) 0; 0 0 -1];  % rotation by Pi around axis that lies in XOY and makes an angle pi*i/n with the axis X
        nextAtom = tmp*opemat;
        cand = cat(1,cand,nextAtom);
       end
      end

      if tooClose
       numIons = numIons - rotRank;
      else
       numIons = numIons - 2*rotRank;
      end
     end
    end

    if isempty(candidate)
      candidate = cand; 
    else
      candidate = cat(1,candidate,cand); 
    end

    if numIons == 0 
      break
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~isempty(findstr(nsym, 'O')) | ~isempty(findstr(nsym, 'o'))   % O (cube rotational) and Oh (full cube) symmetries
  newLattice = [volLat^(1/3) 0 0; 0 volLat^(1/3) 0; 0 0 volLat^(1/3)];
  candidate = [];
  if (numIons ~= 1) & (numIons ~= 6) & (numIons ~= 7) & (numIons ~= 8) & (numIons ~= 9) & (numIons < 12)
    status = ['Impossible to build the cluster with ' num2str(numIons) ' atoms that has symmetry group ' nsym];
    [nothing, nothing] = unix(['echo ' status ' > error_cluster_symmetry']);
    errorS = 1;
  end
  if mod(numIons,2) == 1 % put in the center
      candidate = [0 0 0];
      numIons = numIons - 1;
  end
  n24 = 0;
  numIons1 = numIons;
  while 1
    numIons1 = numIons1 - 24;
    if (numIons1 ~= 6) & (numIons1 ~= 8) & (numIons1 ~= 0) & (numIons1 < 12)
      break;
    else
      n24 = n24 + 1;
    end
  end
  n = 1;
  while n <= n24  % do a 'full' symmetry point
     tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
     while ellipse_mode & (norm(tmp) > 0.5)
       tmp = rand(1,3) - 0.5;
     end
     if (strcmp(nsym, 'Oh') | strcmp(nsym, 'oh')) & (mod(n24,2) == 1) & (n == n24)  % Oh (full cube) symmetry
       w = ceil(rand*2);
       if w == 1 % do only 24 poins :) (by generating a point on mirror plane containing some axis)
         tmp(ceil(rand*3)) = 0;  
       else % do only 24 poins, different way :)  (by generating a point on mirror plane containing some face diagonal)
        w = ceil(rand*3);
        if w == 1
         tmp(1) = tmp(2);  
        elseif w == 2
         tmp(1) = tmp(3);  
        else
         tmp(3) = tmp(2);  
        end
       end
     end
     cand = tmp;
   % rotations around X, Y, Z (O is in the cube center, Z goes up)
     for i = 1 : 3
      angle = i*pi/2;
      opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  % rotate around Z axis
      nextAtom = tmp*opemat;
      cand = cat(1,cand,nextAtom);
     end
     for i = 1 : 3
      angle = i*pi/2;
      opemat = [1 0 0; 0 cos(angle) -sin(angle); 0 sin(angle) cos(angle)];  % rotate around X axis
      nextAtom = tmp*opemat;
      cand = cat(1,cand,nextAtom);
     end
     for i = 1 : 3
      angle = i*pi/2;
      opemat = [cos(angle) 0 sin(angle); 0 1 0; -sin(angle) 0 cos(angle)];  % rotate around Y axis     
      nextAtom = tmp*opemat;
      cand = cat(1,cand,nextAtom);
     end
   % rotations around 2-fold axes in XOY plane (connecting the middle of cube edges) 
     nextAtom = tmp*[0 1 0; 1 0 0; 0 0 -1];
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 -1 0; -1 0 0; 0 0 -1];
     cand = cat(1,cand,nextAtom);
   % rotations around 2-fold axes in YOZ plane (connecting the middle of cube edges) 
     nextAtom = tmp*[-1 0 0; 0 0 1; 0 1 0];
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[-1 0 0; 0 0 -1; 0 -1 0];
     cand = cat(1,cand,nextAtom);
   % rotations around 2-fold axes in ZOX plane (connecting the middle of cube edges) 
     nextAtom = tmp*[0 0 1; 0 -1 0; 1 0 0];
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 0 -1; 0 -1 0; -1 0 0];
     cand = cat(1,cand,nextAtom);    
   % rotations around the main cube diagonals (3-fold axes) 
     nextAtom = tmp*[0 1 0; 0 0 1; 1 0 0]; % +++ diagonal, -120 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 0 1; 1 0 0; 0 1 0]; % -240 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 1 0; 0 0 -1; -1 0 0]; % ++- diagonal, 120 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 0 -1; 1 0 0; 0 -1 0]; % 240 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 -1 0; 0 0 1; -1 0 0]; % +-- diagonal, 120 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 0 -1; -1 0 0; 0 1 0]; % 240 degree     
     cand = cat(1,cand,nextAtom);     
     nextAtom = tmp*[0 -1 0; 0 0 -1; 1 0 0]; % +-+ diagonal, 120 degree
     cand = cat(1,cand,nextAtom);
     nextAtom = tmp*[0 0 1; -1 0 0; 0 -1 0]; % 240 degree     
     cand = cat(1,cand,nextAtom); 
     if (strcmp(nsym, 'Oh') | strcmp(nsym, 'oh')) & ~((mod(n24,2) == 1) & (n == n24))  % Oh (full cube) symmetry
       numIons = numIons - 48;
       n = n + 2;
       cand1 = -1*cand;
       cand = cat(1,cand,cand1); % added inversion
     else
       numIons = numIons - 24;
       n = n + 1;
     end
     if isempty(candidate)
      candidate = cand; 
     else
      candidate = cat(1,candidate,cand); 
     end
     if n > n24 
         break
     end
  end
  if numIons > 0  % numIons can be described as 6x+8y
    n1 = ceil(numIons/8);
    n2 = floor(numIons/6);
    n = n1 + round(rand*(n2-n1)); % numIons = (4n-numIons/2)*6 + (numIons/2-3n)*8
    x = 4*n - round(numIons/2);
    y = round(numIons/2) - 3*n;
    for i = 1 : x   % 6 atoms per go
      tmp = rand - 0.5;
     while abs(tmp) < 0.05
      tmp = rand - 0.5;
     end
      candidate = cat(1, candidate, [0 0 tmp]);
      candidate = cat(1, candidate, [0 0 -1*tmp]);
      candidate = cat(1, candidate, [0 tmp 0]);
      candidate = cat(1, candidate, [0 -1*tmp 0]);
      candidate = cat(1, candidate, [tmp 0 0]);
      candidate = cat(1, candidate, [-1*tmp 0 0]);
    end
    for i = 1 : y  % 8 atoms per go
      tmp = rand - 0.5;
     while abs(tmp) < 0.05 | ((ellipse_mode == 1) & (norm([tmp tmp tmp]) > 0.5))
      tmp = rand - 0.5;
     end
      candidate = cat(1, candidate, [tmp tmp tmp]);
      candidate = cat(1, candidate, [tmp tmp -1*tmp]);
      candidate = cat(1, candidate, [tmp -1*tmp tmp]);
      candidate = cat(1, candidate, [tmp -1*tmp -1*tmp]);
      candidate = cat(1, candidate, [-1*tmp tmp tmp]);
      candidate = cat(1, candidate, [-1*tmp tmp -1*tmp]);
      candidate = cat(1, candidate, [-1*tmp -1*tmp tmp]);
      candidate = cat(1, candidate, [-1*tmp -1*tmp -1*tmp]);
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else    % Cn symmetry
  rotRank = str2num(nsym(2:end));   % ceil(rand(1)*nsym);
  if rotRank > 2
   newLattice = [sqrt(lattice(1,1)*lattice(2,2)) 0 0; 0 sqrt(lattice(1,1)*lattice(2,2)) 0; 0 0 lattice(3,3)];
  end
  candidate = [];
  while 1
    if numIons < rotRank    % put the rest of the atoms on the axis
      cand = [0 0 rand-0.5];
      for i = 2 : numIons
       tmp = [0 0 rand-0.5]; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
       cand = cat(1,cand,tmp); 
      end
      numIons = 0;
    else
     tmp = rand(1,3) - 0.5; % we work in the space [-0,5:0.5;-0,5:0.5;-0,5:0.5] and then add (0.5,0.5,0.5)
     while ellipse_mode & (norm(tmp) > 0.5)
       tmp = rand(1,3) - 0.5;
     end
     if (tmp(1)*newLattice(1,1))^2 + (tmp(2)*newLattice(2,2))^2 < minDistance^2  % too close to the main axis
       numIons = numIons - 1;
       cand = [0 0 tmp(3)];          
     else
      cand = tmp;
      for i = 1 : rotRank-1
       angle = 2*i*pi/rotRank;
       opemat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];  % rotate around Z axis
       nextAtom = tmp*opemat;
       cand = cat(1,cand,nextAtom);
      end
      numIons = numIons - rotRank;      
     end
    end

    if isempty(candidate)
      candidate = cand; 
    else
      candidate = cat(1,candidate,cand); 
    end

    if numIons == 0 
      break
    end
  end
end


AbsoluteCoord = candidate*newLattice;

%unix(['echo cluster > ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo 1.0 >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo ' num2str(newLattice(1,:)) ' >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo ' num2str(newLattice(2,:)) ' >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo ' num2str(newLattice(3,:)) ' >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo ' num2str(nI) ' >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo Direct >> ClusterPOSCAR' num2str(goodPop)]);
%for i = 1 : nI
%  unix(['echo ' num2str(AbsoluteCoord(i,:)) ' >> ClusterPOSCAR' num2str(goodPop)]);
%end

[a,b] = PrincipleAxis(AbsoluteCoord);  %find the principle roation axis

if strcmp(nsym, 'E')
  AbsoluteCoord = AbsoluteCoord*a; % maximal moment of inertia - for Z axis, minimum - for X
  % we will fix the positions so that no atom is outside of the unit cell and cluster basically fits the whole cell
%  mass_center = zeros(1,3);
%  for i = 1 : 3 
%   mass_center(i) = sum(AbsoluteCoord(:,i));
%   AbsoluteCoord(:,i) = AbsoluteCoord(:,i) - mass_center(i); % 'center' the cluster, center of mass should be exactly at [0.5;0.5;0.5]
%  end
  ma = max(AbsoluteCoord);
  mi = min(AbsoluteCoord);
  newLattice = zeros(3,3);
  newLattice(1,1) = ma(1) - mi(1) + 0.02;
  newLattice(2,2) = ma(2) - mi(2) + 0.02;
  newLattice(3,3) = ma(3) - mi(3) + 0.02;
  AbsoluteCoord(:,1) = AbsoluteCoord(:,1) - mi(1) + 0.01; % 'center' the cluster
  AbsoluteCoord(:,2) = AbsoluteCoord(:,2) - mi(2) + 0.01;
  AbsoluteCoord(:,3) = AbsoluteCoord(:,3) - mi(3) + 0.01;
  candidate = (AbsoluteCoord/newLattice) - 0.5;
elseif ~isempty(findstr(nsym, 'O')) | ~isempty(findstr(nsym, 'o')) | ~isempty(findstr(nsym, 'T')) | ~isempty(findstr(nsym, 't')) % cubic lattice
  candidate = candidate*a; % maximal moment of inertia - for Z axis, minimum - for X
  % we will fix the positions so that no atom is outside of the unit cell and cluster basically fits the whole cell
  m = max(max(abs(candidate)));
  candidate = candidate*(0.5/(m + 0.01));
else % 'cylindrical' lattice and inversion
% candidate = candidate*a;
% lat = inv(a)*lat*a; % this way AbsoluteCoord => AbsoluteCoord*a (we rotate lattice and relative coordinates)
  AbsoluteCoord = AbsoluteCoord*a; % maximal moment of inertia - for Z axis, minimum - for X
  mass_center = zeros(1,3);
  for i = 1 : 3 
   mass_center(i) = sum(AbsoluteCoord(:,i))/nI;
   AbsoluteCoord(:,i) = AbsoluteCoord(:,i) - mass_center(i); % 'center' the cluster, center of mass should be exactly at [0.5;0.5;0.5]
  end

  % we will fix the positions so that no atom is outside of the unit cell and cluster basically fits the whole cell
  m = max(abs(AbsoluteCoord));
  newLattice = zeros(3,3);
  newLattice(1,1) = 2 * m(1) + 0.02;
  newLattice(2,2) = 2 * m(2) + 0.02;
  newLattice(3,3) = 2 * m(3) + 0.02;
  candidate = AbsoluteCoord/newLattice;
end


candidate = candidate + 0.5;  

%unix(['echo cluster >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo 1.0 >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo ' num2str(newLattice(1,:)) ' >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo ' num2str(newLattice(2,:)) ' >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo ' num2str(newLattice(3,:)) ' >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo ' num2str(nI) ' >> ClusterPOSCAR' num2str(goodPop)]);
%unix(['echo Direct >> ClusterPOSCAR' num2str(goodPop)]);
%for i = 1 : nI
%  unix(['echo ' num2str(candidate(i,:)) ' >> ClusterPOSCAR' num2str(goodPop)]);
%end

%nAxis = RandInt(1,1,[1,3])
%if nAxis == 1
% transform = [0 0 -1; 0 1 0; 1 0 0];   % to align the rotation axis with the X
%elseif nAxis == 2
% transform = [1 0 0; 0 0 -1; 0 1 0];   % to align the rotation axis with the Y
%else
% transform = [1 0 0; 0 1 0; 0 0 1];    % rotation axis stays as Z
%end

%status = candidate*transform;

