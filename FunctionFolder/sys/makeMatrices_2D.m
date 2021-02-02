function [N, V, dist_matrix, typ_i, typ_j, ho, ht] = makeMatrices_2D(lattice, coordinates, numIons, atomType)
%INPUT:
% lattice, coordinates of a given crystal structure;
% numIons: 4 4
% atomType: [12 8] (Mg O)

%OUTPUT
% N(i) number of atoms in unit cell of the type of i-th atom
% for example if ntyp = 1 1 4, then N = 1 1 4 4 4 4
% if atomType = 10 20 30, then atomType1 = 10 20 30 30 30 30
% dist_matrix;
% typ_i = 1 2 3 3 3 3 (needed for addressing in multicomponent fingerprint)
% typ_i,j replace addr(i,j) = (typ_i-1)*species + typ_j
% ho: The distance from the atom to the lower bound
% ht: The distance from the atom to the upper bound
global ORG_STRUC 

Rmax  = ORG_STRUC.RmaxFing;

species = size(numIons,2);
natom = sum(numIons);

N = zeros(natom,1);
typ_i = zeros(natom,1);
typ_j = 0;
Nfull = 0;
dist_matrix = 0;
V = abs(det(lattice)); 

if sum(numIons)>0
      c = 0;
      for i = 1:species
          for j = 1:numIons(i)
              c = c+1;
              N(c) = numIons(1,i);
              atomType1(c) = atomType(1,i);
              typ_i(c) = i;
          end
      end
      
      
      signum = [1 1 1; -1 1 1; 1 -1 1; 1 1 -1; -1 -1 1; 1 -1 -1; -1 1 -1; -1 -1 -1];
      condition = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 1 1];
      % condition shows what variables shouldn't be equal to zero
      
      vect = zeros(13,3);
      abs_vect = zeros(13,1);
      vect(1,:) = lattice(1,:);
      vect(2,:) = lattice(2,:);
      vect(3,:) = lattice(3,:);
      vect(4,:) = vect(1,:)+vect(2,:);
      vect(5,:) = vect(1,:)-vect(2,:);
      vect(6,:) = vect(1,:)+vect(3,:);
      vect(7,:) = vect(1,:)-vect(3,:);
      vect(8,:) = vect(3,:)+vect(2,:);
      vect(9,:) = vect(3,:)-vect(2,:);
      vect(10,:) = vect(1,:)+vect(2,:)+vect(3,:);
      vect(11,:) = vect(1,:)+vect(2,:)-vect(3,:);
      vect(12,:) = vect(1,:)-vect(2,:)+vect(3,:);
      vect(13,:) = -vect(1,:)+vect(2,:)+vect(3,:);
      
      for i=1:13
        abs_vect(i)=sqrt(dot(vect(i,:),vect(i,:)));  
      end
      lengthX = ceil((Rmax+max(abs_vect))/min(abs_vect))+1; % cells that we copy along x : -lenghtX to lenghtX
      lengthY = lengthX;
      lengthZ = lengthX;
      
      for i = 0:lengthX
        quit_marker_x = 1;     
        for j = 0:lengthY
          quit_marker_y = 1;      
          for k = 0:lengthZ
             if (k > 0)
              continue;
             end
             quit_marker_z = 1;
             for quadrants = 1:8
               if (condition(quadrants,1)*(i == 0)+condition(quadrants,2)*(j == 0)+condition(quadrants,3)*(k == 0)) == 0
                for current_cell = 1:natom
                 distances = zeros(1,natom);
                 ho = zeros(1,natom);
                 ht = zeros(1,natom);
                 XYZcoordinate = coordinates*lattice;
                 Zcoordinate = XYZcoordinate(:,3)' ;
                 marker = 0; % shows if atom current_cell is in range with any of the atoms in unit cell
                 for basic_cell = 1:natom
                     x = coordinates(current_cell,1)+signum(quadrants,1)*i-coordinates(basic_cell,1); % x in lat basis
                     y = coordinates(current_cell,2)+signum(quadrants,2)*j-coordinates(basic_cell,2); % y in lat basis
                     z = coordinates(current_cell,3)+signum(quadrants,3)*k-coordinates(basic_cell,3); % z in lat basis
                     Rij = (x*lattice(1,1)+y*lattice(2,1)+z*lattice(3,1))^2;
                     Rij = Rij + (x*lattice(1,2)+y*lattice(2,2)+z*lattice(3,2))^2;
                     Rij = Rij + (x*lattice(1,3)+y*lattice(2,3)+z*lattice(3,3))^2;
                       ho(basic_cell) = Zcoordinate(basic_cell);
                       ht(basic_cell) = lattice(3,3) - ho(1,basic_cell);                 
                       if Rij < Rmax^2
                       quit_marker_z = 0;
                       quit_marker_y = 0;
                       quit_marker_x = 0;
                       if (marker == 0) && (Nfull >= natom)
                           N = vertcat(N,N(current_cell));
                       end
                       Nfull=Nfull+1-marker;
                       marker = 1;
                       typ_j_coef = typ_i(current_cell);
                       distances(basic_cell) = sqrt(Rij);
                   end
                 end
                  if (i+j+k+current_cell == 1)
                      dist_matrix = distances;
                      typ_j = typ_j_coef;
                  elseif marker
                      dist_matrix = vertcat(dist_matrix, distances);
                      typ_j = vertcat(typ_j, typ_j_coef);
                  end
                end
               end
             end
             if quit_marker_z
                 break;
             end
          end   
          if quit_marker_y
              break;
          end    
        end
        if quit_marker_x
            break;
        end      
      end
end
