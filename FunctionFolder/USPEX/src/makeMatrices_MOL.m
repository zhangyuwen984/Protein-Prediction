function [N, V, dist_matrix, typ_i, typ_j] = makeMatrices_MOL(lattice_vec, atom_positions, ntyp, atom_num)

species = size(ntyp,2);
N = zeros(natom,1);
atom_num1 = zeros(natom,1);
typ_i = zeros(natom,1);
Nfull = 0;
Mol=Judge_Mol_Index();
% N(i) =  number of atoms in unit cell of the type of i-th atom
% for example if ntyp = 1 1 4, then N = 1 1 4 4 4 4
% if atom_num = 10 20 30, then atom_num1 = 10 20 30 30 30 30
% typ_i = 1 2 3 3 3 3 (needed for addressing in multicomponent fingerprint)
% typ_i,j replace addr(i,j) = (typ_i-1)*species + typ_j
% to add - scale_factor is not used yet
A = 0;
for i = 1:species
    for j = 1:species
        A = A+ntyp(i)*ntyp(j)*atom_num(i)*atom_num(j)/natom;
    end
end
c = 0;
for i = 1:species
    for j = 1:ntyp(i)
        c = c+1;
        N(c) = ntyp(1,i);
        atom_num1(c) = atom_num(1,i);
        typ_i(c) = i;
    end
end

V = abs(det(lattice_vec)); 
signum = [1 1 1; -1 1 1; 1 -1 1; 1 1 -1; -1 -1 1; 1 -1 -1; -1 1 -1; -1 -1 -1];
condition = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 1 1];
% condition shows what variables shouldn't be equal to zero
vect = zeros(13,3);
abs_vect = zeros(13,1);
vect(1,:) = lattice_vec(1,:);
vect(2,:) = lattice_vec(2,:);
vect(3,:) = lattice_vec(3,:);
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
       quit_marker_z = 1;
       for quadrants = 1:8
         if (condition(quadrants,1)*(i == 0)+condition(quadrants,2)*(j == 0)+condition(quadrants,3)*(k == 0)) == 0
          for current_cell = 1:natom
           distances = zeros(1,natom);
           Wcoef = zeros(1,natom);
           addr_coef = zeros(1,natom);
           marker = 0; % shows if atom current_cell is in range with any of the atoms in unit cell
           for basic_cell = 1:natom
	       if (~(Mol(basic_cell,4)==Mol(current_cell,4))) || (i+j+k >0)
                  x=atom_positions(current_cell,1)+signum(quadrants,1)*i-atom_positions(basic_cell,1); 
                  y=atom_positions(current_cell,2)+signum(quadrants,2)*j-atom_positions(basic_cell,2);
                  z=atom_positions(current_cell,3)+signum(quadrants,3)*k-atom_positions(basic_cell,3);
                  Rij = (x*lattice_vec(1,1)+y*lattice_vec(2,1)+z*lattice_vec(3,1))^2;
                  Rij = Rij + (x*lattice_vec(1,2)+y*lattice_vec(2,2)+z*lattice_vec(3,2))^2;
                  Rij = Rij + (x*lattice_vec(1,3)+y*lattice_vec(2,3)+z*lattice_vec(3,3))^2;
                  if Rij < Rmax^2
                     quit_marker_z = 0;
                     quit_marker_y = 0;
                     quit_marker_x = 0;
                     if (marker == 0) & (Nfull >= natom)
                         N = vertcat(N,N(current_cell));
                     end
                     Nfull=Nfull+1-marker;
                     marker = 1;
                     Wcoef(basic_cell) = atom_num1(basic_cell)*atom_num1(current_cell);
                     typ_j_coef = typ_i(current_cell);
                     distances(basic_cell) = sqrt(Rij);
                  end
	       end
           end
            if (i+j+k+current_cell == 1)
                W = Wcoef;
                dist_matrix = distances;
                typ_j = typ_j_coef;
            elseif marker
                W = vertcat(W,Wcoef);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mol1=Judge_Mol_Index()

global ORG_STRUC

N=sum(ORG_STRUC.numIons);

Mol=zeros(N,5);
Mol1=zeros(N,5);
count=zeros(ORG_STRUC.atomType);
item=1;
while(item < N)
   for i=1:length(ORG_STRUC.MtypeLIST)
       numatom=length(ORG_STRUC.STDMOL(ORG_STRUC.MtypeLIST(i)).types);
      for j=1:numatom
	for index=1:length(ORG_STRUC.atomType)
            type=ORG_STRUC.atomType(index);
	    if(ORG_STRUC.atomType(ORG_STRUC.STDMOL(ORG_STRUC.MtypeLIST(i)).types(j))==type)
	       count(index)=count(index)+1;
	       Mol(item,:)=[item, index, count(index), i, j];
           item=item+1;
            end
        end
      end
    end
end
item=1
    for index=1:length(ORG_STRUC.atomType)
        for j=1:count(index)
	    for i=1:N
                if(Mol(i,2)==index & (Mol(i,3)==j))
                   Mol1(item,:)=Mol(i,:);
                end
            end
            item=item+1;
	end
    end

