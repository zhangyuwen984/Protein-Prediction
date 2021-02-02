function [newlattice,newcoord,newtypes,neworders,numIons]=cresupercell(lattice, coordinates, atomtype, order , cell1, cell2)
global ORG_STRUC      
%%This program is to covert the cell1 to the cell2
%%The basic idea is like this, let take the example of coversion from 2*2 to 3*3
%%firstly creat supercell 6*6 based on 2*2, and then cut the supercell as 4
%%regions, finally take one region randomly.

newlattice=lattice;
newcoord=[];
neworders=[];
newtypes=[];
numIons=[];
if size(coordinates,1)==0   
    return
else 
   %%%%%1:create supercell
   super(1)=lcm(cell1(1),cell2(1));
   super(2)=lcm(cell1(2),cell2(2));
   num1=super(1)/cell1(1);
   num2=super(2)/cell1(2);
   
   [S_coord, S_lat, N_Size] = Make_SuperCell(coordinates, lattice, [num1, num2], 1);
   S_coord  = S_coord-floor(S_coord);
   S_Natoms = N_Size*size(coordinates, 1);
   tmp = repmat(atomtype, [N_Size,1]);
   S_type = reshape(tmp, [1, S_Natoms]);
   tmp = repmat(order, [N_Size,1]);
   S_order = reshape(tmp, [1, S_Natoms]);

   %%%2: chose the region
   num3=super(1)/cell2(1);
   num4=super(2)/cell2(2);
   %%%%%%%%%%%%%%how many atoms in each region
   max_atom = 0;
   for i=1:num3
       for j=1:num4
           ioncount=0;
           X_max =     i/num3 - 0.0001;
           Y_max =     j/num4 - 0.0001;
           X_min = (i-1)/num3 - 0.0001;
           Y_min = (j-1)/num4 - 0.0001;
           for k=1:S_Natoms
               if     ( X_min < S_coord(k,1) ) && ( S_coord(k,1) < X_max ) ...
                   && (  Y_min< S_coord(k,2) ) && ( S_coord(k,2) < Y_max )
                   ioncount=ioncount+1;
               end
           end
           if ioncount>max_atom
              max_atom = ioncount;
              pick = [i,j];
           end      
       end
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   newlattice(1,:)=S_lat(1,:)/num3;
   newlattice(2,:)=S_lat(2,:)/num4;
   newlattice(3,:)=S_lat(3,:);
   X_max =     pick(1)/num3 - 0.0001;
   Y_max =     pick(2)/num4 - 0.0001;
   X_min = (pick(1)-1)/num3 - 0.0001;
   Y_min = (pick(2)-1)/num4 - 0.0001;
   ioncount = 0;

   for i=1:S_Natoms
       if     ( X_min < S_coord(i,1) ) && ( S_coord(i,1) < X_max ) ...
          && (  Y_min < S_coord(i,2) ) && ( S_coord(i,2) < Y_max )
           ioncount=ioncount+1;
           newcoord(ioncount,:)= S_coord(i,:)*S_lat/newlattice;
           newtypes(ioncount)  = S_type(i);
           neworders(ioncount) = S_order(i);
       end
   end

   for i=1:size(ORG_STRUC.atomType,2)
       numIons(i)=sum(size(find(newtypes==ORG_STRUC.atomType(i)),2));
   end
end

newcoord  = newcoord-floor(newcoord);
