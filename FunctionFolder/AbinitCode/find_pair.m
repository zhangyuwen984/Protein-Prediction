function Pair = find_pair(coor, radii)
%This function is all the atom pairs and constructe the neighbor list
%The bond lenght is estimated by the covalent radii.
%INPUT: coor (n*3 matrix), radii (n*1 matrix)
%OUTPUT: Pair (n*7 matrix)
N_max = 6;
n_atom = length(radii);
%Maximum 6 coordination, 7 gives the coordination number
Pair = zeros(n_atom, N_max); 
for i = 1:n_atom-1
    for j= i+1:n_atom
        dist = norm(coor(i,:)-coor(j,:));
        if dist < 1.2*(radii(i)+radii(j))
           Pair(i,N_max) = Pair(i,N_max) + 1;
           Pair(i, Pair(i,N_max)) = j;
           Pair(j,N_max) = Pair(j,N_max) + 1;
           Pair(j, Pair(j,N_max)) = i;
        end
    end
end

%We assume there is no isolated atom
if n_atom > 1
   for i = 1:n_atom
      if (Pair(i, N_max) ==0) & (N_max > 1)
         disp(['atom_' num2str(i) ' is not connected to any other atom']);
         disp(['Please check your MOL file again. Serious WARNING.... ']);
      end
   end
end
