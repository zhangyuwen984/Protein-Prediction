function format = findformat(Pair)
%Zmatrix needs the i-j-k table to define the bond length/angle/tosional.
%This function is to find the i-j-k table from the resorted neighbor list
%INPUT: Pair(N*3)
%OUTPUT: format (N*3) -- i,j,k table


N_max = 12;
n_atom = size(Pair, 1);
if n_atom == 1
   format(1,:) = [0 0 0];
elseif n_atom == 2
   format(1,:) = [0 0 0];
   format(2,:) = [1 0 0];
else
   format(2,:) = [1 0 0];
   if Pair(1,end) > 1
      format(3,:) = [1 2 0];
   else
      format(3,:) = [2 1 0];
   end
   if n_atom > 3
      for i=4:n_atom
          format(i,1)=Pair(i,1);
          format(i,2)=searchID(i, format(i,1), Pair);
          format(i,3)=searchID(i, [format(i,2), format(i,1)], Pair);
      end
   end
end

%for the i-j-k table, the rule is
%1, the i-j-k must be connected
%2, 0<format(i,:)<i
%3, format(i,1)~=format(i,2)~=format(i,3)
function ID = searchID(Max_num, ID_exist, ID_Data)
for i = 1:length(ID_exist)
    for j = 1:ID_Data(ID_exist(i),end)
        f = ID_Data(ID_exist(i),j);
        if (f>0) & (f<Max_num) & (sum(f== ID_exist)==0)
            ID = f;
            break;
        end
    end
end
