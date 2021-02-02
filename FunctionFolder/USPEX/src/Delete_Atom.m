function [coor, type] = Delete_Atom(coor, type, NDelete, atomType)

%Delete atoms according to the height in c-axis

ID = zeros(1,length(type));

count = 0;
[tmp, rank] = sort(coor(:,3));
for i = 0 : length(type)-1
    ID_tmp = rank(end-i);
    for j = 1:length(atomType)
       if (type(ID_tmp) == atomType(j)) & (NDelete(j) > 0)
          count = count + 1;
          ID(ID_tmp) = 1;
       end
    end
    if count == sum(NDelete)
       break;
    end
end
type(find(ID==1))   = [];
coor(find(ID==1),:) = [];

