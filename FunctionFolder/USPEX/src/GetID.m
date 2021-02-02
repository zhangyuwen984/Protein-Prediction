function [i0,j0]=GetID(atom, MtypeLIST, typesAList)
global ORG_STRUC

item = 0;
i0 = 0;
j0 = 0;
for i=1:length(MtypeLIST)
   for j=1:length(ORG_STRUC.STDMOL(MtypeLIST(i)).types)
       item = item+1;
       for ind=1:length(ORG_STRUC.atomType)
           if typesAList(item)==ORG_STRUC.atomType(ind)
              if item==atom
                 i0 = i;
                 j0 = j;          
                 break
              end
           end
       end
   end
end
