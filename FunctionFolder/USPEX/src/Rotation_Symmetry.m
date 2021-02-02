function [Son,rotNow] = Rotation_Symmetry(father, lat, MtypeLIST)

global ORG_STRUC

warning off all
STDMOL = ORG_STRUC.STDMOL;
numMols= length(father);
lat_6  = latConverter(lat);
rotNow = numMols;

N_Wycc = 0;
Sym = 0;
for i = 1: sum(numMols)
    if ~isempty(father(i).Operation)
       if isequal(father(i).Operation, [1 0 0; 0 1 0; 0 0 1; 0 0 0])
          N_Wycc = N_Wycc + 1;
          Sym(N_Wycc) = 1;
       else
          Sym(N_Wycc) = Sym(N_Wycc) + 1; %[4 4]
       end
    end
end

Son = father;
count = 0;
for i = 1: N_Wycc
   Coords = father(count+1).MOLCOORS;
   format = STDMOL(MtypeLIST(count+1)).format;
   num_optFlags = STDMOL(MtypeLIST(count+1)).num_optFlags;
   [Zmatrix, REF_MOLCOORS] = RotInertia(Coords,format, num_optFlags, MtypeLIST(count+1));

   %%Obtain the rest rotations
   for j = 1:Sym(i)
       R   = Operate_molecule(REF_MOLCOORS, Son(count+j).Operation, lat_6, [1 2 3], [1 2 3]);
       Son(count+j).MOLCOORS= R;
       Son(count+j).ZMATRIX = NEW_coord2Zmatrix(R,format);
   end
   count = count + Sym(i);
end
