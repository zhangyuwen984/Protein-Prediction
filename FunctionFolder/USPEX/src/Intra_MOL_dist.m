function Intra_map = Intra_MOL_dist(MtypeLIST, numIons, STDMOL)

N_type = length(numIons);
N_Mtype = length(MtypeLIST);

N_atom = sum(numIons);
ID = zeros(N_atom,1);
Intra_map = ones(N_atom, N_atom);

count = zeros(N_type, 1);
if N_type > 1
   for i = 2: N_type
       count(i) = count(i-1)+numIons(i-1);
   end
end

for i = 1: N_Mtype
   for j = 1: N_type
       s = sum(STDMOL(MtypeLIST(i)).types==j);
       if s > 0
          ID(count(j)+1: count(j)+s) = i;
          count(j) = count(j)+s;
       end
   end
end

for i = 1:N_atom-1
   for j = i+1:N_atom
       if ID(i)==ID(j)
          Intra_map(i,j) = 0;
          Intra_map(j,i) = 0;
       end
   end
end


