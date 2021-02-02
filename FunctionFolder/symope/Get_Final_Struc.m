function [Lattice, candidate] = Get_Final_Struc(newLattice, candidate, permutationBack)

newLattice1 = latConverter(newLattice);
newLattice2 = newLattice1;
candidate1 = candidate;

for axis = 1 : 3
    newLattice2(axis) = newLattice1(permutationBack(axis));
    candidate(:,axis) = candidate1(:,permutationBack(axis));
end
Lattice = latConverter(newLattice2);
