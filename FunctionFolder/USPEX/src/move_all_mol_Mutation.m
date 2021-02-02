function [new_Coord]= move_all_mol_Mutation(MOLECULES, numMols, LATTICE, order)

% USPEX Version 7.3.5
% Change: introduced; atom mutation based on order; 
% all worst atoms are mutated, the degree of mutation decreases with order
% worst atom is mutated by applying gaussian with sigma = max_sigma, best atom - sigma = 0

global ORG_STRUC

max_sigma = ORG_STRUC.howManyMut;
N = sum(ORG_STRUC.numMols);
new_Lattice = LATTICE;

for ind = 1:sum(numMols)
    new_Coord(ind,:) =  MOLECULES(ind).ZMATRIX(1,:);
end
    new_Coord=new_Coord/new_Lattice;

if length(new_Lattice) == 6
 new_Lattice = latConverter(new_Lattice);
end
temp_potLat = latConverter(new_Lattice);

[junk, ranking] = sort(order); %junk = POP...order(ranking)
r1 = order(ranking(1));
rN = order(ranking(N));

if rN > r1
 for i = 1 : N
%  deviat_dist = randn(3,1)*max_sigma*(i-1)/(N-1);
  rI = order(ranking(i));
  koef = (rN-rI)/(rN-r1);
  deviat_dist = randn(3,1)*max_sigma*koef;

  new_Coord(ranking(i),1) = new_Coord(ranking(i),1) + deviat_dist(1)/temp_potLat(1);
  new_Coord(ranking(i),2) = new_Coord(ranking(i),2) + deviat_dist(2)/temp_potLat(2);
  new_Coord(ranking(i),3) = new_Coord(ranking(i),3) + deviat_dist(3)/temp_potLat(3);

  new_Coord(ranking(i),1) = new_Coord(ranking(i),1) - floor(new_Coord(ranking(i),1));
  new_Coord(ranking(i),2) = new_Coord(ranking(i),2) - floor(new_Coord(ranking(i),2));
  new_Coord(ranking(i),3) = new_Coord(ranking(i),3) - floor(new_Coord(ranking(i),3));
 end
end
