function [new_Coord] = move_all_atom_Mutation(coord, numIons, lattice, order, max_sigma)

% USPEX Version 8.0.0
% Change: removed dependance on POP_STRUC to use in multistage mutation
% atom mutation based on order: all worst atoms are mutated, the degree of mutation decreases with order
% worst atom is mutated by applying gaussian with sigma = max_sigma, best atom - sigma = 0
% 8.7.1 - added renormalization so that sum of all displacements is zero

N = sum(numIons);
new_Coord = coord;

if length(lattice) == 6
 lattice = latConverter(lattice);
end
temp_potLat = latConverter(lattice);

[junk, ranking] = sort(order);  % junk = order(ranking)
r1 = order(ranking(1));
rN = order(ranking(N));

if rN > r1
 % to make sum equal to 0, we divide the summary vector proportionally between displacements and then deduct it
 % we also try to avoid too big summary displacement in the first place
 tmp = 0;
 while tmp < 1000
   sum_dist = zeros(1,3);
   sum_norm = 0;
   deviat_dist = zeros(N,3);
   for i = 1 : N
     rI = order(ranking(i));
     koef = (rN-rI)/(rN-r1);
     deviat_dist(i,:) = randn(3,1)*max_sigma*koef;
     sum_dist = sum_dist + deviat_dist(i,:);
     sum_norm = sum_norm + norm(deviat_dist(i,:));
   end
   if norm(sum_dist) <= 1.5*max_sigma % displacement is small, will be corrected
     break;
   else
     tmp = tmp + 1;
   end
 end
 for i = 1 : N
  deviat_dist(i,:) = deviat_dist(i,:) - sum_dist*norm(deviat_dist(i,:))/sum_norm;
 end 
 for i = 1 : N
  new_Coord(ranking(i),1) = coord(ranking(i),1) + deviat_dist(i,1)/temp_potLat(1);
  new_Coord(ranking(i),2) = coord(ranking(i),2) + deviat_dist(i,2)/temp_potLat(2);
  new_Coord(ranking(i),3) = coord(ranking(i),3) + deviat_dist(i,3)/temp_potLat(3);

  new_Coord(ranking(i),1) = new_Coord(ranking(i),1) - floor(new_Coord(ranking(i),1));
  new_Coord(ranking(i),2) = new_Coord(ranking(i),2) - floor(new_Coord(ranking(i),2));
  new_Coord(ranking(i),3) = new_Coord(ranking(i),3) - floor(new_Coord(ranking(i),3));
 end
else
 new_Coord = coord;
end
