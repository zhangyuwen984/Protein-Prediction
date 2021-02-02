function atoms = sort_order(order, array, sort_type)

% USPEX 7.3.1
% sort type = 'descend' or 'ascend', sorts the number in array (linear tournament) according to order 

L = length(array);
tournament = zeros(L,1);
% the following tournament is linear
tournament(L) = 1;
for loop = 2:L
   tournament(end-loop+1)= tournament(end-loop+2)+loop;
end

[junk, ranking] = sort(order(array)); %junk = POP...order(tempInd(ranking))

if strcmp(sort_type, 'descend')
 tmp = ranking;
 for i = 1:length(ranking)
  ranking(end-i+1) = tmp(i);
 end
end

atom_rank = randperm(max(tournament));
atoms = zeros(length(tournament),1);
chosen = zeros(L,1);
j1 = 1;
for i1 = 1:max(tournament)
  at1 = find (tournament > (atom_rank(i1)-1));
  if chosen(ranking(at1(end))) == 0
      chosen(ranking(at1(end))) = 1;
      atoms(j1) = ranking(at1(end));
      j1 = j1 + 1;
  end
end
