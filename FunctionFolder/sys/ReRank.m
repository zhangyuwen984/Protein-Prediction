function ReRank()

global POP_STRUC

TEMP_STRUC = POP_STRUC;
N = length(POP_STRUC.POPULATION);
N_atom = zeros(N,1);
for i = 1:N
    N_atom(i) = sum(POP_STRUC.POPULATION(i).numIons);
end
[nothing, ranking] = sort(N_atom, 'descend');
for i = 1:N
    POP_STRUC.POPULATION(i) = TEMP_STRUC.POPULATION(ranking(i));
end

