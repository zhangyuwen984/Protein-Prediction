function [fitness] = final_convex_hull(convex_hull, enthalpies, blocks, numIons)

N_P = size(blocks, 1);  % number of different Structures
N_T = size(blocks, 2);  % number of different blocks
N_atomType = size(numIons,2);
fitness = zeros(N_P,1);
convex_hull0 = Transform_Convex_Hull(convex_hull, numIons, 1);  % count by atom
convex_hull_block = convex_hull0(:,1:N_T);
convex_hull_E     = convex_hull0(:,N_T+1);
for i = 1 : N_P
    N_atom  = sum(blocks(i,:)*numIons);
    N_Block = sum(blocks(i,:));
    E = enthalpies(i)*N_atom/N_Block;
    [tmp, samecomp, fitness(i)] = CheckDecomposition(convex_hull_block, ...
        convex_hull_E, blocks(i,:), E);
end
% very important and often neglected step;
% make sure precision is e-4 since it's compared with fitness of that precision!
fitness = round(fitness*10000)/10000;

