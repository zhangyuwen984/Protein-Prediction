function [permutation, permutationBack] = GetPermutation(nsym);

% for monoclinic and orthorhombic groups we do lattice vector exchange/swap
% so we can generate non conventional symmetries too (P2mm, Pm2m, etc)
% they are changed back after Stoke's code generates the 'conventional' structure
% for monoclinic only a and c are swapped, for orthorombic - any combination
if (nsym > 2) & (nsym < 16) % monoclinic
    if rand > 0.5
        permutation = [3 2 1]; % only swap x and z coordinates
        permutationBack = permutation;
    else
        permutation = [1 2 3]; % leave as is
        permutationBack = permutation;
    end
elseif (nsym > 15) & (nsym < 75) % orthorhombic
    ListPerm = [1 2 3; 3 2 1; 1 3 2; 2 1 3; 3 1 2; 2 3 1];  %possible permutation;
    id = RandInt(1,1,6)+1;
    permutation = ListPerm(id,:); % all swaps allowed
    if id == 5
        permutationBack = ListPerm(6,:);  % exchange is not reversible (3 1 2) to (2 3 1)
    elseif id == 6
        permutationBack = ListPerm(5,:);  % exchange is not reversible (3 1 2) to (2 3 1)
    else
        permutationBack = permutation;
    end
else
    permutation = [1 2 3];   % no swaps
    permutationBack = permutation;
end
