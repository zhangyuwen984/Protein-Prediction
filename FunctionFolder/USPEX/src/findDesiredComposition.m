function [numIons, numBlocks] = findDesiredComposition(maxBlocks, blocks, child)

% USPEX 8.4.2 - created
% to finds a composition that requires the least addition/deleting of atoms from child
% one should be able to compose numIons out of blocks
% numBlocks - full amount of blocks from both parents, determines which atoms could be used to make a child

% blocks = ORG_STRUC.numIons;    just FYI
maxAtoms = maxBlocks*blocks;
maxAdded = maxAtoms - child; % how many atoms one could, possibly, add
Nb = size(blocks,1); % N of blocks
Nt = size(blocks,2); % N of atom types

tmp = 1; % calculate how many different block combinations we have to check with exhaustive search
for i = 1 : Nb
 min1 = max(max(maxAtoms));
 for j = 1 : Nt
  if blocks(i,j) > 0
   if min1 > maxAtoms(j)/blocks(i,j)
    min1 = floor(maxAtoms(j)/blocks(i,j));
   end
  end
 end
 % now 'min1' says us how many times block 'j' can be fitted into maxAtoms
 if min1 > 0
  tmp = tmp*min1;
 end
end

% now we do 'greedy' algorithm for big amount of combinations and exhaustive search for small one
% 'greedy' algorithm won't give the best answer, but it should be good in most cases
% we repeat greedy algorithm 10 times and choose the best answer

if tmp > 0   % greedy algorithm
 bestGreed = sum(maxAtoms);
 for g = 1 : 20

   tolerance = round(rand*2);    % maximum amount of atoms that could be added to child for any specific type

   blockN = zeros(1,Nb); % specifies the number of blocks in the composition found by algorithm
   child_tmp = child;
   crutch = rand(Nb,1);
   [junk, blockOrder] = sort(crutch);
   for i = 1 : Nb
     block = blocks(blockOrder(i), 1:end);
     ind = 1;
     min1 = sum(child_tmp);
     while 1
       if min(child_tmp - ind*block) < -1*tolerance
         break;
       end
% take into account maxAtoms, sometimes we can't add atoms at all for varcomp, since all atoms of specific type could be ALREADY in the child
       if min(maxAdded + (child_tmp - ind*block)) < 0
         break;
       end
       if min1 > sum(abs(child_tmp - ind*block))
         min1 = sum(abs(child_tmp - ind*block));
         blockN(blockOrder(i)) = ind;
       end
       ind = ind + 1;
     end
     child_tmp = child_tmp - blockN(blockOrder(i))*block;
   end
   if bestGreed > sum(abs(child - blockN*blocks))
     bestGreed = sum(abs(child - blockN*blocks));
     numBlocks = blockN;
   end
 end
 numIons = numBlocks*blocks;
else  % exhaustive search
end
