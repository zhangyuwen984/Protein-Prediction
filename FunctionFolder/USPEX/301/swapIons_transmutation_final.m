function [goodSon, goodFather, transNow, numIons, numBlocks] = swapIons_transmutation_final(Ind_No)

% USPEX Version 7.2.4
% ORG_STRUC.maxIons and ORG_STRUC.minIons are not activated yet, not sure if they are needed at all
% the algorithm: 
% 1. Determine how many blocks are removed and how many blocks are added
% 2. Determine which atoms are being transmutated and make them 'grey' (atom Type = 0)
% 3. Assign types to grey atoms, if there is excess of them - delete, if we need extra atoms to add - add them randomly
global POOL_STRUC
global ORG_STRUC

numBlocks = POOL_STRUC.POPULATION(Ind_No).numBlocks;
addBlocks = zeros(1, length(numBlocks));  % how many blocks to add
delBlocks = addBlocks;                    % how many blocks to remove
father = POOL_STRUC.POPULATION(Ind_No).COORDINATES;
goodFather = father;
numIons = POOL_STRUC.POPULATION(Ind_No).numIons;

% we set the max fraction of atoms to be transmutated, not the set number of it, thus we calculate the set number first
% howManyTrans == 0 means the user doesn't want to think for himself (we told him so in INPUT.txt) 

howManyTrans = round(ORG_STRUC.howManyTrans*sum(numBlocks));
if howManyTrans < 1 | isnan( howManyTrans )
    howManyTrans = 1;
end

% would be kind of boring always to swap the same amount of ions, so we introduce an element of randomness

transNow = RandInt(1,1,[1,howManyTrans]);

% sometimes transition doesn't exist! f.e.: 5 0 0 and only allowed transmutation is 2 <=> 3
% check whether structure can be transmuted, i.e. exist non zero block that is listed in allowed transmutations
if ~(ORG_STRUC.specificTrans == 0)
 validStr = 0;
 for i = 1 : length(numBlocks)
   if (~isempty(find(i == ORG_STRUC.specificTrans))) & (numBlocks(i) > 0)
    validStr = 1;
   end
 end
 if validStr == 0
  transNow = 0;
  goodSon = [];
 end
end

for any = 1 : transNow

  if ORG_STRUC.specificTrans == 0
     constrain = 1;
     while constrain  
     % choose an block randomly
      whichTypeFirst = RandInt(1,1,[1,length(numBlocks)]);
     % TO DO: check if there is a transmutation that satisfies constrains
      if numBlocks(whichTypeFirst) > 0
        constrain = 0;
      end
     end
  else
  % TO DO: add minIon and maxIon contrains!
     specified = 1;
     while specified
       whichTypeFirst = RandInt(1,1,[1,length(numBlocks)]);
       if (~isempty(find(whichTypeFirst == ORG_STRUC.specificTrans))) & (numBlocks(whichTypeFirst) > 0)
          [row, col] = find(whichTypeFirst == ORG_STRUC.specificTrans);
          acceptedTrans = [];
          for counting = 1 : length(row)
             if col(counting)-0.1 > 1
               acceptedTrans = cat(1,acceptedTrans,ORG_STRUC.specificTrans(row(counting),1:col(counting)-1));
             end
             if col(counting)+0.1 < size(ORG_STRUC.specificTrans,2)
                acceptedTrans = cat(1,acceptedTrans,ORG_STRUC.specificTrans(row(counting),col(counting)+1:end));
             end
          end
         specified = 0;
       end
     end
  end

  % choose a second block, different from the first one. 


  match = 1;
  while match
    whichTypeSecond = RandInt(1,1,[1,length(numBlocks)]);
    if ORG_STRUC.specificTrans == 0
      if (whichTypeFirst ~= whichTypeSecond) % &(numBlocks(whichTypeSecond) < ORG_STRUC.maxIons(ind))
        match = 0;
      end
    else
      if ~isempty(find(whichTypeSecond == acceptedTrans))
        match = 0;
      end
    end
  end

  delBlocks(whichTypeFirst) = delBlocks(whichTypeFirst) + 1;
  addBlocks(whichTypeSecond) = addBlocks(whichTypeSecond) + 1;
  numBlocks(whichTypeFirst) = numBlocks(whichTypeFirst) - 1;
  numBlocks(whichTypeSecond) = numBlocks(whichTypeSecond) + 1;

end

if transNow > 0   % valid transmutation existed
 addAtoms = addBlocks*ORG_STRUC.numIons; % how many atoms to add for every type
 delAtoms = delBlocks*ORG_STRUC.numIons; % how many atoms to remove from every type
 for i = 1 : length(numIons)
  if delAtoms(i) > numIons(i)
    addAtoms(i) = addAtoms(i) - (delAtoms(i) - numIons(i));
    delAtoms(i) = numIons(i);
  end
 end
 toDelete = zeros(1, sum(delAtoms));  % nomera atomov kotorye stanut 'serymi'
 toAdd = zeros(1, sum(addAtoms));     % atom TYPES(!) that will be added
 atomTypes = zeros(1, sum(numIons));  % atom types for final structure
 c = 0;
 for i = 1 : length(numIons)
  dummy = randperm(numIons(i)) + c;
  c = c + numIons(i);
  i1 = sum(delAtoms(1:i-1));
  toDelete(i1+1 : i1+delAtoms(i)) = dummy(1:delAtoms(i));
  i1 = sum(addAtoms(1:i-1));
  toAdd(i1+1 : i1+addAtoms(i)) = i;
  i1 = sum(numIons(1:i-1));
  atomTypes(i1+1 : i1+numIons(i)) = i;
 end
 % make atoms 'grey'
 atomTypes(toDelete(1:end)) = 0;
 % fill 'grey' atoms with the types we have to add
 dummy = randperm(length(toAdd));
 i1 = 1;
 for i = 1 : sum(numIons)
  if (atomTypes(i) == 0) & (i1 <= length(toAdd))
    atomTypes(i) = toAdd(dummy(i1));
    i1 = i1 + 1;
  end
 end

 % make a child and add atoms if needed
 numIons = numBlocks*ORG_STRUC.numIons;   % new POP_STRUC numIons
 goodSon = zeros(sum(numIons), 3);
 g = 1;
 if length(toDelete) >= length(toAdd)    % we just form the offspring ignoring atoms with atomType = 0
  for i = 1 : length(numIons)
   for j = 1 : length(atomTypes)
     if atomTypes(j) == i
       goodSon(g, 1:3) = father(j, 1:3);
       g = g + 1;
     end
   end
  end
 else   % we have to add a few atoms randomly
  extraAtoms = rand(length(toAdd) - length(toDelete),3);
  extraTypes = toAdd(dummy(i1:end)); % those that were left out
  for i = 1 : length(numIons)
   for j = 1 : length(atomTypes)
     if atomTypes(j) == i
       goodSon(g, 1:3) = father(j, 1:3);
       g = g + 1;
     end
   end
   for j = 1 : length(extraTypes)
     if extraTypes(j) == i
       goodSon(g, 1:3) = extraAtoms(j, 1:3);
       g = g + 1;
     end
   end
  end
 end
end
