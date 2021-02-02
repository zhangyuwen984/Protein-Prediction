function [goodSon, swapNow] = swapIons_mutation_final(father,numIons,order)

% Completely Rewritten

global ORG_STRUC

specificSwaps = ORG_STRUC.specificSwaps;
 howManySwaps = ORG_STRUC.howManySwaps;
     atomType = ORG_STRUC.atomType;

ionChange = zeros(1,length(numIons)+1);
for ind = 1:length(numIons)
  ionChange(ind+1) = sum(numIons(1:ind));
end    

% would be kind of boring always to swap the same amount of ions, so we
% introduce an element of randomness
swapNow = RandInt(1,1,[1,howManySwaps]);
atom1_old = 0;
atom2_old = 0;

for any = 1:swapNow
    
    if specificSwaps==0  %randomly pick
        whichType = RandInt(1,2,[1,length(atomType)]);
    else
        ID = RandInt(1,1,[1,size(specificSwaps, 1)]);
        whichType = specificSwaps(ID,:);
    end
    if numIons(whichType(1))==0 | numIons(whichType(2))==0
        swapNow = 0;
        break;
    else
       order1 = order(ionChange(whichType(1))+1:ionChange(whichType(1)+1));
       order2 = order(ionChange(whichType(2))+1:ionChange(whichType(2)+1));
       [atoms1,ranking1] = sort(order1, 'ascend');
       [atoms2,ranking2] = sort(order2, 'ascend');
       ranking1   = ranking1 + ionChange(whichType(1));
       ranking2   = ranking2 + ionChange(whichType(2));
   
       if rand()>0.5
          atom1  = ranking1(Roulette(length(order1)));
          atom2  = ranking2(Roulette(length(order2)));
       else
          atom1  = ranking1(1); %just chose the atom with lowest order
          atom2  = ranking2(1);
       end
   
       if atom1 ~= atom1_old  %make sure we don't select the same atoms twice
          % swap the coordinates
          tempFather = father(atom1,:) ;
          father(atom1,:) = father(atom2,:);
          father(atom2,:) = tempFather;
          atom1_old = atom1;
          atom2_old = atom2;
       end
   end
end

goodSon = father;


function select = Roulette(N)
tournament = zeros(N,1);
tournament(N) = 1;
for loop = 2:N
    tournament(end-loop+1) = tournament(end-loop+2) + loop^2;
end
tmp = find (tournament>RandInt(1,1,[0,max(tournament)-1]));
select = tmp(end);
