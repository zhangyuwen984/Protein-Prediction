function [newcoor, newsymbol, newcharge, newtypes, newbound, newPair] = resortXYZ(coor, symbol, charge, types, bound, Pair, first, radii);
%This function is to resort the atoms according to the neighbor list
%The atoms are located based on the sequence in the neighbor list.
%INPUT: 
%coor (N*3) -- atomic coordinate;
%radii(N*1) -- covalant radii;
%symbol(N*1)-- atomic symbol (C, N, O, .etc)
%Pair(N*7) --  The neighbor list (1-6 records neighbor atom, 7 gives the coord. number)
%The idea is that we build a graph according to the neighbor list.
%And then include the atoms layer by layer.
%  1    
% 3 4
%2 9 8
%.....
N_max = size(Pair,2);
newcoor = coor;
newList = [];
newsymbol = symbol;
n_atom = size(coor,1);
%to choose 1st atom from which direction
if first == 0
   dist = zeros(3,1);
   dist(1) = max(coor(:,1)) - min(coor(:,1));
   dist(2) = max(coor(:,2)) - min(coor(:,2));
   dist(3) = max(coor(:,3)) - min(coor(:,3));
   [tmp1, axes]=max(dist);
   [tmp, first]=min(coor(:,axes));
end

newList(1) = first;  %target table -- we start to add the atom from 1
count = 1;       %how many atoms in the table --- we start to count the atom from 1
layer_num   = 1; %how many atoms in the current layer -- we start from 1st layer
for i = 1: n_atom  %Number of layers; maximum layer = n_atom
    count0 = count;
    for j = (count-layer_num+1):count  %Number of atoms in the current layer
        ID = newList(j);
        for k = 1:Pair(ID,N_max) %Number of net atoms for the given atom
            [count, newList] = updateList(newList, Pair(ID,k));
        end
    end
    layer_num = count - count0;
    if layer_num ==0 %If the new layer does not have any new atoms, we stop
       disp('Search is complete');
       break;
    end
end

if count < n_atom
   disp('The atoms are not fully connected');
   quit
end
%update everything
for i = 1:n_atom 
     newcoor(i,:) =   coor(newList(i),:);
   newcharge(i,:) = charge(newList(i));
    newtypes(i,:) =  types(newList(i));
    if ~isempty(symbol)
        newsymbol(i)   = symbol(newList(i));
    end
    newradii(i)   =  radii(newList(i));
    newbound(i)   =  bound(newList(i));
end
%update the Pair as well
newPair = find_pair(newcoor, newradii)

%If the atom is not in the list, add it
function [count, newList]=updateList(newList, ID)
new = 1;
for k = 1:length(newList)
    if ID == newList(k)
       new = 0;
       break;
    end
end
if new == 1
   newList=[newList; ID];
end
count = length(newList);
