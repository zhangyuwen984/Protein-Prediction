function [numIons, offspring,potentialLattice,fracFrac,dimension,offset,fracLattice] = heredity_cluster(parent1,parent2,parentlat1,parentlat2, order1, order2, fracFrac)
% This function looks terribly complicated - but don't worry, it's all quite straight forward.
% Keep in mind that we are dealing with a periodical problem (periodicity : 1)
% => example (coordinate) : 2.41 = 1.41 = 0.41 = -0.59 = - 1.59...

global POP_STRUC
global ORG_STRUC
ordering_on = ORG_STRUC.ordering;
numIons = ORG_STRUC.numIons;

fracFrac = 0.49999999;
fracLattice = fracFrac;

% We select randomly a dimension used for the spacial criteria of the heredity
dimension = RandInt(1,1,[1,3]);
%%%%% CHANGE 8.0.6 - rotate both cluster for the same angle, around the same axis
% rotate around main inertia axis only
offset = pi*(rand(3,1) - 0.5);
[parentlat1, parent1] = rotateCluster(parent1, parentlat1, 0, 0, offset(1,1));
[parentlat2, parent2] = rotateCluster(parent2, parentlat2, 0, 0, offset(1,1));

if rand(1) > ORG_STRUC.percSliceShift   % rotate around main inertia axis only
  % done above
else % general rotation
  [parentlat1, parent1] = rotateCluster(parent1, parentlat1, offset(2,1), offset(3,1), 0);
  [parentlat2, parent2] = rotateCluster(parent2, parentlat2, offset(2,1), offset(3,1), 0);
end

% find the indices of all ions located in the respective fraction of each structure
whichIons1 = find(parent1(:,dimension)<fracFrac);
whichIons2 = find(parent2(:,dimension)>=fracFrac);

correlation_coefficient = ORG_STRUC.correlation_coefficient;
cor_dir = ORG_STRUC.cor_dir;

% NOTE: unlike for crystals, we rotate both clusters in the same way and choose the most ordered part of the random one
ord_parent = rand; % < 0.5 - 1st parent order matters, >= 0.5 - 2nd -//-

if ord_parent < 0.5 
 nAt = size(parent1,1);
 Lchar = 0.5*power(abs(det(parentlat1))/nAt, 1/3); % characteristic length
else
 nAt = size(parent2,1);
 Lchar = 0.5*power(abs(det(parentlat2))/nAt, 1/3); % characteristic length
end

% for clusters:  L = Lchar*round((6*pi*pi*N)^1/3) and we divide the Nmax by 2
Nmax = 0.5*round((6*pi*pi*nAt)^1/3);
Nslabs = round(Nmax/(1+(Nmax-1)*(cos(correlation_coefficient*pi/2))^2));
%disp(['Number of cuts = ' num2str(Nslabs)])

% chooses the more/less ordered slab, depending on degree of (anti)correlation
if ordering_on
 if length(whichIons1) > 0
  ord1 = sum(order1(whichIons1))/length(whichIons1);
 elseif cor_dir <= 0
  ord1 = 0;
 else
  ord1 = 1;
 end
 if length(whichIons2) > 0
  ord2 = sum(order2(whichIons2))/length(whichIons2);
 elseif cor_dir <= 0
  ord2 = 0;
 else
  ord2 = 1;
 end

 lat1_candidate = parentlat1;
 lat2_candidate = parentlat2;
 parent1_candidate = parent1;
 parent2_candidate = parent2;

 for i = 1 : Nslabs
  offset_extra = pi*(rand(3,1) - 0.5);
  [lat1_tmp, parent1_tmp] = rotateCluster(parent1, parentlat1, offset_extra(1,1),offset_extra(2,1),offset_extra(3,1));
  [lat2_tmp, parent2_tmp] = rotateCluster(parent2, parentlat2, offset_extra(1,1),offset_extra(2,1),offset_extra(3,1));
  whichIons1_extra = find(parent1_tmp(:,dimension)<fracFrac);
  whichIons2_extra = find(parent2_tmp(:,dimension)>=fracFrac);
  if length(whichIons1_extra) > 0
   ord1_extra = sum(order1(whichIons1_extra))/length(whichIons1_extra);
  elseif cor_dir <= 0
   ord1_extra = 0;
  else
   ord1_extra = 1;
  end
  if length(whichIons2_extra) > 0
   ord2_extra = sum(order2(whichIons2_extra))/length(whichIons2_extra);
  elseif cor_dir <= 0
   ord2_extra = 0;
  else
   ord2_extra = 1;
  end
  if ord_parent < 0.5 % check order of the 1st cluster
    if ((ord1 > ord1_extra) & (cor_dir <= 0)) | ((ord1 < ord1_extra) & (cor_dir == 1))
   %  do nothing, worse cut
    else
     whichIons1 = whichIons1_extra;
     whichIons2 = whichIons2_extra;
     ord1 = ord1_extra;
     ord2 = ord2_extra;
     lat1_candidate = lat1_tmp;
     lat2_candidate = lat2_tmp;
     parent1_candidate = parent1_tmp;
     parent2_candidate = parent2_tmp;
    end
  else % choose 2nd cluster as parent
    if ((ord2 > ord2_extra) & (cor_dir <= 0)) | ((ord2 < ord2_extra) & (cor_dir == 1))
   %  do nothing, worse cut
    else
     whichIons1 = whichIons1_extra;
     whichIons2 = whichIons2_extra;
     ord1 = ord1_extra;
     ord2 = ord2_extra;
     lat1_candidate = lat1_tmp;
     lat2_candidate = lat2_tmp;
     parent1_candidate = parent1_tmp;
     parent2_candidate = parent2_tmp;
    end
  end
 end

 parentlat1 = lat1_candidate;
 parentlat2 = lat2_candidate;
 parent1 = parent1_candidate;
 parent2 = parent2_candidate;

end

% prepare the variable offspring. This needs to be done due to the command 'cat' used later on
offspring = zeros(0,3);
ionCount = zeros(length(ORG_STRUC.atomType),1);

for ind = 1 : length(numIons)
  ionChange(ind) = sum(numIons(1:ind));
end
ionCh1 = zeros(1,length(ionChange)+1);
ionCh1(2:end) = ionChange;
ionCh2 = ionCh1;

for ind = 1:length(ORG_STRUC.atomType)
    % this loop is for all the different ion types.

    candidates = [];

    % find the indices of the ions for the type of ion done in this turn of
    % the loop
    smaller1 = find(whichIons1<=ionCh1(ind+1));
    smaller2 = find(whichIons2<=ionCh2(ind+1));

    ions1  = find(whichIons1(smaller1)>ionCh1(ind));
    ions2  = find(whichIons2(smaller2)>ionCh2(ind));

    % total number of this type of ion
    ionCount(ind,1) = length(ions1) + length(ions2);

    % then we have to fill up with further ions
    if ionCount(ind,1) < numIons(ind)

        % calculate how many ions (of this type) are missing
        howmany = numIons(ind) - ionCount(ind,1);
        % determine how many ions are to be filled into each fraction probabilistically.
        % Linear probability of each ion to be of a fraction (linear with the size of the fraction)
        suppleIons_1 = length(find(rand(howmany,1)<fracFrac));
        suppleIons_2 = howmany-suppleIons_1;
        
        % for each fraction we choose the determined amount of ions. The new ions are chosen (randomly) 
        % from the structure the ions in this fraction weren't originally chosen from
        if suppleIons_1 > 0
           tempInd = find(parent2(ionCh1(ind)+1:ionCh1(ind+1),dimension)<fracFrac);
           
           if cor_dir <= 0  % low fitness (=good) <=> high order
             atoms = sort_order(order2, tempInd, 'descend');
           else % low fitness (=good) <=> low order
             atoms = sort_order(order2, tempInd, 'ascend');
           end
           candidates = parent2(ionCh2(ind) + tempInd(atoms(1:suppleIons_1)) , :);
        end
        if suppleIons_2 > 0
            tempInd = find(parent1(ionCh1(ind)+1:ionCh1(ind+1),dimension)>fracFrac);
            if cor_dir <= 0  % low fitness (=good) <=> high order
              atoms = sort_order(order1, tempInd, 'descend');
            else % low fitness (=good) <=> low order
              atoms = sort_order(order1, tempInd, 'ascend');
            end
            candidates = cat(1,candidates,parent1(ionCh1(ind) + tempInd(atoms(1:suppleIons_2)) , :));
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if on the other hand the number of ions we have is too big, we need to eliminate a few
    elseif  ionCount(ind,1)>numIons(ind)

        howmany = ionCount(ind,1) - numIons(ind);
        % we determine the amount to be eliminated from each fraction
        % probabilisticly (probability of each ion to be of a certain 
        %fraction linear with the size of this fraction)
       
        deleteIons_1 = length(find(rand(howmany,1)<fracFrac));
        deleteIons_2 = howmany - deleteIons_1;

        % quick safeguarding against trouble
        if length(ions1) < deleteIons_1
            oops = deleteIons_1 - length(ions1);
            deleteIons_1 = length(ions1);
            deleteIons_2 = deleteIons_2+oops;
        elseif length(ions2) < deleteIons_2
            oops = deleteIons_2 - length(ions2);
            deleteIons_2 = length(ions2);
            deleteIons_1 = deleteIons_1+oops;
        end
        % now we delete the surplus ions (randomly within each fraction) 
        for xy = 1:deleteIons_1
          ions1(RandInt(1,1,[1,length(ions1)])) = [];
        end
        for xy = 1:deleteIons_2
          ions2(RandInt(1,1,[1,length(ions2)])) = [];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fill the selected ions coordinates into the offspring
    addOn = cat(1,parent1(whichIons1(smaller1(ions1(:))) , :), parent2(whichIons2(smaller2(ions2(:))) , :));
    addOn = cat(1,addOn,candidates);
    offspring = cat(1,offspring,addOn);
end

  temp_potLat = fracFrac*parentlat1 + (1-fracFrac)*parentlat2;
  volLat = det(temp_potLat);
  % scale the lattice to the volume we asume it approximately to be
  latVol = det(parentlat1)*fracFrac + det(parentlat2)*(1-fracFrac);
  ratio = latVol/volLat;
  temp_potLat = latConverter(temp_potLat);
  temp_potLat(1:3)= temp_potLat(1:3)*(ratio)^(1/3);
  potentialLattice = latConverter(temp_potLat);
