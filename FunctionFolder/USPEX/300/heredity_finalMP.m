function [numIons, offspring,potentialLattice,fracFrac,dimension,offset,fracLattice] = heredity_finalMP(parents, dimension)

% USPEX Version 7.4.2
% Change: order added

% The heredity was designed to allow structural heritage. 
% Keep in mind that we are dealing with a periodical problem (periodicity : 1)
% => example (coordinate) : 2.41 = 1.41 = 0.41 = -0.59 = - 1.59...

global POP_STRUC
global ORG_STRUC

ordering_on = ORG_STRUC.ordering;

% find optimal number of parents

% following constrains help to counter the lattice, that you can't optimally split. F.e. let mediane for split be 2A, 2 parents
% gave 5.2A lattice vector along the split size and 3 parents gave 4.7A lattice vector. In this case both variables listed below
% will be equal to one and we choose the best of these non optimal splits. Thus avoiding a possible infinite cycle

soft_const_low = 0; % =1 when you could satisfy constrain for the first time, but split will be above (mediane + mediane/2)
soft_const_high = 0; % =1 when you could satisfy constrain and are under (mediane - mediane/2) split size for the first time
minSlice = ORG_STRUC.minSlice;
maxSlice = ORG_STRUC.maxSlice;
medSlice = (maxSlice+minSlice)/2;
optim = 1000000; % deviation of best split from optimal
tempN = 0;


flag = 0;
Np1 = 1;


while flag == 0
  
  Np1 = Np1 + 0.5;
  Np = round(Np1);
  % fracLat determines the weights of the parent lattices in child lattice

  n = Np; 
  tmp = 0;
  fracLat = zeros(1, Np);
  for i = 1 : (Np-1)
   r = rand(1);
   fracLat(i) = tmp + (1-tmp)*(1/(2*n) + r/n); % equivalent to random(1/n - 1/2n, 1/n + 1/2n)
   tmp = fracLat(i);
   n = n-1;
  end
  fracLat(Np) = 1;

  if ~ORG_STRUC.constLattice

      %%%% represantation problem!
      temp_potLat = fracLat(1)*POP_STRUC.POPULATION(parents(1)).LATTICE;
      for i = 2 : Np 
        temp_potLat = temp_potLat + (fracLat(i)-fracLat(i-1))*POP_STRUC.POPULATION(parents(i)).LATTICE;
      end;
      volLat = det(temp_potLat);
      if sign(volLat) == -1
           temp_potLat = -1*temp_potLat;
      end
    
      % scale the lattice to the volume we asume it approximately to be

      latVol = fracLat(1)*POP_STRUC.POPULATION(parents(1)).Vol(end);
      for i = 2 : Np 
        latVol = latVol + (fracLat(i)-fracLat(i-1))*POP_STRUC.POPULATION(parents(i)).Vol(end);
      end;
      ratio = latVol/volLat;
      temp_potLat = latConverter(temp_potLat);
      temp_potLat(1:3)= temp_potLat(1:3)*(ratio)^(1/3);
      potentialLattice = latConverter(temp_potLat);
    
  else
     potentialLattice = ORG_STRUC.lattice;
     temp_potLat = ORG_STRUC.lattice;
  end

  % convert lattice to convinient format if it's 3x3 matrix
  if size(temp_potLat,1) == 3 
     temp_potLat = latConverter(temp_potLat);
  end

  if (temp_potLat(dimension) > (Np-0.6)*medSlice) & (temp_potLat(dimension) < (Np+0.6)*medSlice)
    flag = 1; % Np found!
    break;
  elseif (temp_potLat(dimension) < Np*maxSlice) & (temp_potLat(dimension) > Np*minSlice) % hard constrains OK, but not optimal split
    optim1 = abs(temp_potLat(dimension)-Np*medSlice);
    if optim1 < optim
      optim = optim1;
      tempN = Np;
      tempLat = temp_potLat;
    end
% now we add additional safeguard to reduce the number of 'quits' - try to vary Np if we can't satisfy constrains
  elseif (temp_potLat(dimension) > Np*maxSlice) & (Np < length(parents))
    if (temp_potLat(dimension) < (Np+1)*maxSlice) & (temp_potLat(dimension) > (Np+1)*minSlice) % hard constrains OK, but not optimal split
     optim1 = abs(temp_potLat(dimension)-(Np+1)*medSlice);
     if optim1 < optim
      optim = optim1;
      tempN = Np+1;
      tempLat = temp_potLat;
     end
    end
  elseif (temp_potLat(dimension) < Np*minSlice) & (Np > 2)
    if (temp_potLat(dimension) < (Np-1)*maxSlice) & (temp_potLat(dimension) > (Np-1)*minSlice) % hard constrains OK, but not optimal split
     optim1 = abs(temp_potLat(dimension)-(Np-1)*medSlice);
     if optim1 < optim
      optim = optim1;
      tempN = Np-1;
      tempLat = temp_potLat;
     end
    end
  end

  if (Np == length(parents)) & (tempN == 0)
       error = 'Not enough parent to satisfy slice constraints'
       dimension
       minSlice
       maxSlice
       temp_potLat
       fracLat
       parents
      for i=1:Np 
%        POP_STRUC.POPULATION(parents(i)).LATTICE
      end
       save ([ POP_STRUC.resFolder '/EEERRROOOR_heredity.txt'],'error');
       quit;
  end

  if (Np == length(parents)) 
   Np = tempN;
   temp_potLat = tempLat;
   potentialLattice = latConverter(temp_potLat);
   flag = 1; % best non optimal Np found
  end

end % of while


% fracFrac determines the splits, where spatial fraction of the each parent will be taken.
% Thus: 0.25, 0.75, 1 means one spatial fourth versus half versus one fourth.

n = Np; 
tmp = 0;
fracFrac = zeros(1, Np);
for i = 1 : (Np-1)
  r = rand(1);
  fracAdd = (1-tmp)*(1/(2*n) + r/n); % equivalent to random(1/n - 1/2n, 1/n + 1/2n);
  if fracAdd*temp_potLat(dimension) > maxSlice
   fracAdd = maxSlice/temp_potLat(dimension);
  elseif fracAdd*temp_potLat(dimension) < minSlice
   fracAdd = minSlice/temp_potLat(dimension);
  end
  fracFrac(i) = tmp + fracAdd;
  tmp = fracFrac(i);
  n = n-1;
end
fracFrac(Np)=1;

% in this way the most extreme difference between different slices is 1/Np
% check if the last slice is good and change other slices to make it good if needed


% In the following lines we simply change from colomn vectors to matrices,
% where each colomn represents the coordinates on one lattice vector
% (it makes life easier)

nIonsMax = 0;
for i = 1 : Np
 if sum(POP_STRUC.POPULATION(parents(i)).numIons) > nIonsMax 
   nIonsMax = sum(POP_STRUC.POPULATION(parents(i)).numIons);
 end
end

parent = -1*ones(nIonsMax,3,Np);

for i2 = 1 : 3
 for i3 = 1 : Np
  for i1 = 1 : sum(POP_STRUC.POPULATION(parents(i3)).numIons)
   parent(i1,i2,i3) = POP_STRUC.POPULATION(parents(i3)).COORDINATES(i1,i2);
  end
 end
end

if rand(1) > ORG_STRUC.percSliceShift
 offset = rand(Np,1);
 if ORG_STRUC.manyParents == 3  % all structures representing a parent should have same offsets
  offset(end) = offset(1);
  for i = 2 : floor(Np/2)
   offset(2*i-1) = offset(1);
   offset(2*i) = offset(2);
  end
 end
 for i = 1 : Np
    parent(:,dimension,i) = parent(:,dimension,i)+offset(i,1);
    parent(:,dimension,i) = parent(:,dimension,i) - floor(parent(:,dimension,i));
 end
else
 % The offset allows
 offset = rand(3*Np,1);
 if ORG_STRUC.manyParents == 3  % all structures representing a parent should have same offsets
  offset(end) = offset(3);
  offset(end-1) = offset(2);
  offset(end-2) = offset(1);
  for i = 2 : floor(Np/2)
   offset(6*i-5) = offset(1);
   offset(6*i-4) = offset(2);
   offset(6*i-3) = offset(3);
   offset(6*i-2) = offset(4);
   offset(6*i-1) = offset(5);
   offset(6*i) = offset(6);
  end
 end
 for i = 1 : Np
    parent(:,1,i)= parent(:,1,i)+offset((i-1)*3+1,1);
    parent(:,2,i)= parent(:,2,i)+offset((i-1)*3+2,1);
    parent(:,3,i)= parent(:,3,i)+offset((i-1)*3+3,1);
 end
 % IMPORTANT :  The reason for this translation is, that the localization of substructures change by inheritance (inhibiting convergence)
 parent = parent - floor(parent);
end


% for varcomp - fill the unused atoms with -1
for i2 = 1 : 3
 for i3 = 1 : Np
  for i1 = sum(POP_STRUC.POPULATION(parents(i3)).numIons)+1 : nIonsMax
   parent(i1,i2,i3) = -1;
  end
 end
end

% Nslabs = ORG_STRUC.Nslabs;
% take into account (anti)correlation of order vs energy:

j = 1;
FvsO = zeros(1,2);
for i = 1 : length(POP_STRUC.POPULATION)
  if (POP_STRUC.fitness(i) < 99999) & (POP_STRUC.DoneOrder(i) > 0)
    a_o = sum(POP_STRUC.POPULATION(i).order)/sum(POP_STRUC.POPULATION(i).numIons);
    if j == 1
     FvsO(1,1) = a_o;
     FvsO(1,2) = POP_STRUC.fitness(i);
    else
     FvsO = vertcat(FvsO, [a_o POP_STRUC.fitness(i)]);
    end
    j = j + 1;
  end
end

if exist('corrcoef')
    corr_tmp = corrcoef(FvsO);
else
    corr_tmp = 0;
end
correlation_coefficient = corr_tmp(1,2);
disp(' ')
%disp(['Correlation coefficient = ' num2str(correlation_coefficient)])
cor_dir = sign(corr_tmp(1,2)); % > 0 - correlation, < 0 - anticorrelation

% determine the thickness of the unit cell while cutting in a given direction (not equal to the lattice vector for non-cubic cells!)
L_cut = zeros(1,Np);
Lchar = zeros(1,Np);
Nslabs = zeros(1,Np);
for i = 1 : Np
  lat1 = POP_STRUC.POPULATION(parents(i)).LATTICE;
  if dimension == 1 
   L_cut(i) = abs(det(lat1)/norm(cross(lat1(2,:),lat1(3,:))));
  elseif dimension == 2
   L_cut(i) = abs(det(lat1)/norm(cross(lat1(1,:),lat1(3,:))));
  else
   L_cut(i) = abs(det(lat1)/norm(cross(lat1(1,:),lat1(2,:))));
  end
  Lchar(i) = 0.5*power(abs(det(lat1))/sum(POP_STRUC.POPULATION(parents(i)).numIons), 1/3); % characteristic length - approximate 'radius' of the atom in the cell = 0.5*(V/N)^1/3
  Nslabs(i) = round(L_cut(i)/(Lchar(i)+(L_cut(i)-Lchar(i))*(cos(corr_tmp(1,2)*pi/2))^2));
end
disp(['Number of slabs = ' num2str(Nslabs)])
%number_of_slabs = Nslabs

% find the indices of all ions located in the right respective fraction of
% each structure

ord = zeros(1,Np);

whichIons = zeros(nIonsMax,Np);
whichIonsN = zeros(Np,1);
whichI = find((parent(:,dimension,1)) < fracFrac(1) & (parent(:,dimension,1) > -0.5));
whichIonsN(1,1) = length(whichI);
for j = 1:length(whichI)
   whichIons(j,1) = whichI(j);
end
if ordering_on 
 if whichIonsN(1,1) > 0
   ord(1) = sum(POP_STRUC.POPULATION(parents(1)).order(whichI))/whichIonsN(1,1);
 elseif cor_dir > 0
   ord(1) = 1;
 end
end
for i = 2 : Np 
   whichI = find(parent(:,dimension,i)<fracFrac(i) & parent(:,dimension,i)>=fracFrac(i-1));
   whichIonsN(i,1) = length(whichI);
  for j = 1:length(whichI)
    whichIons(j,i) = whichI(j);
  end
  if ordering_on 
    if whichIonsN(i,1) > 0
      ord(i) = sum(POP_STRUC.POPULATION(parents(i)).order(whichI))/whichIonsN(i,1);
    elseif cor_dir > 0
      ord(i) = 1;
    end
  end
end

% chooses the more ordered slab out of Nslabs(i) for all parents, don't do it for ORG_STRUC.manyParents == 3
% it is easier to make a general Nslabs loop and only perform the slicing when Nslabs(j) <= i
if ordering_on & (ORG_STRUC.manyParents < 3)
 for i = 1 : max(Nslabs)     % ORG_STRUC.Nslabs
  offset_extra = rand(Np,1);
  for j = 1 : Np
   if Nslabs(j) <= i
     parent(:,dimension,j) = parent(:,dimension,j) + offset_extra(j);
     parent(:,dimension,j) = parent(:,dimension,j) - floor(parent(:,dimension,j));
   end
  end
  % for varcomp - fill the unused atoms with -1s
  for i2 = 1 : 3
   for i3 = 1 : Np
    for i1 = sum(POP_STRUC.POPULATION(parents(i3)).numIons)+1 : nIonsMax
     parent(i1,i2,i3) = -1;
    end
   end
  end
  ord_extra = zeros(Np,1);
  whichIons_extra = zeros(nIonsMax,Np);
  whichIonsN_extra = zeros(Np,1);
% ----------------- 1st parent
  if Nslabs(1) <= i
    whichI = find((parent(:,dimension,1)) < fracFrac(1) & (parent(:,dimension,1) > -0.5));
    whichIonsN_extra(1,1) = length(whichI);
    for j = 1 : length(whichI)
      whichIons_extra(j,1) = whichI(j);
    end
    if whichIonsN_extra(1,1) > 0
      ord_extra(1) = sum(POP_STRUC.POPULATION(parents(1)).order(whichI))/whichIonsN_extra(1,1);
    elseif cor_dir > 0
      ord_extra(1) = 1;
    end
  end
% ------------------------- other parents
  for i1 = 2 : Np 
   if Nslabs(i1) <= i
     whichI = find(parent(:,dimension,i1)<fracFrac(i1) & parent(:,dimension,i1)>=fracFrac(i1-1));
     whichIonsN_extra(i1,1) = length(whichI);
     for j = 1 : length(whichI)
       whichIons_extra(j,i1) = whichI(j);
     end
     if whichIonsN_extra(i1,1) > 0
       ord_extra(i1) = sum(POP_STRUC.POPULATION(parents(i1)).order(whichI))/whichIonsN_extra(i1,1);
     elseif cor_dir > 0
       ord_extra(i1) = 1;
     end
   end
  end
% ----------------------------------------
  for j = 1 : Np
   if Nslabs(j) <= i
    if ((ord(j) > ord_extra(j)) & (cor_dir <= 0)) | ((ord(j) < ord_extra(j)) & (cor_dir == 1))
     parent(:,dimension,j) = parent(:,dimension,j) - offset_extra(j);
     parent(:,dimension,j) = parent(:,dimension,j) - floor(parent(:,dimension,j));
     % for varcomp - fill the unused atoms with -1s
     for i2 = 1 : 3
      for i3 = 1 : Np
       for i1 = sum(POP_STRUC.POPULATION(parents(i3)).numIons)+1 : nIonsMax
        parent(i1,i2,i3) = -1;
       end
      end
     end
    else
     whichIons(j,:) = whichIons_extra(j,:);
     whichIonsN(j) = whichIonsN_extra(j);
     ord(j) = ord_extra(j);
    end
   end
  end
 end
end

% prepare the variable offspring. This needs to be done due to the command
% 'cat' used later on
offspring = zeros(0,3);


ionCount = zeros(length(ORG_STRUC.atomType),1);
ionCh = zeros(Np, length(ORG_STRUC.atomType)+1);

for i = 1 : Np
  ionCh(i, :) = ORG_STRUC.ionCh;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we have to fix the correct number of atoms/blocks in the offspring %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 desired_comp = ORG_STRUC.numIons;


for ind = 1:length(ORG_STRUC.atomType)
    % this loop is for all the different ion types.

    % find the indices of the ions for the type of ion done in this turn of the loop
    ions = zeros(nIonsMax,Np);
    smaller = zeros(nIonsMax,Np);
    ionsN = zeros(Np,1);
   for i = 1:Np 
     smaller1 = find(whichIons(1:whichIonsN(i,1),i) <= ionCh(i, ind+1));
     ion      = find(whichIons(smaller1,i) > ionCh(i,ind));
     ionsN(i,1) = length(ion);
     for j = 1:length(smaller1)
       smaller(j,i) = smaller1(j);
     end
     for j = 1:length(ion)
       ions(j,i) = ion(j);
     end
     for j = length(ion)+1:nIonsMax
       ions(j,i) = 0;
     end
   end
    % total number of this type of ion
    ionCount(ind,1) = sum(ionsN);
      
    candidates = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if the number of ions we have is smaller than the number we need,
    % then we have to fill up with further ions
    if ionCount(ind,1) < desired_comp(ind)

        % calculate how many ions (of this type) are missing
        howmany = desired_comp(ind) - ionCount(ind,1);

        % determine how many ions are to be filled into each fraction probabilistically.
        % Linear probability of each ion to be of a fraction (linear with the size of the fraction)
        % for each fraction we choose the determined amount of ions. The new ions are chosen 
        % (randomly) from the structure the ions in this fraction weren't originally chosen from

        filled = zeros(Np, desired_comp(ind));
        for xy = 1:howmany
          flag = 1;
counter = 0;
          while flag
counter = counter+1;
if counter > 2000
ins_par
tempInd
ins_ind
atoms
filled
fracFrac
coord1
quit
end
             ins_par = RandInt(1,1,[1,Np]);
             if ordering_on
               tempInd = (ionCh(ins_par, ind)+1:ionCh(ins_par, ind+1));
               atoms = sort_order(POP_STRUC.POPULATION(parents(ins_par)).order, tempInd, 'descend');
               ins_ind = tempInd(atoms(1))-ionCh(ins_par, ind);
               for i = 1 : length(tempInd)
                if filled(ins_par, tempInd(atoms(i))-ionCh(ins_par, ind)) == 0
                  ins_ind = tempInd(atoms(i))-ionCh(ins_par, ind);
                  break;
                end
               end
             else 
               ins_ind = RandInt(1,1,[1,desired_comp(ind)]);
             end
             if filled(ins_par, ins_ind) == 0
               flag = 0;
             end;
             coord1 = parent(ionCh(ins_par, ind)+ins_ind,dimension,ins_par);
             if ins_par == 1
               if coord1 < fracFrac(1)
                 flag = 1;
               end  
             elseif (coord1 >= fracFrac(ins_par-1)) & (coord1 < fracFrac(ins_par))
               flag = 1;
             end;
             if (ORG_STRUC.manyParents == 3) & (round(ins_par/2) == floor(ins_par/2)) % even ins_par
               for i = 1 : floor(Np/2)
                if (coord1 >= fracFrac(2*i-1)) & (coord1 < fracFrac(2*i))
                 flag = 1;
                end
               end
             end
             if (ORG_STRUC.manyParents == 3) & (round(ins_par/2) == floor(ins_par/2) + 1) % odd ins_par
               if coord1 < fracFrac(1)
                 flag = 1;
               end  
               fracFrac(end+1) = 1;
               for i = 1 : floor(Np/2)
                if (coord1 >= fracFrac(2*i)) & (coord1 < fracFrac(2*i+1))
                 flag = 1;
                end
               end
             end
             filled(ins_par, ins_ind) = 1;
          end
          filled(ins_par, ins_ind) = 1;
          candidates = cat(1,candidates,parent(ionCh(ins_par, ind)+ins_ind,:,ins_par));
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if on the other hand the number of ions we have is too big, we need to eliminate a few
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif ionCount(ind,1) > desired_comp(ind)

        howmany = ionCount(ind,1) - desired_comp(ind);
      
    % now we delete the surplus ions (randomly within each fraction) 
         ionsN1 = ionsN;

         if ordering_on
             order = zeros(1,ionCount(ind,1));
             ord_ind = 1;
             for i = 1 : Np
              for j = 1 : ionsN(i,1)
               order(ord_ind) = POP_STRUC.POPULATION(parents(i)).order(whichIons(smaller(ions(j,i),i),i));
               ord_ind = ord_ind + 1;
              end
             end
             tournament = zeros(ionCount(ind,1),1);
             % the following tournament is linear
             tournament(ionCount(ind,1)) = 1;
             for loop = 2:ionCount(ind,1)
                tournament(end-loop+1)= tournament(end-loop+2)+loop;
             end
             [junk, ranking] = sort(order); %junk = order(ranking)
             atom_rank = randperm(max(tournament));
             atoms = zeros(length(tournament),1);
             chosen = zeros(ionCount(ind,1),1);
             j1 = 1;
             for i1 = 1:max(tournament) 
               at1 = find (tournament > (atom_rank(i1)-1));
               if chosen(ranking(at1(end))) == 0
                  chosen(ranking(at1(end))) = 1;
                  atoms(j1) = ranking(at1(end));
                  j1 = j1 + 1;
               end
             end
         end

         for xy = 1:howmany
           flag = 1;
           while flag
            if ordering_on
              del_ind = atoms(end-xy+1); % we want worst to be choosen more often
            else
              del_ind = RandInt(1,1,[1,ionCount(ind,1)]);
            end
            parent_ind = 1;
            while del_ind > ionsN(parent_ind,1)
              del_ind = del_ind - ionsN(parent_ind,1);
              parent_ind = parent_ind + 1;
            end
            if ions(del_ind, parent_ind) > 0
             flag = 0;
            end  
           end

           ions(del_ind, parent_ind) = 0;
         end

    end

    addOn = candidates;
    % fill the selected ions coordinates into the offspring
    for in1 = 1:Np
      for in2 = 1:ionsN(in1)
       if ions(in2, in1) > 0
         addOn = cat(1,addOn, parent(whichIons(smaller(ions(in2,in1),in1),in1),:,in1));
       end; 
      end;
    end;

    offspring = cat(1,offspring,addOn);
end

  numIons = ORG_STRUC.numIons;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% % substract the offset (thus returning to the true values)
% offspring(:,dimension) = offspring(:,dimension) - offset;
% % adjust to the [0,1] boundaries
% offspring(:,dimension) = offspring(:,dimension)  - floor(offspring(:,dimension) );
%
% HOWEVER: substracting the offset is not necessary, since it changes
% nothing regarding the structure. If we leave the offset as it is, then we
% leave the substructures at their changed position. I asume this to be
% good against convergence.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% change back to a colomn vector
%quit
