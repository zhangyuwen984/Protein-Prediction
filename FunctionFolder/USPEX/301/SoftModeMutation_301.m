function SoftModeMutation_301(Ind_No)

% USPEX Version 9.4.0
% this is a mutation along the soft mode
% Last updated by Qiang Zhu (2013/11/26)

global POOL_STRUC
global POP_STRUC
global ORG_STRUC
global OFF_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  CREATING Mutants by atom positions mutation  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodAtomMutant = 0;
goodMutLattice = 0;
%goodComposition = 0;

structureFailed = 0;

while goodAtomMutant + goodMutLattice ~= 2
   tryTime = 0;
   
   ind = chooseGoodComposition(ORG_STRUC.tournament, POOL_STRUC.POPULATION);
   if (tryTime>200) | (ind<0)
          USPEXmessage(516,'',0);
		  Random_301(Ind_No);
          return;
    end
        
   %while ~goodComposition
   %    toMutate = find(ORG_STRUC.tournament>RandInt(1,1,[0,max(ORG_STRUC.tournament)-1]));
   %    ind = toMutate(end);
   %    numBlocks = POOL_STRUC.POPULATION(ind).numBlocks;
   %    goodComposition = CompositionCheck(numBlocks);
% 	   tryTime = tryTime + 1;
% 	   if tryTime > 200
%           USPEXmessage(516,'',0);
% 		  Random_301(Ind_No);
%           return;
% 	   end
%    end

   numIons = POOL_STRUC.POPULATION(ind).numIons;
    N = sum(numIons); % N is used for Grey Fp calculation
   lat = POOL_STRUC.POPULATION(ind).LATTICE;
   coor = POOL_STRUC.POPULATION(ind).COORDINATES;
   vol = det(lat);
   flag = 0;
   maxSoftMode = N;  % after this soft mode number is reached, we start stochastic softmutation
  for j = 1 : length(POP_STRUC.SOFTMODEParents)
   if (sum(abs(numIons - POP_STRUC.SOFTMODEParents(j).numIons)) > 0.5)  %different composition
     continue
   end
%%%%%%%%%%%%%% If the chosen struc has been used to do mutation, move to the next soft mode
     vol1 = det(POP_STRUC.SOFTMODEParents(j).lattice);
  if (vol-vol1)/vol < 0.03
    [freq, eigvector] = calcSoftModes(ORG_STRUC.NvalElectrons, ORG_STRUC.valences, numIons, lat, coor);
    freq = diag(freq);
    [freq, IX] = sort(freq);

    L_char = (det(lat)/N)^(1/3); % characteristic length - edge of a cube of volume for one atom
    ff = 0;
    non_zero = 0;
    last_good = POP_STRUC.SOFTMODEParents(j).Softmode_Fre;
    maxSoftModeFreq = freq(maxSoftMode);% after this soft mode frequency is reached, we start stochastic softmutation

%%%%%%%%%%% Now we will check the degeneracy of frequencies using different criteria; fingerprint criterion used only for non-stochastic softmutation
%%%%%%%%%%% Note that fingerprint distance is shorter, to speed the things up
      if ORG_STRUC.doFing
       if (last_good <= maxSoftModeFreq)
         [MUT_LAT0, MUT_COORD0, deviation0] = move_along_SoftMode_Mutation(coor, numIons, lat, eigvector(:,IX(POP_STRUC.SOFTMODEParents(j).Softmode_num)), 1);
          [Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(MUT_LAT0, MUT_COORD0, N, 1);
         [order0, fingerprint0, atom_fing0] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, N);
       end
      end


    for i = 1 : length(freq)
      if freq(i) < last_good
        continue;
      end

      quite_different = 0;
      if ORG_STRUC.doFing
       if (freq(i) <= maxSoftModeFreq)
         [MUT_LAT0, MUT_COORD0, deviation0] = move_along_SoftMode_Mutation(coor, numIons, lat, eigvector(:,IX(i)), 1);
           [Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(MUT_LAT0, MUT_COORD0, N, 1);
           [order1, fingerprint1, atom_fing1] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, N);
         if (cosineDistance(fingerprint0, fingerprint1, ORG_STRUC.weight) > ORG_STRUC.toleranceFing)
           quite_different = 1;
         end
       else % stochastic softmutaiton regime
        non_zero = maxSoftMode + 1;
        break;
       end
      elseif freq(i) > (1.05)*last_good
         quite_different = 1;     
      end

      if quite_different == 1
         non_zero = i;
         break
      end
    end
          
    if non_zero == 0  % (the structure can not produce softmutation any more)
       Mutation_301(Ind_No);
       logic1 = 1;
       logic2 = 1;
    else
       opposite_degenerate = 0;
       current_freq = non_zero;
       good_freq = non_zero;
       loop = 0;
       logic1 = 0;
       logic2 = 0;

       for f = good_freq : non_zero + round((3*N-non_zero)/1.5)  % about middle of the freq spectra if we ignore the degeneracy
         if loop==1
           break
         else
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SoftMutation: direction 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          for i = 0 : 10
           if f > maxSoftMode   % stochastic softmutation
             stoch_vect = (2*rand-1)*eigvector(:,IX(1));
             for st = 2 : maxSoftMode
               stoch_vect = stoch_vect + (2*rand-1)*eigvector(:,IX(st));
             end
             [MUT_LAT, MUT_COORD, deviation] = move_along_SoftMode_Mutation(coor, numIons, lat, stoch_vect, 1-i/21);
           else
             [MUT_LAT, MUT_COORD, deviation] = move_along_SoftMode_Mutation(coor, numIons, lat, eigvector(:,IX(f)), 1-i/21);
           end
           goodAtomMutant = distanceCheck(MUT_COORD, MUT_LAT, numIons, ORG_STRUC.minDistMatrice);
           structureFailed = structureFailed + 1 - goodAtomMutant;
           if structureFailed >= 100
             %disp('Distance check failed more than 100 times in a single softmutation.')
             USPEXmessage(513,'',0); 
             structureFailed = 0;
           end
           goodMutLattice = 1; % since we don't change it
           if goodAtomMutant + goodMutLattice == 2

% check the degeneracy of moving in opposite direction, using fingerprints
            if ORG_STRUC.doFing & (f <= maxSoftMode) 
               [Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(MUT_LAT, MUT_COORD, N, 1);
               [order0, fingerprint0, atom_fing0] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, N);
               [MUT_LAT1, MUT_COORD1, deviation1] = move_along_SoftMode_Mutation(coor, numIons, lat, -1*eigvector(:,IX(f)), 1-i/21);
               [Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(MUT_LAT1, MUT_COORD1, N, 1);
               [order1, fingerprint1, atom_fing1] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, N);
             if (cosineDistance(fingerprint0, fingerprint1, ORG_STRUC.weight) <= ORG_STRUC.toleranceFing)
               opposite_degenerate = 1;
             end              
            end

              logic1 = 1;  
              OFF_STRUC.POPULATION(Ind_No).COORDINATES = MUT_COORD;
              OFF_STRUC.POPULATION(Ind_No).LATTICE = MUT_LAT;
              OFF_STRUC.POPULATION(Ind_No).howCome = 'softmutate';
              info_parents = struct('parent', {},'mut_degree', {},'mut_mode',{},'mut_fre',{},'enthalpy',{});
              info_parents(1).parent = num2str(POOL_STRUC.POPULATION(ind).Number);
              info_parents.mut_degree = deviation;
              info_parents.mut_mode = f;
              info_parents.mut_fre = freq(f);
              info_parents.enthalpy = POOL_STRUC.POPULATION(ind).enthalpy;
              disp(['Structure ' num2str(Ind_No) ' generated by softmutation of a structure which has been selected']);

              OFF_STRUC.POPULATION(Ind_No).Parents = info_parents;
              OFF_STRUC.POPULATION(Ind_No).numIons = numIons;
              OFF_STRUC.POPULATION(Ind_No).numBlocks = POP_STRUC.SOFTMODEParents(j).numBlocks;
              POP_STRUC.SOFTMODEParents(j).Softmode_Fre = freq(f);
              POP_STRUC.SOFTMODEParents(j).Softmode_num = f;
              loop = 1;
              break;
           elseif (i == 21) & (f == non_zero + round((3*N-non_zero)/2))
              break;         
           end
          end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SoftMutation: opposite direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          for i = 0 : 10
           if (f > maxSoftMode) & (logic1 == 1)  % stochastic softmutation was successfully done
             break;
           elseif opposite_degenerate == 1 % opposite direction is degenerate using fingerprints criteria
             break;
           elseif f > maxSoftMode    % stochastis softmutation failed in one direction - thus try it once more for the sake of it
             stoch_vect = (2*rand-1)*eigvector(:,IX(1));
             for st = 2 : maxSoftMode
               stoch_vect = stoch_vect + (2*rand-1)*eigvector(:,IX(st));
             end
             [MUT_LAT, MUT_COORD, deviation] = move_along_SoftMode_Mutation(coor, numIons, lat, stoch_vect, 1-i/21);
           else
             [MUT_LAT, MUT_COORD, deviation] = move_along_SoftMode_Mutation(coor, numIons, lat,-1*eigvector(:,IX(f)), 1-i/21);
           end
           goodAtomMutant = distanceCheck(MUT_COORD, MUT_LAT, numIons, ORG_STRUC.minDistMatrice);
           structureFailed = structureFailed + 1 - goodAtomMutant;
           if structureFailed >= 100
             %disp('Distance check failed more than 100 times in a single softmutation.')
             USPEXmessage(513,'',0);
             structureFailed = 0;
           end
           goodMutLattice = 1; % since we don't change it
           if goodAtomMutant + goodMutLattice == 2
              logic2 = 1;
              info_parents = struct('parent', {},'mut_degree', {},'mut_mode',{},'mut_fre',{}, 'enthalpy',{});
              info_parents(1).parent = num2str(POOL_STRUC.POPULATION(ind).Number);
              info_parents.mut_degree = deviation;
              info_parents.mut_mode = f;
              info_parents.mut_fre = freq(f);
              info_parents.enthalpy = POOL_STRUC.POPULATION(ind).enthalpy;
              if logic1 == 1
                OFF_STRUC.POPULATION(end+1).COORDINATES = MUT_COORD;
                OFF_STRUC.POPULATION(end).LATTICE = MUT_LAT;
                OFF_STRUC.POPULATION(end).Parents = info_parents;
                OFF_STRUC.POPULATION(end).numIons = numIons;
                OFF_STRUC.POPULATION(end).numBlocks = POP_STRUC.SOFTMODEParents(j).numBlocks;
                OFF_STRUC.POPULATION(end).howCome = 'softmutate';
              else
                OFF_STRUC.POPULATION(Ind_No).COORDINATES = MUT_COORD;
                OFF_STRUC.POPULATION(Ind_No).LATTICE = MUT_LAT;
                OFF_STRUC.POPULATION(Ind_No).Parents = info_parents;
                OFF_STRUC.POPULATION(Ind_No).numIons = numIons;
                OFF_STRUC.POPULATION(Ind_No).numBlocks = POP_STRUC.SOFTMODEParents(j).numBlocks;
                OFF_STRUC.POPULATION(Ind_No).howCome = 'softmutate';
              end
              
              POP_STRUC.SOFTMODEParents(j).Softmode_Fre = freq(f);
              POP_STRUC.SOFTMODEParents(j).Softmode_num = f;
              disp(['Structure ' num2str(Ind_No) ' generated by softmutation of a structure which has been selected']);
              loop = 1;
              break;
           elseif (i == 21) & (f == non_zero + round((3*N-non_zero)/2))
              break;
           end
          end

         end  %if loop==1
         good_freq = good_freq + 1;
       end  %% Outer loop, freq++ 
    end     %% Loop: nonzero
    flag = 1;
    break;    %since we already run out of all the choices, exit 
  end        % Loop: cosine_dist
 end          % Loop: comparison with parent database

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%% if the chosen struc is a new struc %%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
 if flag==0
   [freq, eigvector] = calcSoftModes(ORG_STRUC.NvalElectrons, ORG_STRUC.valences, numIons, lat, coor);
   freq = diag(freq);
   [freq, IX] = sort(freq);
   for i = 1:length(freq)        % first non zero element, should usually be 4       
     if freq(i) > 0.0000001
       non_zero = i;
       break;
     end
   end
   non_zero = 4; % try to remove only 3 main accoustic modes - experimental
   current_freq = non_zero;
   good_freq = non_zero;
  
   L_char = (det(lat)/N)^(1/3); % characteristic length - edge of a cube of volume for one atom
   ff = 0;
   loop = 0; 
   logic1 = 0;
   logic2 = 0;
   opposite_degenerate = 0;

   for f = good_freq : non_zero + round((3*N-non_zero)/1.5)  % about middle of the freq spectra if we ignore the degeneracy
     if loop==1
       break
     else
       for i = 0 : 10
         [MUT_LAT, MUT_COORD, deviation] = move_along_SoftMode_Mutation(coor, numIons, lat, eigvector(:,IX(f)), 1-i/21);
         goodAtomMutant = distanceCheck(MUT_COORD, MUT_LAT, numIons, ORG_STRUC.minDistMatrice);
         structureFailed = structureFailed + 1 - goodAtomMutant;
         if structureFailed >= 100
           %disp('Distance check failed more than 100 times in a single softmutation.')
           USPEXmessage(513,'',0);
           structureFailed = 0;
         end

         goodMutLattice = 1; % since we don't change it

         if goodAtomMutant + goodMutLattice == 2

% check the degeneracy of moving in opposite direction, using fingerprints
            if ORG_STRUC.doFing & (f <= maxSoftMode) 
               [Ni, V,  dist_matrix, typ_i, typ_j] = makeMatrices(MUT_LAT, MUT_COORD, N, 1);
               [order0, fingerprint0, atom_fing0] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, N);
               [MUT_LAT1, MUT_COORD1, deviation1] = move_along_SoftMode_Mutation(coor, numIons, lat, -1*eigvector(:,IX(f)), 1-i/21);
               [Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(MUT_LAT1, MUT_COORD1, N, 1);
               [order1, fingerprint1, atom_fing1] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, N);
             if (cosineDistance(fingerprint0, fingerprint1, ORG_STRUC.weight) <= ORG_STRUC.toleranceFing)
               opposite_degenerate = 1;
             end              
            end

             logic1 = 1;
             POP_STRUC.SOFTMODEParents(end+1).lattice = MUT_LAT;
             POP_STRUC.SOFTMODEParents(end).coordinates = coor;
             POP_STRUC.SOFTMODEParents(end).Softmode_Fre = freq(f);
             POP_STRUC.SOFTMODEParents(end).Softmode_num = f;
             POP_STRUC.SOFTMODEParents(end).numIons = numIons;
             POP_STRUC.SOFTMODEParents(end).numBlocks = POOL_STRUC.POPULATION(ind).numBlocks;
             OFF_STRUC.POPULATION(Ind_No).COORDINATES = MUT_COORD;
             OFF_STRUC.POPULATION(Ind_No).LATTICE = MUT_LAT;
             ff = f;
             info_parents = struct('parent', {},'mut_degree', {},'mut_mode',{},'mut_fre',{}, 'enthalpy', {});
             info_parents(1).parent = num2str(POOL_STRUC.POPULATION(ind).Number);
             info_parents.mut_degree = deviation;
             info_parents.mut_mode = f;
             info_parents.mut_fre = freq(f);
             info_parents.enthalpy = POOL_STRUC.POPULATION(ind).enthalpy;
             disp(['Structure ' num2str(Ind_No) ' generated by softmutation of a new structure']);
             OFF_STRUC.POPULATION(Ind_No).Parents = info_parents;
             OFF_STRUC.POPULATION(Ind_No).numIons = numIons;
             OFF_STRUC.POPULATION(Ind_No).numBlocks = POOL_STRUC.POPULATION(ind).numBlocks;
             OFF_STRUC.POPULATION(Ind_No).howCome = 'softmutate';
             loop = 1;
             break;
         elseif (i == 21) & (f == non_zero + round((3*N-non_zero)/2))
            break;
         end
       end

       if opposite_degenerate == 1 % opposite direction is degenerate using fingerprints criteria
          break;
       end

       for i = 0 : 10
         [MUT_LAT, MUT_COORD, deviation] = move_along_SoftMode_Mutation(coor, numIons, lat, -1*eigvector(:,IX(f)), 1-i/21);
         goodAtomMutant = distanceCheck(MUT_COORD, MUT_LAT, numIons, ORG_STRUC.minDistMatrice);
         structureFailed = structureFailed + 1 - goodAtomMutant;
         if structureFailed >= 100
           %disp('Distance check failed more than 100 times in a single softmutation.')
           USPEXmessage(513,'',0);
           structureFailed = 0;
         end
         goodMutLattice = 1; % since we don't change it
         if goodAtomMutant + goodMutLattice == 2
             logic2 = 1;
             ff = f;
             info_parents = struct('parent', {},'mut_degree', {},'mut_mode',{},'mut_fre',{}, 'enthalpy', {});
             info_parents(1).parent = num2str(POOL_STRUC.POPULATION(ind).Number);
             info_parents.mut_degree = deviation;
             info_parents.mut_mode = f;
             info_parents.mut_fre = freq(f);
             info_parents.enthalpy = POOL_STRUC.POPULATION(ind).enthalpy;
             disp(['Structure ' num2str(Ind_No) ' generated by softmutation of a new structure']);
             if logic1 == 1
              OFF_STRUC.POPULATION(end+1).COORDINATES = MUT_COORD;
              OFF_STRUC.POPULATION(end).LATTICE = MUT_LAT;
              OFF_STRUC.POPULATION(end).Parents = info_parents;
              OFF_STRUC.POPULATION(end).numIons = numIons;
              OFF_STRUC.POPULATION(end).numBlocks = POOL_STRUC.POPULATION(ind).numBlocks;
              OFF_STRUC.POPULATION(end).howCome = 'softmutate';
             else
              OFF_STRUC.POPULATION(Ind_No).COORDINATES = MUT_COORD;
              OFF_STRUC.POPULATION(Ind_No).LATTICE = MUT_LAT;
              OFF_STRUC.POPULATION(Ind_No).Parents = info_parents;
              OFF_STRUC.POPULATION(Ind_No).numIons = numIons;
              OFF_STRUC.POPULATION(Ind_No).numBlocks = POOL_STRUC.POPULATION(ind).numBlocks;
              OFF_STRUC.POPULATION(Ind_No).howCome = 'softmutate';
              POP_STRUC.SOFTMODEParents(end+1).lattice = MUT_LAT;
              POP_STRUC.SOFTMODEParents(end).coordinates = coor;
              POP_STRUC.SOFTMODEParents(end).Softmode_Fre = freq(f);
              POP_STRUC.SOFTMODEParents(end).Softmode_num = f;
              POP_STRUC.SOFTMODEParents(end).numIons = numIons;
              POP_STRUC.SOFTMODEParents(end).numBlocks = POOL_STRUC.POPULATION(ind).numBlocks;
             end
             loop = 1;
             break;
         elseif (i == 21) & (f == non_zero + round((3*N-non_zero)/2))
            break;
         end
       end
       good_freq = good_freq+1;
     end  
   end
 end    
 
 if (logic1 == 0) & (logic2 == 0)
    Mutation_301(Ind_No);
 end

 goodAtomMutant = 1;
 goodMutLattice = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   END creating mutants   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
