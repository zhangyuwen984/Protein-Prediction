function SoftModeMutation_300(OFF_ID)

% USPEX v10
global POP_STRUC
global ORG_STRUC
global OFF_STRUC

%------------------------------------------------------
%-------Step1: Set up the inputs
%------------------------------------------------------
[POP_ID, numIons, lat, coor, SOFT_ID] = ChooseParent();

if isempty(ORG_STRUC.maxAt)
   maxAt = sum(numIons);
else
   maxAt = ORG_STRUC.maxAt;
end
maxCellIncrease = 3;
[freq, eigvector, supercells_sizes] = calcSoftModes_varcomp(ORG_STRUC.NvalElectrons, ...
                    ORG_STRUC.valences, lat, coor, numIons, maxAt, maxCellIncrease);


%------------------------------------------------------
%------Step2: calculate which softmodes to start
%------------------------------------------------------
%maxCell increases up to 3
last = POP_STRUC.SOFTMODEParents(SOFT_ID).Softmode_num;
good_freq = ObtainFreq(lat, coor, numIons, freq, eigvector, supercells_sizes, last);

%------------------------------------------------------
%------Step3: try soft mutation until it succeeds.
%------------------------------------------------------
for f = good_freq : length(freq)
    good1 = 0;
    good2 = 0;
    POP_STRUC.SOFTMODEParents(SOFT_ID).Softmode_num = f;  %!!!!!IMPORTANT.....
    [eigV, coord0, lat0, numIons0] = calcEigenvectorK(eigvector(:,f), supercells_sizes(f,:), lat, coor, numIons);
    superCellSize = supercells_sizes(f,:);
    % opposite direction is degenerate using fingerprints criteria
    opposite_degenerate = checkDegenerate(coord0, numIons0, lat0, eigV);
    %---SoftMutation: direction 1  
    [good1, MUT_LAT, MUT_COORD, deviation] = Do_soft_mutation(coord0, numIons0, lat0, eigV(:,1));
    if good1
        SaveOFF(POP_ID, OFF_ID, MUT_COORD, MUT_LAT, numIons0, freq, f, ...
                MUT_COORD-coord0, deviation, superCellSize);
    end
    %---SoftMutation: direction 2
    if opposite_degenerate == 0
       [good2, MUT_LAT, MUT_COORD, deviation] = Do_soft_mutation(coord0, numIons0, lat0, -1*eigV(:,1));
       if good2
          if good1 == 1
             ToWrite = length(OFF_STRUC.POPULATION)+1;
          else
             ToWrite = OFF_ID;
          end
          SaveOFF(POP_ID, ToWrite, MUT_COORD, MUT_LAT, numIons0, freq, f, ...
                  MUT_COORD-coord0, deviation, superCellSize);
       end
    end
%---------------------------------------------------------------------------------
    if max([good1, good2])>0 %Succeed
       break
    elseif f==length(freq)
       disp('Out of soft modes.........OOPPPS, Softmutation');
       Mutation_300(OFF_ID);
       % In principle we need to do stochastic modes
       % but here we do ordinary coormutation.....
    end
end

%%%%%%%%%%%%%    
function [POP_ID, numIons, lat, coor, SOFT_ID] = ChooseParent()
global POP_STRUC
global ORG_STRUC
goodParent=0;
while ~goodParent
    toMutate = find(ORG_STRUC.tournament>RandInt(1,1,[0,max(ORG_STRUC.tournament)-1]));
    ind = POP_STRUC.ranking(toMutate(end));
    if (POP_STRUC.POPULATION(ind).Enthalpies(end) < 9999) & (sum(POP_STRUC.POPULATION(ind).numIons) > 1) ...
        & (CompositionCheck(POP_STRUC.POPULATION(ind).numIons/ORG_STRUC.numIons))
        goodParent = 1;
        POP_ID = ind;
        f1 = POP_STRUC.POPULATION(ind).FINGERPRINT;
    end
end

SOFT_ID = 0;
for j = 1 : length(POP_STRUC.SOFTMODEParents)
    dist_ij = cosineDistance(f1, POP_STRUC.SOFTMODEParents(j).fingerprint, ORG_STRUC.weight);
    if dist_ij < ORG_STRUC.toleranceFing
       disp(['The parent has already been softmutated.......']);
       SOFT_ID = j;
       numIons = POP_STRUC.SOFTMODEParents(j).numIons;
       lat     = POP_STRUC.SOFTMODEParents(j).lattice;
       coor    = POP_STRUC.SOFTMODEParents(j).coordinates;
       break;
    end
end

if SOFT_ID == 0
    numIons = POP_STRUC.POPULATION(POP_ID).numIons;
    lat     = POP_STRUC.POPULATION(POP_ID).LATTICE;
    coor    = POP_STRUC.POPULATION(POP_ID).COORDINATES;
    SOFT_ID = length(POP_STRUC.SOFTMODEParents) + 1;
    POP_STRUC.SOFTMODEParents(SOFT_ID).numIons = numIons;
    POP_STRUC.SOFTMODEParents(SOFT_ID).lattice = lat;
    POP_STRUC.SOFTMODEParents(SOFT_ID).coordinates = coor;
    POP_STRUC.SOFTMODEParents(SOFT_ID).Softmode_num = 0;  %1st count 
    POP_STRUC.SOFTMODEParents(SOFT_ID).fingerprint = f1;
end
%%%%%%%%%%%%%%
function SaveOFF(POP_ID, OFF_ID, MUT_COORD, MUT_LAT, numIons0, freq, f, vec, deviation, superCellSize)
%Save the necessary items if the mutation is accepted
global POP_STRUC
global OFF_STRUC

OFF_STRUC.POPULATION(OFF_ID).COORDINATES = MUT_COORD;
OFF_STRUC.POPULATION(OFF_ID).LATTICE     = MUT_LAT;
OFF_STRUC.POPULATION(OFF_ID).numIons     = numIons0;
OFF_STRUC.POPULATION(OFF_ID).howCome     = 'softmutate';

info_parents = struct('parent', {},'mut_degree', {}, ...
     'mut_vec',{},'mut_mode',{},'mut_fre',{},'enthalpy',{});
info_parents(1).parent  = num2str(POP_STRUC.POPULATION(POP_ID).Number);
N_atom = sum(POP_STRUC.POPULATION(POP_ID).numIons);
info_parents.enthalpy   = POP_STRUC.POPULATION(POP_ID).Enthalpies(end)/N_atom;
info_parents.mut_degree = deviation;
info_parents.mut_mode   = f;
info_parents.mut_vec    = vec;
info_parents.mut_fre    = freq(f);
OFF_STRUC.POPULATION(OFF_ID).Parents = info_parents;
disp(['Structure ' num2str(OFF_ID) ' softmutated in freq = ' num2str(freq(f),'%6.3f') ...
      ' (mode # ' num2str(f) ')']);
disp(['Mut_max   = ' num2str(max(deviation),'%6.3f') ', Mut_average = ' num2str(mean(deviation),'%6.3f')]);
if prod(superCellSize) > 1
   disp(['SuperCell = ' num2str(superCellSize) '   Parent: ' info_parents.parent ]);
end
disp(['  ']);

%%%%%%%%%%%%%%%
function goodfreq = ObtainFreq(lat, coor, numIons, freq, eigvector, supercells_sizes, lastfreq);
%%%%%%%%%%% Check the degeneracy of frequencies using different criteria; 
global POP_STRUC
global OFF_STRUC
global ORG_STRUC

weight        = ORG_STRUC.weight;
atomType      = ORG_STRUC.atomType;
toleranceFing = ORG_STRUC.toleranceFing;
if lastfreq > 0
   [eigV, coord0, lat0, numIons0]     = calcEigenvectorK(eigvector(:,lastfreq), supercells_sizes(lastfreq,:), lat, coor, numIons);
   [MUT_LAT0, MUT_COORD0, deviation0] = move_along_SoftMode_Mutation(coord0, numIons0, lat0, eigV, 1);
else %The structure itself --- to distinguish accoustic mode
    MUT_LAT0 = lat;
    MUT_COORD0 = coor;
    numIons0 = numIons;
end
[Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(MUT_LAT0, MUT_COORD0, numIons0, atomType);
[order1, f1, atom_fing1]           = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons0);
goodfreq = lastfreq;
for i = lastfreq+1 : length(freq)
   goodfreq = goodfreq + 1;
   if abs(freq(goodfreq)) > 5.0e-4   %we don't use zeros frequency
      %%fing for the next one
       [eigV, coord0, lat0, numIons0]      = calcEigenvectorK(eigvector(:,i), supercells_sizes(i,:), lat, coor, numIons);
       [MUT_LAT0, MUT_COORD0, deviation0]  = move_along_SoftMode_Mutation(coord0, numIons0, lat0, eigV, 1);
       [Ni, V, dist_matrix, typ_i, typ_j]  = makeMatrices(MUT_LAT0, MUT_COORD0, numIons0, atomType);
       [order1, f2, atom_fing1]            = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons0);
       dist = cosineDistance(f1, f2, weight);
   else
       dist = 0;
   end

   if (dist > 2*toleranceFing)
      break
   else
      %disp('Find a degenerate mode, skip it.....')
   end
end

%%%%%%%%%%%%%%%%%%%%
function degenerate = checkDegenerate(coord0, numIons0, lat0, eigV)
global ORG_STRUC
weight        = ORG_STRUC.weight;
atomType      = ORG_STRUC.atomType;
tolerance     = ORG_STRUC.toleranceFing;

[MUT_LAT0, MUT_COORD0, deviation0] = move_along_SoftMode_Mutation(coord0, numIons0, lat0, eigV, 1);
[Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(MUT_LAT0, MUT_COORD0, numIons0, atomType);
[order1, f1, atom_fing1]           = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons0);

[MUT_LAT0, MUT_COORD0, deviation0] = move_along_SoftMode_Mutation(coord0, numIons0, lat0, -1*eigV, 1);
[Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(MUT_LAT0, MUT_COORD0, numIons0, atomType);
[order1, f2, atom_fing1]           = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons0);

dist = cosineDistance(f1, f2, weight);
if dist < tolerance
   degenerate = 1;
   %disp('Opposite direction leads to same structure, skip it.....')
else
   degenerate = 0;
end

%%%%%%%%%%%%%%%%%%%%%%
function  [goodAtomMutant, MUT_LAT, MUT_COORD, deviation] = Do_soft_mutation(coord0, numIons0, lat0, eigV)
global ORG_STRUC
for i = 0 : 10
   [MUT_LAT, MUT_COORD, deviation] = move_along_SoftMode_Mutation(coord0, numIons0, lat0, eigV(:,1), 1-i/13);
   [MUT_LAT, strainMatrix]= lattice_Mutation(MUT_LAT, 0.15);  %lattice mutation
   goodAtomMutant = distanceCheck(MUT_COORD, MUT_LAT, numIons0, ORG_STRUC.minDistMatrice);
   if goodAtomMutant
      break
   end
end

