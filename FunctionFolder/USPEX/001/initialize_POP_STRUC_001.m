function initialize_POP_STRUC_001()

% USPEX Version 9.4.0
% modulization

global ORG_STRUC
global POP_STRUC
global CLUSTERS

POP_STRUC = struct('POPULATION',{}, 'SOFTMODEParents',{}, 'SOFTMUTATED',{}, 'resFolder', {},'generation',{},...
 'DoneOrder',{}, 'bodyCount', {}, 'ranking',{},'bad_rank',{}, 'convex_hull',{}, 'fitness', {});

POP_STRUC(1).POPULATION = struct('COORDINATES', {}, 'INIT_COORD', {}, 'LATTICE', {}, 'INIT_LAT', {}, 'INIT_numIons', {}, 'numIons',{}, ...
'Vol',{},'struc_entr',{},'order',{},'FINGERPRINT',{},'K_POINTS',{},'Step',{}, 'Enthalpies',{}, 'Error',{},'Done',{},'ToDo',{},'Parents',{},...
'S_order',{}, 'howCome',{},'JobID',{},'Folder',{}, 'Number',{}, 'symg', {});
POP_STRUC.POPULATION(1) = QuickStart(POP_STRUC.POPULATION);

POP_STRUC(1).SOFTMUTATED = struct('FINGERPRINT',{}, 'mutatedAt', {}, 'fallBack', {});
POP_STRUC(1).SOFTMODEParents=struct('lattice',{},'coordinates',{},'fingerprint',{},'eignFre',{},'eignVec',{},'Softmode_Fre',{},'Softmode_num',{},'numIons',{});

POP_STRUC.generation = 1;
POP_STRUC.bodyCount = 0;
POP_STRUC.bad_rank = 0;

%-----001mode-----creating CLUSTERS structure-----------------------------------------------------------------------
if size(ORG_STRUC.numIons,2) == 1      % one-component clusters
        CLUSTERS.number_compositions = ORG_STRUC.numIons(2,1) - ORG_STRUC.numIons(1,1) + 1;
        for i = 1:CLUSTERS.number_compositions
            CLUSTERS.composition(i).bestEnthalpy = 10000;
            CLUSTERS.composition(i).numIons = ORG_STRUC.numIons(1,1) + i - 1;
        end

elseif size(ORG_STRUC.numIons,2) == 2  % two-component clusters
        for i = 1 : ORG_STRUC.numIons(2,1) - ORG_STRUC.numIons(1,1) + 3
	    for j = 1 : ORG_STRUC.numIons(2,2) - ORG_STRUC.numIons(1,2) + 3
            	CLUSTERS.composition(i,j).bestEnthalpy = 10000;
            	CLUSTERS.composition(i,j).numIons = [i + ORG_STRUC.numIons(1,1) - 2, j + ORG_STRUC.numIons(1,2) - 2];
	    end;
        end
end
safesave ('CLUSTERS.mat', CLUSTERS)
%-------------------------------------------------------------------------------------------------------------------

%create good initial population. Every individual fulfills hard constraints.
    goodPop = 1;
    newSym = 1;
    badSymmetry = 0;
    nsym = 0;
    sym_coef = 1;
    failedDist = 0;
    minDistMatrice = ORG_STRUC.minDistMatrice;
tic
    while goodPop < ORG_STRUC.initialPopSize + 0.5
        
        failedTime = toc;
        if (failedDist > 10000) || (failedTime > 300)
         if minDistMatrice(1,1) > 0.8*ORG_STRUC.minDistMatrice(1,1)
           if failedTime > 300
             disp('WARNING! Can not generate a structure after 5 minutes. The minimum distance threshold will be lowered by 0.1.');
           else
             disp('WARNING! Can not generate a structure after 10000 tries. The minimum distance threshold will be lowered by 0.1.');
           end
           disp('Please check your IonDistances parameter.');
           disp(' ');
           minDistMatrice = 0.9*minDistMatrice;
           failedDist = 0;
           tic
         else
             disp('Could not generate a structure after 30000 tries or 15 minutes.');
             disp('Please check the input files. The calculation has to stop.');
             disp('Remember they should be much smaller than the real interatomic distances');
             disp('but not too small for pseudopotential overlap errors to kill interatomic repulsion.');
             quit;
         end
        end

        errorS = 0;

%        startLat = rand(6,1);
%        startLat(4:6) = (pi/2);
%        check_startLat = latConverter(startLat);
%        volLat = det(check_startLat);
%        ratio = ORG_STRUC.latVolume/volLat;
%        startLat(1:3) = startLat(1:3)*(ratio)^(1/3);
%        lat = latConverter(startLat);

        %numIons = ORG_STRUC.numIons;

if size(ORG_STRUC.numIons,2) == 1      % one-component clusters
   numIons = ORG_STRUC.numIons(1,1) + round((ORG_STRUC.numIons(2,1) - ORG_STRUC.numIons(1,1))*(goodPop-1)/(ORG_STRUC.initialPopSize-1));
elseif size(ORG_STRUC.numIons,2) == 2  % two-component clusters
   numComp = (ORG_STRUC.numIons(2,1) - ORG_STRUC.numIons(1,1) + 1) * (ORG_STRUC.numIons(2,2) - ORG_STRUC.numIons(1,2) + 1);
   curComp = 1 + round((numComp-1)*(goodPop-1)/(ORG_STRUC.initialPopSize-1));
   numIons(1) = ORG_STRUC.numIons(1,1) + fix((curComp - 1) / (ORG_STRUC.numIons(2,2) - ORG_STRUC.numIons(1,2) + 1));
   numIons(2) = ORG_STRUC.numIons(1,2) - 1 + (curComp - (numIons(1) - ORG_STRUC.numIons(1,1)) *(ORG_STRUC.numIons(2,2) - ORG_STRUC.numIons(1,2) + 1));
end
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 	%calculation of lat - 001mode
        startLat = rand(6,1);
        startLat(4:6) = (pi/2);
        check_startLat = latConverter(startLat);
        volLat = det(check_startLat);

        lat1 = 0;
        for it = 1 : length(ORG_STRUC.atomType)
            lat1 = lat1 + numIons(it)*ORG_STRUC.latVolume(it);
        end

        ratio = lat1/volLat;
        startLat(1:3) = startLat(1:3)*(ratio)^(1/3);
        lat = latConverter(startLat);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        POP_STRUC.POPULATION(goodPop).numIons = numIons;
        sym_coef = ORG_STRUC.sym_coef;
        if badSymmetry > 1000    % 100k clusters were generated in 110.7561 sec, so about 1k in 1 sec
          if sqrt(lat(1,1)*lat(2,2)) < lat(3,3)
            coef = (0.65 + 0.4*rand);
          else
            coef = (0.95 + 0.4*rand);
          end
          lat(3,3) = lat(3,3)*coef;
          lat(1,1) = lat(1,1)/sqrt(coef);
          lat(2,2) = lat(2,2)/sqrt(coef);
        elseif badSymmetry > 150
           badSymmetry = 0;
           newSym = 1;      % change the symmetry if can't generate the cluster
        end
        badSymmetry = badSymmetry + 1;
        if newSym
         tmp = ceil(rand*size(ORG_STRUC.nsymN,1));
         nsym = ORG_STRUC.nsym(ORG_STRUC.nsymN(tmp,1):ORG_STRUC.nsymN(tmp,2));
         newSym = 0;
        end
        [candidate, lat, errorS] = symope_001(nsym, numIons, lat, minDistMatrice*ORG_STRUC.sym_coef);


        if errorS == 1
          goodBad = 0;
        else
          goodBad = distanceCheck(candidate, lat, numIons, minDistMatrice*sym_coef);
          if goodBad
             goodBad = checkConnectivity(candidate, lat, numIons);
          end
        end

        if goodBad
          POP_STRUC.POPULATION(goodPop).LATTICE = lat;
          POP_STRUC.POPULATION(goodPop).COORDINATES = candidate;
tic
          failedDist = 0;
          newSym = 1;

          disp(['Cluster ' num2str(goodPop) ' built with the symmetry ' nsym ' composition ' num2str(numIons)])
          [lat, candidate] = makeCluster(lat, candidate, ORG_STRUC.vacuumSize(1));
          POP_STRUC.POPULATION(goodPop).LATTICE = lat;
          POP_STRUC.POPULATION(goodPop).COORDINATES = candidate;
          POP_STRUC.POPULATION(goodPop).howCome = '  Random  ';

          goodPop = goodPop + 1;
        else
          failedDist = failedDist + 1;
        end
    end

%%%%%%%%%%%%%%%%%%%% SEEDING %%%%%%%%%%%%%%%%%
    pick_Seeds();
if ORG_STRUC.doFing
 pickAntiSeeds();
end
Start_POP_001();
