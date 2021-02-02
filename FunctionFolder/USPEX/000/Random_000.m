function Random_000(Ind_No)

% implemented - USPEX Version 8.5.0
global ORG_STRUC
global OFF_STRUC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CREATING random structures using space groups provided %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

goodStructure = 0;

newSym = 1;
badSymmetry = 0;
failedDist = 0;
minDistMatrice = ORG_STRUC.minDistMatrice;
badLat = 0;
badStruc = 0;
tic

while goodStructure ~= 1

  failedTime = toc;
  if (failedDist > 10000) || (failedTime > 300)
   if minDistMatrice(1,1) > 0.8*ORG_STRUC.minDistMatrice(1,1)
     if failedTime > 300
       %disp('WARNING! Can not generate a structure after 5 minutes. The minimum distance threshold will be lowered by 0.1.');
       USPEXmessage(505,'',0);
     else
       %disp('WARNING! Can not generate a structure after 10000 tries. The minimum distance threshold will be lowered by 0.1.');
       USPEXmessage(506,'',0);
     end
     disp('Please check your IonDistances parameter.');
     disp(' ');
     minDistMatrice = 0.9*minDistMatrice;
     failedDist = 0;
     tic
   else
      disp('Could not generate a structure after 30000 tries or 15 minutes.');
      disp('Please check the input files. The calculation has to stop.');
      disp('Possible reasons:  unreasonably big IonDistances.');
      quit;
   end
  end

  errorS = 0;
  startLat = rand(6,1);
  startLat(4:6) = (pi/2);
  check_startLat = latConverter(startLat);
  volLat = det(check_startLat);
  ratio = ORG_STRUC.latVolume/volLat;
  startLat(1:3) = startLat(1:3)*(ratio)^(1/3); 
  lat = latConverter(startLat);
    
  numIons = ORG_STRUC.numIons;

  sym_coef = ORG_STRUC.sym_coef;
  if badSymmetry > 500    % 100k clusters were generated in 110.7561 sec, so about 1k in 1 sec
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
  [candidate, lat, errorS] = symope_000(nsym, numIons, lat, minDistMatrice*sym_coef);

  if errorS == 0
    goodStructure = distanceCheck(candidate, lat, numIons, minDistMatrice*sym_coef);
    if goodStructure == 1
       goodStructure = checkConnectivity(candidate, lat, numIons);  %we also need to check the connectivity
    end
    badStruc = badStruc + 1 - goodStructure;
    if badStruc > 1000
      %disp('Structure failed 1000 times during random structure generation.')
      USPEXmessage(504,'',0);
      disp(' ')
      badStruc = 0;
    end
  else
    goodStructure = 0;
  end

  if goodStructure == 1
      OFF_STRUC.POPULATION(Ind_No).COORDINATES = candidate;
      OFF_STRUC.POPULATION(Ind_No).LATTICE = lat;
      OFF_STRUC.POPULATION(Ind_No).Parents = [];
      OFF_STRUC.POPULATION(Ind_No).numIons = numIons;
      OFF_STRUC.POPULATION(Ind_No).howCome = '  Random  ';

      [lat, candidate] = makeCluster(lat, candidate, ORG_STRUC.vacuumSize(1));
      OFF_STRUC.POPULATION(Ind_No).COORDINATES = candidate;
      OFF_STRUC.POPULATION(Ind_No).LATTICE = lat;
      disp(['Structure ' num2str(Ind_No) ' generated with random symmetry: '  num2str(nsym)]);
  else
    failedDist = failedDist + 1;
  end

end
