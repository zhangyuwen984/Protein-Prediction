function [candidate_2D, lat_2D, candidate, lat] = Random_Init_M200(Ind_No, numIons)

% implemented - USPEX Version 8.5.0
global ORG_STRUC

newSym = 1;
badSymmetry = 0;
goodBad = 0;
minDistMatrice = ORG_STRUC.minDistMatrice;
thickness      = ORG_STRUC.thicknessS;

while ~goodBad
    if (badSymmetry > 15) 
      badSymmetry = 0;
      newSym = 1;      % change the symmetry group if can't generate the crystal
    end

      badSymmetry = badSymmetry + 1;

    if newSym
      tmp = find(ORG_STRUC.nsym > 0);
      nsym = tmp(ceil(rand*length(tmp))); % pick a random group from those specified by user
      newSym = 0;
    end
    
     cd([ORG_STRUC.homePath '/CalcFoldTemp'])
     [candidate_2D, lat_2D, errorS] = symope_2D(nsym, numIons, ORG_STRUC.latVolume/(thickness+2.0), minDistMatrice, thickness);
     cd(ORG_STRUC.homePath)
   
    if errorS == 0
        [lat,candidate] = make2D(lat_2D, candidate_2D, ORG_STRUC.vacuumSize(1));
        goodBad = distanceCheck(candidate, lat, numIons, ORG_STRUC.minDistMatrice);
        if goodBad == 1
           goodBad = checkConnectivity(candidate, lat, numIons);
        end
    end
  
    if goodBad
        disp(['Structure ' num2str(Ind_No) ' generated with plane group: '  num2str(nsym)]);
    end
end
