function [candidate_2D, lat_2D, candidate, lat] = Random_Init_M210(Ind_No, numMols)

% implemented - USPEX Version 8.5.0
% This is only for the example of Boron
% Needs to be improved later
global ORG_STRUC

newSym = 1;
badSymmetry = 0;
Molecules = struct('MOLCOORS',{},'ZMATRIX',{},'ID',{},'MOLCENTER',{});

CenterMat = ORG_STRUC.CenterminDistMatrice-1;
minDistMatrice = ORG_STRUC.minDistMatrice;
thickness2 = ORG_STRUC.thicknessS;
thickness1 = ORG_STRUC.thicknessB;
goodBad = 0;

while ~goodBad
   [typesAList, MtypeLIST, numIons] = GetPOP_MOL(numMols);

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
     [coor, lat_2D, numSites, Operation, errorS] = symope_2D_MOL(nsym, numMols, ...
                   1.2*ORG_STRUC.latVolume/(thickness2), CenterMat+1, thickness1);
     cd(ORG_STRUC.homePath)

     if errorS == 0
           [lat_2D,coor] = make2D(lat_2D, coor, thickness2);              
           for item=1:30
               Molecules = [];
               [typesAList, MtypeLIST, numIons] = GetPOP_MOL(numMols);
               P = [1 2 3]; PB = [1,2,3];
               Molecules = GetOrientation(coor, lat_2D, numSites(:,1), Operation, MtypeLIST, nsym, P, PB);
               goodBad = newMolCheck(Molecules,lat_2D, MtypeLIST, ORG_STRUC.minDistMatrice+0.4);  %change distance
               if goodBad
                   disp(['Structure ' num2str(Ind_No) '  generated from Random symmetry']);
                   cattedCoors = [];
                   for j = 1: sum(numMols)
                      cattedCoors = cat(1,cattedCoors,Molecules(j).MOLCOORS);
                   end
                   saveded = cattedCoors/lat_2D;
                   newCoords = zeros(0,3);
                   for m = 1: length(ORG_STRUC.atomType)
                     s = find(typesAList == ORG_STRUC.atomType(m));
                     newCoords = cat(1,newCoords, saveded(s,:));
                   end

                   candidate_2D = newCoords;
                   [lat,candidate] = make2D(lat_2D, candidate_2D, ORG_STRUC.vacuumSize(1));              
                   break
               end
           end
     end %if errorS
end %end while
