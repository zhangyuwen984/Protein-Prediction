function ANS = Random_200(Ind_No, Vacuum, cell)

global ORG_STRUC
global OFF_STRUC

minDistMatrice = ORG_STRUC.minDistMatrice;

%---------Step1: to find the cell size and lattice
bulk_lat    =ORG_STRUC.bulk_lat;
bulk_pos    =ORG_STRUC.bulk_pos;
bulk_atyp   =ORG_STRUC.bulk_atyp;
bulk_numIons=ORG_STRUC.bulk_ntyp;
if prod(cell)>1
   bigcell= cell;
   bigcell(3)=1;
  [bulk_lat, bulk_pos, bulk_atyp, bulk_numIons]=supercell(bulk_lat, bulk_pos, bulk_atyp, bulk_numIons, bigcell);
end
Lat = latConverter(bulk_lat);
Lat(3) = ORG_STRUC.thicknessS;
sur_lat = latConverter(Lat);
sur_numIons = ORG_STRUC.numIons*prod(cell);


%%%%Step 2: to get the atomic positions
goodBad = 0;
Groups = ORG_STRUC.nsym;
ini_lat = latConverter(Lat);

while ~goodBad 
     tmp = find(Groups > 0);
     nsym = tmp(ceil(rand*length(tmp)));
     
     cd([ORG_STRUC.homePath '/CalcFoldTemp'])
     [sur_candidate, sur_lat, errorS] = symope_GB(nsym, sur_numIons, ini_lat, minDistMatrice);
     cd(ORG_STRUC.homePath)
     
     if errorS == 0
        [lat, candidate, numIons, typesAList, chanAList] = ...
        makeSurface(sur_lat, sur_candidate, sur_numIons, bulk_lat, bulk_pos, bulk_atyp, bulk_numIons, Vacuum);
        goodBad = distanceCheck(candidate, lat, numIons, minDistMatrice);
     else
        Groups(nsym)=0;
     end
     if goodBad
        disp(['Structure ' num2str(Ind_No) ' generated randomly']);
        disp(['composition: ' num2str(sur_numIons) ';   cell size: ' num2str(cell(1)) '*' num2str(cell(2))]);
        ANS.candidate     = candidate;
        ANS.numIons       = numIons;
        ANS.lat           = lat;
        ANS.typesAList    = typesAList;
        ANS.chanAList     = chanAList;
        ANS.sur_lat       = sur_lat;
        ANS.sur_candidate = sur_candidate;
        ANS.sur_numIons   = sur_numIons;
        ANS.bulk_lat      = bulk_lat;
        ANS.bulk_pos      = bulk_pos;
        ANS.bulk_atyp     = bulk_atyp;
        ANS.bulk_numIons  = bulk_numIons;
     
        break
     end
end
