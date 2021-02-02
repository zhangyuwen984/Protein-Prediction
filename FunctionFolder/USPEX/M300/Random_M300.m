function Random_M300(Ind_No)

% implemented - USPEX Version 8.5.0

global ORG_STRUC
global OFF_STRUC

bulk_lat    =ORG_STRUC.bulk_lat;
bulk_pos    =ORG_STRUC.bulk_pos;
bulk_atyp   =ORG_STRUC.bulk_atyp';
bulk_numIons=ORG_STRUC.bulk_ntyp;
Lat         = latConverter(bulk_lat);
Lat(3)      = ORG_STRUC.thicknessS;
GB_numIons=ORG_STRUC.numIons;
item = 0;
for i = 1:length(GB_numIons)
    GB_atyp(1,item+1:item+GB_numIons(i))=ORG_STRUC.atomType(i);
    item = item + GB_numIons(i);
end
numIons = GB_numIons + bulk_numIons;

goodBad = 0;
while ~goodBad
       tmp = find(ORG_STRUC.nsym > 0);
       nsym = tmp(ceil(rand*length(tmp)));
       lat1 = latConverter(Lat);
      if sum(ORG_STRUC.splitInto)<2
         cd([ORG_STRUC.homePath '/CalcFoldTemp']);
         [GB_candidate, GB_lat, errorS] = symope_GB(nsym, GB_numIons, lat1, ORG_STRUC.minDistMatrice);
         cd(ORG_STRUC.homePath)
      else
         splitInto = ORG_STRUC.splitInto(ceil(length(ORG_STRUC.splitInto)*rand));
         [GB_candidate, GB_lat, split, errorS] = splitBigCell_M300(lat1, splitInto, GB_numIons, nsym, ORG_STRUC.minDistMatrice);
      end

     if errorS == 0
        [lat,candidate, typesAList, chanAList] = makeGB(numIons, GB_lat, GB_candidate, GB_atyp, bulk_lat, bulk_pos, bulk_atyp, ORG_STRUC.vacuumSize(1));
        goodBad = distanceCheck(candidate, lat, numIons, ORG_STRUC.minDistMatrice);
     else
        goodBad =0;
     end

  if goodBad
      OFF_STRUC.POPULATION(Ind_No).Bulk_LATTICE    =bulk_lat;
      OFF_STRUC.POPULATION(Ind_No).Bulk_COORDINATES=bulk_pos;
      OFF_STRUC.POPULATION(Ind_No).Bulk_typesAList =bulk_atyp;
      OFF_STRUC.POPULATION(Ind_No).Bulk_numIons    =bulk_numIons;
      OFF_STRUC.POPULATION(Ind_No).COORDINATES = candidate;
      OFF_STRUC.POPULATION(Ind_No).numIons = numIons;
      OFF_STRUC.POPULATION(Ind_No).LATTICE = lat;
      OFF_STRUC.POPULATION(Ind_No).typesAList = typesAList;
      OFF_STRUC.POPULATION(Ind_No).GB_LATTICE = GB_lat;
      OFF_STRUC.POPULATION(Ind_No).GB_COORDINATES = GB_candidate;
      OFF_STRUC.POPULATION(Ind_No).GB_numIons = GB_numIons;
      OFF_STRUC.POPULATION(Ind_No).GB_typesAList = GB_atyp;
      OFF_STRUC.POPULATION(Ind_No).chanAList = chanAList;
      OFF_STRUC.POPULATION(Ind_No).howCome = ' Random ';
      if sum(ORG_STRUC.splitInto)>=2
         disp(['split into: ' num2str(split) '  subcells']);
      end

      disp(['Structure ' num2str(Ind_No) ' with plane group ' num2str(nsym) ' generated randomly']);
  end
end
