function Random_M200(Ind_No)

% implemented - USPEX Version 8.5.0

global ORG_STRUC
global OFF_STRUC

if ~isempty(ORG_STRUC.numMols)
   [candidate_2D, lat_2D, candidate, lat] = Random_Init_M210(Ind_No, ORG_STRUC.numMols);
else
   [candidate_2D, lat_2D, candidate, lat] = Random_Init_M200(Ind_No, ORG_STRUC.numIons);
end

OFF_STRUC.POPULATION(Ind_No).LATTICE = lat;
OFF_STRUC.POPULATION(Ind_No).LATTICE_2D = lat_2D;
OFF_STRUC.POPULATION(Ind_No).COORDINATES = candidate;
OFF_STRUC.POPULATION(Ind_No).COORDINATES_2D = candidate_2D;
OFF_STRUC.POPULATION(Ind_No).howCome = '  Random  ';
OFF_STRUC.POPULATION(Ind_No).numIons = ORG_STRUC.numIons;

