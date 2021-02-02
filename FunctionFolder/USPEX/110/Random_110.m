function Random_110(Ind_No)

% implemented - USPEX Version 9.3.7
global ORG_STRUC
global OFF_STRUC

[Molecules, lat, MtypeLIST, typesAList, numIons] = Random_Init_110(Ind_No, ORG_STRUC.numMols);

OFF_STRUC.POPULATION(Ind_No).LATTICE = lat;
OFF_STRUC.POPULATION(Ind_No).MOLECULES = Molecules;
OFF_STRUC.POPULATION(Ind_No).numIons = numIons;
OFF_STRUC.POPULATION(Ind_No).numMols = ORG_STRUC.numMols;
OFF_STRUC.POPULATION(Ind_No).MtypeLIST = MtypeLIST;
OFF_STRUC.POPULATION(Ind_No).typesAList = typesAList;
OFF_STRUC.POPULATION(Ind_No).howCome = '  Random  ';
