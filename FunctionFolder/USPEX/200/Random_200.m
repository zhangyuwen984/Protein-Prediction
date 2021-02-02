function Random_200(Ind_No)

global ORG_STRUC
global OFF_STRUC

CellList = findcell(ORG_STRUC.reconstruct);
Vacuum = ORG_STRUC.vacuumSize(1);
ID   = ceil(rand()*size(CellList,1));
cell =CellList(ID,:);
ANS  = Random_Init_200(Ind_No, Vacuum, cell);

OFF_STRUC.POPULATION(Ind_No).COORDINATES         = ANS.candidate;
OFF_STRUC.POPULATION(Ind_No).numIons             = ANS.numIons;
OFF_STRUC.POPULATION(Ind_No).LATTICE             = ANS.lat;
OFF_STRUC.POPULATION(Ind_No).typesAList          = ANS.typesAList;
OFF_STRUC.POPULATION(Ind_No).chanAList           = ANS.chanAList;
OFF_STRUC.POPULATION(Ind_No).Surface_LATTICE     = ANS.sur_lat;
OFF_STRUC.POPULATION(Ind_No).Surface_COORDINATES = ANS.sur_candidate;
OFF_STRUC.POPULATION(Ind_No).Surface_numIons     = ANS.sur_numIons;
OFF_STRUC.POPULATION(Ind_No).Bulk_LATTICE        = ANS.bulk_lat;
OFF_STRUC.POPULATION(Ind_No).Bulk_COORDINATES    = ANS.bulk_pos;
OFF_STRUC.POPULATION(Ind_No).Bulk_typesAList     = ANS.bulk_atyp;
OFF_STRUC.POPULATION(Ind_No).Bulk_numIons        = ANS.bulk_numIons;

OFF_STRUC.POPULATION(Ind_No).cell = cell;
OFF_STRUC.POPULATION(Ind_No).howCome = '  Random  ';
