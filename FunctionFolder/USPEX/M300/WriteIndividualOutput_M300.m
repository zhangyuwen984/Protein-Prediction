function WriteIndividualOutput_M300(Ind_No)

% USPEX 9.4.0
% We choose to update the files requiring I/O files first
% In case sometimes I/O errors in file system
% Lastly updated by Qiang Zhu (2013/11/22)

global POP_STRUC
global ORG_STRUC

resFolder = POP_STRUC.resFolder;
atomType = ORG_STRUC.atomType;
Count   = POP_STRUC.POPULATION(Ind_No).Number;
symg    = POP_STRUC.POPULATION(Ind_No).symg;
lattice = POP_STRUC.POPULATION(Ind_No).LATTICE;
coor    = POP_STRUC.POPULATION(Ind_No).COORDINATES;
numIons = POP_STRUC.POPULATION(Ind_No).numIons;
order   = POP_STRUC.POPULATION(Ind_No).order;

Write_POSCAR(atomType, Count, symg, numIons, lattice, coor);

[a,b]=unix([' cat POSCAR      >> ' resFolder '/gatheredPOSCARS']);
%[a,b]=unix([' cat POSCAR_order >>' resFolder '/gatheredPOSCARS_order']);

INIT_numIons=POP_STRUC.POPULATION(Ind_No).INIT_numIons;
INIT_LAT    =POP_STRUC.POPULATION(Ind_No).INIT_LAT;
INIT_COORD  =POP_STRUC.POPULATION(Ind_No).INIT_COORD;
Write_POSCAR(atomType, Count, symg, INIT_numIons, INIT_LAT, INIT_COORD);
[a,b]=unix([' cat POSCAR      >> ' resFolder '/gatheredPOSCARS_unrelaxed']);

update_USPEX_INDIVIDUAL(POP_STRUC.POPULATION(Ind_No), resFolder, ...
                        POP_STRUC.generation, ORG_STRUC.atomType);
WriteOUTPUT(Count, resFolder);
