function WriteIndividualOutput_201(Ind_No)

% USPEX 9.4.0
% We choose to update the files requiring I/O files first
% In case sometimes I/O errors in file system
% Lastly updated by Qiang Zhu (2013/11/22)

global POP_STRUC
global ORG_STRUC

atomType  = ORG_STRUC.atomType;
resFolder = POP_STRUC.resFolder;

POP_STRUC.POPULATION(Ind_No).symg = 0;
symg    = POP_STRUC.POPULATION(Ind_No).symg;
count       =POP_STRUC.POPULATION(Ind_No).Number;
INIT_numIons=POP_STRUC.POPULATION(Ind_No).INIT_numIons;
INIT_LAT    =POP_STRUC.POPULATION(Ind_No).INIT_LAT;
INIT_COORD  =POP_STRUC.POPULATION(Ind_No).INIT_COORD;
Write_POSCAR(atomType, count, symg, INIT_numIons, INIT_LAT, INIT_COORD);
[nothing, nothing] = unix([' cat POSCAR      >> ' POP_STRUC.resFolder '/gatheredPOSCARS_unrelaxed']);

lattice = POP_STRUC.POPULATION(Ind_No).LATTICE;
coor    = POP_STRUC.POPULATION(Ind_No).COORDINATES;
numIons = POP_STRUC.POPULATION(Ind_No).numIons;

count1 = 0;
order = zeros(sum(numIons),1);
Surface_order   = POP_STRUC.POPULATION(Ind_No).Surface_order;
for i=1:size(coor,1)
    if POP_STRUC.POPULATION(Ind_No).chanAList(i)==1
       count1=count1+1;
       order(i) = Surface_order(count1);
    end
end
POP_STRUC.POPULATION(Ind_No).order = order;

Write_POSCAR(atomType, count, symg, numIons, lattice, coor);
Write_POSCAR_order(atomType, count, symg, numIons, lattice, coor, order);
[nothing, nothing] = unix([' cat POSCAR      >> ' POP_STRUC.resFolder '/gatheredPOSCARS']);
[nothing, nothing] = unix([' cat POSCAR_order >>' POP_STRUC.resFolder '/gatheredPOSCARS_order']);

update_USPEX_INDIVIDUAL(POP_STRUC.POPULATION(Ind_No), resFolder, ...
                        POP_STRUC.generation, atomType);
WriteOUTPUT(count, resFolder);
