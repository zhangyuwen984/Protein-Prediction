function MakeupSurface(whichInd,lattice,coor,type, numIons)

%This routine is used to compensate the atom 
%if some atoms are removed during relaxation

global ORG_STRUC
global POP_STRUC 

Step         = POP_STRUC.POPULATION(whichInd).Step;
cell         = POP_STRUC.POPULATION(whichInd).cell;
numatom      = ORG_STRUC.numIons*prod(cell) - numIons; 
vacuum       = ORG_STRUC.vacuumSize(Step);
atomType     = ORG_STRUC.atomType;

%Delete
NDelete = numatom;
NDelete(find(NDelete>0)) = 0;
NDelete = -1*NDelete;

%Add
NAdd    = numatom;
NAdd(find(NAdd<0)) = 0;

disp(['Delete atoms: ' num2str(NDelete) ' Add atom: ' num2str(NAdd) ]);
if sum(NDelete) > 0
   [coor, type] = Delete_Atom(coor, type, NDelete, atomType);
end
newCoords    = Add_Atom(lattice, coor, type, NAdd, atomType);

surnumIons = ORG_STRUC.numIons*prod(cell);

bulk_lattice = POP_STRUC.POPULATION(whichInd).Bulk_LATTICE;
bulk_pos     = POP_STRUC.POPULATION(whichInd).Bulk_COORDINATES;
bulk_atyp    = POP_STRUC.POPULATION(whichInd).Bulk_typesAList;
bulk_numIons = POP_STRUC.POPULATION(whichInd).Bulk_numIons;          


[lat,candidate,numIons,typesAList,chanAList] = makeSurface...
(lattice,newCoords,surnumIons,bulk_lattice,bulk_pos,bulk_atyp,bulk_numIons, vacuum);

POP_STRUC.POPULATION(whichInd).Surface_COORDINATES = newCoords;
POP_STRUC.POPULATION(whichInd).Surface_LATTICE     = lattice;
POP_STRUC.POPULATION(whichInd).Surface_numIons     = surnumIons;
POP_STRUC.POPULATION(whichInd).COORDINATES         = candidate;
POP_STRUC.POPULATION(whichInd).numIons             = numIons;
POP_STRUC.POPULATION(whichInd).LATTICE             = lat;
POP_STRUC.POPULATION(whichInd).typesAList          = typesAList;
POP_STRUC.POPULATION(whichInd).chanAList           =chanAList;
